################################################################################
# Script: 07_le_li_decomposition.R
# Description: Disaggregate ANs and baseline deaths using pclm (ungroup)
#              and estimate changes in LE and LI (Aburto/Lloyd methodology)
################################################################################

library(data.table)
library(dplyr)
library(ungroup)
library(DemoDecomp)
library(parallel)

# Create output directory
dir.create("results/le_li_decomposition", recursive = TRUE, showWarnings = FALSE)

# Load data
ans_data <- fread("results/le_li_input/le_li_input_ans.csv")

# We have age groups: 65-74, 75-84, 85+ 
# Needs to map to age intervals for pclm
# 65-74: [65, 75) -> width 10
# 75-84: [75, 85) -> width 10
# 85+:   [85, 100+] -> width 15 (assuming 100 as open interval end for pclm)

ans_data[agegroup == "65-74", `:=`(age_start = 65, age_width = 10)]
ans_data[agegroup == "75-84", `:=`(age_start = 75, age_width = 10)]
ans_data[agegroup == "85+",    `:=`(age_start = 85, age_width = 15)] # Open interval

# Groups that uniquely identify a "population" to disaggregate
# We need to disaggregate:
# 1. Total deaths (baseline)
# 2. AN for each range
# 3. Population (to get rates)

# Let's first disaggregate baseline deaths and population per (cntr_name, decade, ssp, gcm)
# Note: population and deaths_baseline are common across ranges in the input file.

unique_pops <- unique(ans_data[, .(cntr_name, decade, ssp, gcm)])

p_disaggregate <- function(row_idx) {
  pop_info <- unique_pops[row_idx]
  sub_data <- ans_data[cntr_name == pop_info$cntr_name & 
                       decade == pop_info$decade & 
                       ssp == pop_info$ssp & 
                       gcm == pop_info$gcm]
  
  # Prepare for pclm
  # We need unique age groups (65-74, 75-84, 85+)
  age_agg <- unique(sub_data[, .(agegroup, age_start, age_width, pop, death_baseline)])
  setorder(age_agg, age_start)
  
  x <- age_agg$age_start
  n <- age_agg$age_width
  
  # Disaggregate Population
  pclm_pop <- tryCatch({
    pclm(x = x, y = age_agg$pop, nlast = 15)$y
  }, error = function(e) return(rep(NA, 35))) 
  
  # Disaggregate Baseline Deaths
  pclm_deaths <- tryCatch({
    pclm(x = x, y = age_agg$death_baseline, nlast = 15)$y
  }, error = function(e) return(rep(NA, 35)))
  
  # Create a base table for single ages 65-99
  base_single <- data.table(
    cntr_name = pop_info$cntr_name,
    decade = pop_info$decade,
    ssp = pop_info$ssp,
    gcm = pop_info$gcm,
    age = 65:99, # 35 intervals
    pop = as.vector(pclm_pop),
    death_baseline = as.vector(pclm_deaths)
  )
  
  # Now disaggregate AN for each range
  ranges <- unique(sub_data$range)
  range_results <- lapply(ranges, function(r) {
    an_agg <- sub_data[range == r]
    setorder(an_agg, age_start)
    pclm_an <- tryCatch({
      pclm(x = x, y = an_agg$an, nlast = 15)$y
    }, error = function(e) return(rep(0, 35)))
    
    dt <- copy(base_single)
    dt[, range := r]
    dt[, an := as.vector(pclm_an)]
    return(dt)
  })
  
  return(rbindlist(range_results))
}

# Run disaggregation in parallel
num_cores <- min(detectCores(), 8)
cat("Running disaggregation on", num_cores, "cores...\n")
results_list <- mclapply(1:nrow(unique_pops), p_disaggregate, mc.cores = num_cores)

# Check for errors in results_list
errors <- sapply(results_list, inherits, "try-error")
if (any(errors)) {
    cat("Error in some workers:\n")
    print(results_list[errors][1])
}

ans_single_age <- rbindlist(results_list)

cat("Disaggregation finished. Rows in ans_single_age:", nrow(ans_single_age), "\n")
if (nrow(ans_single_age) > 0) {
    cat("Columns in ans_single_age:", paste(names(ans_single_age), collapse=", "), "\n")
    cat("Unique ranges in ans_single_age:", paste(unique(ans_single_age$range), collapse=", "), "\n")
    print(head(ans_single_age))
    print(any(is.na(ans_single_age$pop)))
}

# Verify results and filter NAs
# Now it should have values correctly
ans_single_age <- ans_single_age[!is.na(pop)]

# Calculate mx
ans_single_age[, mx_total := death_baseline / pop]
ans_single_age[, mx_range := an / pop]

# Reshape to wide range for decomposition
# We need causes: ExtrCold, ModCold, ModHeat, ExtrHeat, and "rest"
# Fill missing ranges with 0 (some future scenarios might not have ExtrCold)
ans_wide <- dcast(ans_single_age, cntr_name + decade + ssp + gcm + age + pop + death_baseline + mx_total ~ range, value.var = "mx_range", fill = 0)

cat("Wide data rows:", nrow(ans_wide), "\n")

# Ensure all 4 columns exist even if fill=0 didn't create them
for (col in c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")) {
  if (!(col %in% names(ans_wide))) {
    ans_wide[, (col) := 0]
  }
}

# mx_rest = mx_total - sum(all temperature mx)
ans_wide[, rest := mx_total - (ExtrCold + ModCold + ModHeat + ExtrHeat)]
# Ensure no negative mortality
ans_wide[rest < 0, rest := 0]

# Melt back to cause format compatible with Lloyd's code
ans_long <- melt(ans_wide, 
                 id.vars = c("cntr_name", "decade", "ssp", "gcm", "age", "pop", "mx_total"),
                 measure.vars = intersect(c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat", "rest"), names(ans_wide)),
                 variable.name = "cause", value.name = "mx_cause")

cat("Long data rows:", nrow(ans_long), "\n")
cat("Long data unique causes:", paste(unique(ans_long$cause), collapse=", "), "\n")

# ------------------------------------------------------------------------------
# Decomposition Logic (Aburto/Lloyd)
# ------------------------------------------------------------------------------

# Life expectancy function for 65+
life_expectancy_from_mx_65plus <- function(mx, nx = rep(1, length(mx))) {
    # Close at age 100+ -> Inf
    # We have ages 65-99, so 35 ages.
    nx[length(nx)] <- 1/mx[length(nx)] # Crude way to close but robust for 99
    px <- exp(-mx * nx)
    lx <- head(cumprod(c(1, px)), -1)
    # Simple midpoint for single age
    Lx <- lx * exp(-mx * 0.5) 
    Tx <- rev(cumsum(rev(Lx)))
    ex <- Tx / lx
    return(ex[1]) # Return e65
}

# Wrapper for horiuchi
le_wrapper <- function(mx_vec, age_vec) {
    # horiuchi passes causes then ages: mx_c1_a65, mx_c1_a66... mx_c2_a65...
    n_age <- length(age_vec)
    n_cause <- length(mx_vec) / n_age
    mx_mat <- matrix(mx_vec, nrow = n_age, ncol = n_cause, byrow = FALSE)
    mx_total <- rowSums(mx_mat)
    life_expectancy_from_mx_65plus(mx_total)
}

# ------------------------------------------------------------------------------
# Run decomposition loop
# ------------------------------------------------------------------------------
cat("Running LE decomposition...\n")

# Need to compare 2010 vs 2090 for each (cntr_name, ssp, gcm)
comparison_groups <- unique(ans_long[, .(cntr_name, ssp, gcm)])

process_decomposition <- function(grp_idx) {
    grp <- comparison_groups[grp_idx]
    
    # 2010 Data
    d1 <- ans_long[cntr_name == grp$cntr_name & ssp == grp$ssp & gcm == grp$gcm & decade == 2010]
    setorder(d1, cause, age)
    
    # 2090 Data
    d2 <- ans_long[cntr_name == grp$cntr_name & ssp == grp$ssp & gcm == grp$gcm & decade == 2090]
    setorder(d2, cause, age)
    
    if (nrow(d1) == 0 || nrow(d2) == 0) return(NULL)
    
    # Check that ages and causes match
    if (!all(unique(d1$age) == unique(d2$age)) || !all(unique(d1$cause) == unique(d2$cause))) return(NULL)
    
    # Rates vectors
    r1 <- d1$mx_cause
    r2 <- d2$mx_cause
    
    age_vec <- unique(d1$age)
    
    # Total change in e65
    e65_1 <- le_wrapper(r1, age_vec)
    e65_2 <- le_wrapper(r2, age_vec)
    diff_total <- e65_2 - e65_1
    
    # Horiuchi decomposition
    decomp <- tryCatch({
        horiuchi(le_wrapper, r1, r2, N = 20, age_vec = age_vec)
    }, error = function(e) return(rep(NA, length(r1))))
    
    # Aggregate contributions by cause
    res_dt <- data.table(
        age = d1$age,
        cause = d1$cause,
        contribution = decomp
    )
    
    cause_contrib <- res_dt[, .(contribution = sum(contribution)), by = cause]
    
    # Add metadata
    cause_contrib[, `:=`(
        cntr_name = grp$cntr_name,
        ssp = grp$ssp,
        gcm = grp$gcm,
        e65_2010 = e65_1,
        e65_2090 = e65_2,
        diff_total = diff_total
    )]
    
    return(cause_contrib)
}

results_decomp <- mclapply(1:nrow(comparison_groups), process_decomposition, mc.cores = num_cores)
le_decomp_final <- rbindlist(results_decomp)

# Save results
fwrite(le_decomp_final, "results/le_li_decomposition/le_decomposition_results.csv")

cat("LE decomposition complete. Results saved to results/le_li_decomposition/le_decomposition_results.csv\n")
    causes <- unique(d2020$cause)
    
    # Horiuchi decomposition
    decomp_results <- horiuchi(func = le_wrapper, pars1 = pars1, pars2 = pars2, N = 20, age_vec = age_vec)
    
    # Reshape results
    res_dt <- data.table(
        cntr_name = scen$cntr_name,
        ssp = scen$ssp,
        gcm = scen$gcm,
        age = rep(age_vec, length(causes)),
        cause = rep(causes, each = length(age_vec)),
        contribution = decomp_results
    )
    
    return(res_dt)
}

le_results_list <- mclapply(1:nrow(unique_scenarios), decomp_worker, mc.cores = num_cores)
le_decomp_final <- rbindlist(le_results_list)

# Save results
fwrite(le_decomp_final, "data/le_decomposition_results.csv")

cat("LE decomposition complete. Results saved to data/le_decomposition_results.csv\n")
