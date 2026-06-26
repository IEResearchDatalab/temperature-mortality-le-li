################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 05_Disaggregate.R
#   PCLM disaggregation of ANs from wide age groups to single years.
#
# For each city × year × GCM × SSP × temperature range, the 5 grouped
# ANs (20-44, 45-64, 65-74, 75-84, 85+) are disaggregated to single
# ages 20-100 using PCLM (Penalized Composite Link Model).
#
# Changes from the original 03_Disaggregate.R:
#   - age groups are explicitly ordered before PCLM
#   - grouped totals are checked before and after disaggregation
#   - missing age groups trigger a warning instead of silent skip
#   - uncertainty low/high bounds are also disaggregated (approximate)
#   - open-age handling is explicit
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")

library(ungroup)
library(foreach)
library(doParallel)

#---------------------------
# Load annual ANs
#---------------------------

annual_ans <- fread(file.path(path_out, "ans_annual_all_cities.csv"))

# Filter to only the cities we have schedules for (handles demo mode)
city_mx <- readRDS(file.path(path_out, "city_mortality_schedules.rds"))
city_codes_available <- unique(city_mx$city_code)
annual_ans <- annual_ans[city_code %in% city_codes_available]

city_gcm_ssp <- unique(annual_ans[, .(city_code, gcm, ssp)])
n_combos <- nrow(city_gcm_ssp)

cat(sprintf("Disaggregating %d city-GCM-SSP combinations...\n", n_combos))

#---------------------------
# PCLM function with explicit age ordering and validation
#---------------------------

disaggregate_an <- function(an_values, age_starts, nlast = 16,
  age_labs_ordered = c("20-44", "45-64", "65-74", "75-84", "85+"),
  kr = 2, deg = 3) {

  if (any(is.na(an_values)) || any(an_values < 0)) {
    return(NULL)
  }
  if (sum(an_values) == 0) {
    return(rep(0, 81))
  }

  M <- pclm(x = age_starts, y = an_values, nlast = nlast,
    out.step = 1, control = list(lambda = NA, kr = kr, deg = deg,
      max.iter = 200, tol = 1e-8, opt.method = "BIC"))
  result <- as.numeric(fitted(M))
  result[result < 0] <- 0

  # Check that grouped total is preserved
  total_before <- sum(an_values)
  total_after <- sum(result)
  if (abs(total_after - total_before) / max(total_before, 1) > 0.05) {
    warning(sprintf("PCLM total mismatch: %.2f vs %.2f", total_before, total_after))
  }

  result
}

#---------------------------
# Age group definitions (explicitly ordered)
#---------------------------

age_starts <- c(20, 45, 65, 75, 85)
age_labs_ordered <- c("20-44", "45-64", "65-74", "75-84", "85+")
single_age_range <- 20:100
n_single <- length(single_age_range)

#---------------------------
# Setup parallel backend
#---------------------------

cl <- makeCluster(ncores)
registerDoParallel(cl)

dir.create("temp", showWarnings = FALSE)
writeLines(c(""), "temp/log_pclm.txt")
cat(as.character(as.POSIXct(Sys.time())), file = "temp/log_pclm.txt", append = TRUE)

#---------------------------
# Parallel loop
#---------------------------

all_single <- foreach(i = seq_len(n_combos),
  .packages = c("data.table", "ungroup"),
  .export  = c("disaggregate_an", "age_starts", "age_labs_ordered",
               "single_age_range", "n_single", "annual_ans"),
  .combine = function(x, y) rbind(x, y, fill = TRUE),
  .errorhandling = "pass") %dopar% {

  cc <- city_gcm_ssp$city_code[i]
  gc <- city_gcm_ssp$gcm[i]
  ss <- city_gcm_ssp$ssp[i]

  if (i %% 50 == 0)
    cat("\n", "i = ", i, cc, gc, as.character(Sys.time()), "\n",
      file = "temp/log_pclm.txt", append = TRUE)

  city_data <- annual_ans[city_code == cc & gcm == gc & ssp == ss]

  results_list <- list()
  combos <- unique(city_data[, .(year, temp_range)])

  for (j in seq_len(nrow(combos))) {
    yr <- combos$year[j]
    tr <- combos$temp_range[j]

    subset_dt <- city_data[year == yr & temp_range == tr]

    # Match to ordered age groups
    subset_dt <- subset_dt[match(age_labs_ordered, age_group)]
    an_vec <- subset_dt$an_est
    an_low_vec <- subset_dt$an_low
    an_hi_vec <- subset_dt$an_hi

    # Check for missing age groups
    missing_groups <- age_labs_ordered[is.na(an_vec)]
    if (length(missing_groups) > 0) {
      warning(sprintf("%s %s %s %s yr %s: missing age groups %s",
        cc, gc, ss, tr, yr, paste(missing_groups, collapse = ",")))
      next
    }

    # Disaggregate point estimate
    single_est <- disaggregate_an(an_vec, age_starts)
    if (is.null(single_est)) next

    # Disaggregate confidence bounds (approximate: same PCLM weights applied to bounds)
    if (any(!is.na(an_low_vec))) {
      an_low_vec[is.na(an_low_vec)] <- 0
      weights <- single_est / sum(single_est) * sum(an_vec)
      single_low <- an_low_vec[1] * (weights / an_vec[1])
      single_low[is.nan(single_low) | is.infinite(single_low)] <- 0
      single_hi <- an_hi_vec[1] * (weights / an_vec[1])
      single_hi[is.nan(single_hi) | is.infinite(single_hi)] <- 0
    } else {
      single_low <- rep(NA_real_, n_single)
      single_hi <- rep(NA_real_, n_single)
    }

    results_list[[length(results_list) + 1]] <- data.table(
      city_code = cc, gcm = gc, ssp = ss,
      year = yr, age = single_age_range, temp_range = tr,
      AN = single_est,
      AN_low = single_low,
      AN_hi = single_hi)
  }

  if (length(results_list) > 0)
    rbindlist(results_list)
  else
    NULL
}

stopCluster(cl)

#---------------------------
# Save
#---------------------------

fwrite(all_single, file.path(path_out, "ans_single_age_all_cities.csv"))
cat(sprintf("Single-age AN: %s rows\n",
  format(nrow(all_single), big.mark = ",")))