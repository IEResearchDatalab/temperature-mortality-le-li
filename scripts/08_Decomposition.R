################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 08_Decomposition.R
#   Decompose annual changes in e65 and lifespan inequality (SD above 65)
#   by age and temperature range using the linear integral (Horiuchi-style)
#   decomposition.
#
# Method:
#   For each consecutive year-pair (t, t+1):
#   1. Build midpoint mortality schedule (average of t and t+1)
#   2. Compute total change in e65 and SD
#   3. For each age x ≥ 65 and each cause c, compute:
#        Contribution_c(x) = (∂I/∂μ_c(x)) * Δμ_c(x)
#      where ∂I/∂μ_c(x) is the numerical sensitivity evaluated at
#      the midpoint schedule, and Δμ_c(x) is the change in cause-
#      specific mortality from t to t+1.
#   4. Sum annual contributions into reporting periods (2030, 2050, 2090).
#
# Input: analysis_dataset_all_cities.csv
# Output: decomposition_all_cities.csv
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")

#---------------------------
# Load analysis dataset
#---------------------------

analysis_all <- fread(file.path(path_out, "analysis_dataset_all_cities.csv"))

#---------------------------
# Cause labels
#---------------------------

cause_cols <- c("residual_deaths", "extreme_cold", "moderate_cold",
  "moderate_heat", "extreme_heat")
cause_labels <- c("residual", "extreme_cold", "moderate_cold",
  "moderate_heat", "extreme_heat")

#---------------------------
# Life-table function (65+, returns e65 and SD)
#---------------------------

compute_e65_sd <- function(mx_vec, age_vec, ax = 0.5) {
  keep <- age_vec >= 65
  mx <- mx_vec[keep]
  age <- age_vec[keep]
  ord <- order(age)
  mx <- mx[ord]
  age <- age[ord]
  n <- c(diff(age), Inf)

  # qx
  qx <- (n * mx) / (1 + n * (1 - ax) * mx)
  qx[length(qx)] <- 1

  # lx, dx
  lx <- numeric(length(qx))
  lx[1] <- 100000
  for (i in 2:length(qx)) lx[i] <- lx[i - 1] * (1 - qx[i - 1])
  dx <- lx * qx

  # Lx, Tx, ex
  Lx <- numeric(length(qx))
  for (i in 1:(length(qx) - 1)) Lx[i] <- n[i] * (lx[i] - (1 - ax) * dx[i])
  Lx[length(qx)] <- if (mx[length(qx)] > 0) lx[length(qx)] / mx[length(qx)] else 0
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  ex[is.nan(ex)] <- 0

  e65_val <- ex[age == 65]

  # SD of age at death above 65
  lx_start <- lx[1]
  mean_age <- sum(dx * (age + ax)) / lx_start
  sd_val <- sqrt(sum(dx * (age + ax - mean_age)^2) / lx_start)

  list(e65 = e65_val, sd = sd_val)
}

#---------------------------
# Decomposition function for one year-pair
#---------------------------

decompose_year_pair <- function(data_t1, data_t2, eps_rel = 1e-4) {
  # data_t1, data_t2: data.tables with columns age, population, and cause_cols
  # Returns: data.table with age, cause, delta_e65, delta_sd

  ages <- data_t1$age
  age65_idx <- which(ages >= 65)
  ages_65plus <- ages[age65_idx]

  # Cause-specific mx at t1 and t2
  mx_t1 <- do.call(cbind, lapply(cause_cols, function(col) data_t1[[col]] / data_t1$population))
  mx_t2 <- do.call(cbind, lapply(cause_cols, function(col) data_t2[[col]] / data_t2$population))
  colnames(mx_t1) <- cause_labels
  colnames(mx_t2) <- cause_labels

  # Total mx
  mx_total_t1 <- rowSums(mx_t1)
  mx_total_t2 <- rowSums(mx_t2)

  # Midpoint schedule
  mx_mid <- (mx_total_t1 + mx_total_t2) / 2
  mx_cause_mid <- (mx_t1 + mx_t2) / 2

  # Baseline e65 and SD at midpoint
  baseline <- compute_e65_sd(mx_mid, ages)

  # Delta mx per cause
  delta_mx <- mx_t2 - mx_t1

  # Results
  result_list <- list()

  for (i in seq_along(age65_idx)) {
    ai <- age65_idx[i]
    age_val <- ages[ai]
    mx_total_at_age <- mx_mid[ai]

    for (c_idx in seq_along(cause_labels)) {
      delta <- delta_mx[ai, c_idx]
      if (abs(delta) < 1e-15) next

      # Numerical sensitivity: perturb total mx at this age by eps
      eps <- max(mx_total_at_age * eps_rel, 1e-10)
      mx_perturbed <- mx_mid
      mx_perturbed[ai] <- mx_perturbed[ai] + eps
      perturbed <- compute_e65_sd(mx_perturbed, ages)

      # Sensitivity = (I_perturbed - I) / eps  [note: only cause c at age x changes]
      # But the perturbation affects total mx at age x, which is the sum of all causes.
      # The marginal effect of changing cause c specifically at age x is:
      #   dI/dμ_c(x) = dI/dμ_total(x) * dμ_total(x)/dμ_c(x) = dI/dμ(x) * 1
      # So we perturb total mx and attribute the sensitivity to each cause.

      sens_e65 <- (perturbed$e65 - baseline$e65) / eps
      sens_sd <- (perturbed$sd - baseline$sd) / eps

      result_list[[length(result_list) + 1]] <- data.table(
        age = age_val,
        cause = cause_labels[c_idx],
        delta_e65 = sens_e65 * delta,
        delta_sd = sens_sd * delta)
    }
  }

  rbindlist(result_list)
}

#---------------------------
# Run decomposition per city-GCM-SSP
#---------------------------

combos <- unique(analysis_all[, .(city_code, gcm, ssp)])
n_combos <- nrow(combos)

cat(sprintf("Decomposition for %d city-GCM-SSP combos...\n", n_combos))

all_decomp <- list()

for (i in seq_len(n_combos)) {
  cc <- combos$city_code[i]
  gc <- combos$gcm[i]
  ss <- combos$ssp[i]

  if (i %% 50 == 0)
    cat(sprintf("  %d / %d: %s %s %s\n", i, n_combos, cc, gc, ss))

  series <- analysis_all[city_code == cc & gcm == gc & ssp == ss]
  setorder(series, year, age)

  years <- sort(unique(series$year))

  # Process consecutive year-pairs
  for (j in seq_len(length(years) - 1)) {
    yr_t1 <- years[j]
    yr_t2 <- years[j + 1]

    data_t1 <- series[year == yr_t1]
    data_t2 <- series[year == yr_t2]

    if (nrow(data_t1) == 0 || nrow(data_t2) == 0) next

    dc <- decompose_year_pair(data_t1, data_t2)

    # Total change in e65 and SD (for validation)
    mx_t1_total <- data_t1$deaths / data_t1$population
    mx_t2_total <- data_t2$deaths / data_t2$population
    lt_t1 <- compute_e65_sd(mx_t1_total, data_t1$age)
    lt_t2 <- compute_e65_sd(mx_t2_total, data_t2$age)
    total_delta_e65 <- lt_t2$e65 - lt_t1$e65
    total_delta_sd <- lt_t2$sd - lt_t1$sd

    dc[, `:=`(
      city_code = cc,
      gcm = gc,
      ssp = ss,
      year_t1 = yr_t1,
      year_t2 = yr_t2,
      total_delta_e65 = total_delta_e65,
      total_delta_sd = total_delta_sd
    )]

    all_decomp[[length(all_decomp) + 1]] <- dc
  }
}

decomp_all <- rbindlist(all_decomp, fill = TRUE)

#---------------------------
# Aggregate annual contributions to reporting periods
#---------------------------

# Define reporting periods
decomp_all[, period := fcase(
  year_t1 >= 2020 & year_t1 <= 2029, "2030",
  year_t1 >= 2040 & year_t1 <= 2049, "2050",
  year_t1 >= 2080 & year_t1 <= 2089, "2090",
  default = NA_character_
)]
decomp_all <- decomp_all[!is.na(period)]

# Sum annual contributions within period
decomp_period <- decomp_all[,
  .(delta_e65 = sum(delta_e65),
    delta_sd = sum(delta_sd),
    n_years = .N / .GRP),
  by = .(city_code, gcm, ssp, period, age, cause)]

#---------------------------
# Save
#---------------------------

fwrite(decomp_all, file.path(path_out, "decomposition_annual_all_cities.csv"))
fwrite(decomp_period, file.path(path_out, "decomposition_period_all_cities.csv"))

cat(sprintf("Annual decomposition: %s rows\n", format(nrow(decomp_all), big.mark = ",")))
cat(sprintf("Period decomposition: %s rows\n", format(nrow(decomp_period), big.mark = ",")))