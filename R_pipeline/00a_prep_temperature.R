#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality / life expectancy pipeline
#
# R Pipeline Step 00a: Observed temperature series and per-city thresholds
#   Builds data/prep_data.RData, consumed by 01_attribution.R:
#     - obs_data:   ERA5-Land daily mean temperature per city (the series used
#                   by Masselot et al. 2023 to estimate the ERFs), 1990-2019
#     - thresholds: city_results.csv (one row per city x age group) plus the
#                   2.5th / 97.5th percentiles of obs_data that delimit the
#                   extreme temperature ranges (Lloyd et al. 2024 definitions)
#     - cities:     vector of URAU codes
#   Inputs are read from data/ only (see README, "Data").
#
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

message("\n[00a] Building observed temperature series and thresholds...")

#----- Parameters

# Reference period for the extreme-range percentiles. 1990-2019 is the full
# ERA5-Land series used to estimate the ERFs (Masselot et al. 2023); it is the
# same window that defines the ERF knots and boundaries in 01_attribution.R.
threshold_years <- c(1990L, 2019L)
extreme_probs <- c(p2_5 = 0.025, p97_5 = 0.975)

out_file <- "data/prep_data.RData"
check_dir <- "results/checks"
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
checks_file <- file.path(check_dir, "00a_prep_temperature_checks.csv")

#----- Load inputs

era5 <- as.data.table(read_parquet("data/era5series.gz.parquet"))
city_results <- fread("data/city_results.csv")

obs_data <- era5[, .(URAU_CODE, date = as.IDate(date), tmean_obs = as.numeric(era5landtmean))]
setorder(obs_data, URAU_CODE, date)

#----- Per-city extreme-range percentiles

pct <- obs_data[year(date) %between% threshold_years, .(
  p2_5 = as.numeric(quantile(tmean_obs, extreme_probs[["p2_5"]], na.rm = TRUE)),
  p97_5 = as.numeric(quantile(tmean_obs, extreme_probs[["p97_5"]], na.rm = TRUE))
), by = URAU_CODE]

thresholds <- merge(city_results, pct, by = "URAU_CODE", all.x = TRUE, sort = TRUE)
setkey(thresholds, URAU_CODE)
cities <- sort(unique(city_results$URAU_CODE))

#----- Invariant checks

n_days <- obs_data[, .N, by = URAU_CODE]
checks <- data.table(
  check_name = c("all_cities_have_obs", "obs_finite", "thresholds_complete", "p2_5_below_p97_5", "obs_period"),
  status = c(
    if (setequal(cities, unique(obs_data$URAU_CODE))) "PASS" else "FAIL",
    if (all(is.finite(obs_data$tmean_obs))) "PASS" else "FAIL",
    if (!anyNA(thresholds[, .(p2_5, p97_5, mmt)])) "PASS" else "FAIL",
    if (all(thresholds$p2_5 < thresholds$p97_5)) "PASS" else "FAIL",
    if (all(range(year(obs_data$date)) == threshold_years)) "PASS" else "FAIL"
  ),
  value = c(
    sprintf("%d cities", length(cities)),
    sprintf("%d rows", nrow(obs_data)),
    sprintf("%d rows", nrow(thresholds)),
    sprintf("min gap = %.2f C", min(thresholds$p97_5 - thresholds$p2_5)),
    paste(range(year(obs_data$date)), collapse = "-")
  )
)
fwrite(checks, checks_file)
if (any(checks$status == "FAIL")) {
  stop(sprintf("00a_prep_temperature.R failed %d check(s); see %s", sum(checks$status == "FAIL"), checks_file), call. = FALSE)
}

save(obs_data, thresholds, cities, file = out_file)
message("Saved ", out_file, " (", length(cities), " cities, ", nrow(obs_data), " daily rows)")
