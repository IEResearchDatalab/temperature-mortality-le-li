#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# R Pipeline Part 00a: Observed temperature series and per-city thresholds
#   Builds data/prep_data.RData, used by 01_attribution.R:
#     - obs_data:   ERA5-Land daily mean temperature per city, 1990-2019 (the
#                   series used by Masselot et al. 2023 to estimate the ERFs)
#     - thresholds: city_results.csv (one row per city x age group) with
#                   - p2_5 / p97_5: percentiles delimiting extreme cold / heat
#                     (Lloyd et al. 2024), over `threshold_years`
#                   - mmt: Masselot et al. (2025) rule, argmin of the ERF over
#                     percentiles 25-99 (the 2023 value is kept as `mmt_2023`)
#     - cities:     URAU codes
#
################################################################################

source("R_pipeline/00_pkg_params.R")

message("\n[00a] Building observed temperature series and thresholds...")

out_file <- "data/prep_data.RData"
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

#----- MMT, Masselot et al. (2025) 03_attribution.R:
#   tper  = percentiles of the full observed series on the `predper` grid
#   basis = bs(degree 2, knots at 10/75/90th percentiles, Bound = range(tper))
#   mmt   = tper value minimising the ERF among percentiles 25-99

coefs <- fread("data/coefs.csv")
mmt_tbl <- obs_data[, {
  tper <- quantile(tmean_obs, predper / 100)
  argvar <- list(fun = varfun, degree = vardegree, knots = tper[paste0(varper, ".0%")], Bound = range(tper))
  bper <- suppressWarnings(do.call(onebasis, c(list(x = tper), argvar)))
  ind <- tper >= tper["25.0%"] & tper <= tper["99.0%"]
  cc <- coefs[URAU_CODE == .BY$URAU_CODE]
  list(agegroup = cc$agegroup,
       mmt = vapply(seq_len(nrow(cc)), function(i) {
         b <- as.numeric(cc[i, .(b1, b2, b3, b4, b5)])
         as.numeric(tper[ind][which.min(drop(bper[ind, ] %*% b))])
       }, numeric(1)))
}, by = URAU_CODE]

thresholds <- merge(city_results, pct, by = "URAU_CODE", all.x = TRUE, sort = TRUE)
setnames(thresholds, "mmt", "mmt_2023")
thresholds <- merge(thresholds, mmt_tbl, by = c("URAU_CODE", "agegroup"), all.x = TRUE, sort = TRUE)
setkey(thresholds, URAU_CODE)
cities <- sort(unique(city_results$URAU_CODE))

#----- Invariant checks

checks <- data.table(
  check_name = c("all_cities_have_obs", "obs_finite", "thresholds_complete", "p2_5_below_p97_5", "obs_period", "mmt_within_25_99"),
  status = c(
    if (setequal(cities, unique(obs_data$URAU_CODE))) "PASS" else "FAIL",
    if (all(is.finite(obs_data$tmean_obs))) "PASS" else "FAIL",
    if (!anyNA(thresholds[, .(p2_5, p97_5, mmt)])) "PASS" else "FAIL",
    if (all(thresholds$p2_5 < thresholds$p97_5)) "PASS" else "FAIL",
    if (all(range(year(obs_data$date)) == threshold_years)) "PASS" else "FAIL",
    if (all(is.finite(thresholds$mmt))) "PASS" else "FAIL"
  ),
  value = c(
    sprintf("%d cities", length(cities)),
    sprintf("%d rows", nrow(obs_data)),
    sprintf("%d rows", nrow(thresholds)),
    sprintf("min gap = %.2f C", min(thresholds$p97_5 - thresholds$p2_5)),
    paste(range(year(obs_data$date)), collapse = "-"),
    sprintf("%d of %d city-age MMTs differ from the 2023 value by > 0.5 C", sum(abs(thresholds$mmt - thresholds$mmt_2023) > 0.5), nrow(thresholds))
  )
)
fwrite(checks, checks_file)
if (any(checks$status == "FAIL")) {
  stop(sprintf("00a_prep_temperature.R failed %d check(s); see %s", sum(checks$status == "FAIL"), checks_file), call. = FALSE)
}

save(obs_data, thresholds, cities, file = out_file)
message("Saved ", out_file, " (", length(cities), " cities, ", nrow(obs_data), " daily rows)")
