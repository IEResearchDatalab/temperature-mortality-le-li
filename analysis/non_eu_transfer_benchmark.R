#!/usr/bin/env Rscript

################################################################################
# Leave-one-country-out benchmark for non-EU ERF coefficient transfer
#
# This is an analysis-layer precursor, not part of the five-step production
# pipeline. Every target country's cities are held out together. Candidate
# rules use only other countries' coefficient vectors. Nearest-analogue rules
# use geography or ERA5 temperature-distribution features, never the target
# outcome or target coefficients, to select donors.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(dlnm)
  library(splines)
  library(ggplot2)
})

coefs_file <- Sys.getenv("COEFS_FILE", "data/coefs.csv")
city_file <- Sys.getenv("CITY_FILE", "data/city_results.csv")
era5_file <- Sys.getenv("ERA5_FILE", "data/era5series.gz.parquet")
out_dir <- Sys.getenv(
  "OUT_DIR",
  "results/europe/auxiliary/non_eu_transfer_benchmark_20261008"
)
age_groups <- strsplit(
  Sys.getenv("AGE_GROUPS", "65-74,75-84,85+"),
  ",",
  fixed = TRUE
)[[1L]]
age_groups <- trimws(age_groups)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(coefs_file, city_file, era5_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Missing required input(s): ", paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

coef_cols <- paste0("b", 1:5)
required_city_cols <- c(
  "URAU_CODE", "LABEL", "CNTR_CODE", "cntr_name", "region", "lon", "lat"
)

coefs <- fread(coefs_file)
city_results <- fread(city_file)
if (!all(c("URAU_CODE", "agegroup", coef_cols) %in% names(coefs))) {
  stop("Coefficient file does not contain the expected city, age and b1-b5 columns.", call. = FALSE)
}
if (!all(required_city_cols %in% names(city_results))) {
  stop("City metadata file is missing required geographic fields.", call. = FALSE)
}

meta <- unique(city_results[, ..required_city_cols])
if (meta[, anyDuplicated(URAU_CODE)] > 0L) {
  stop("City metadata is not unique by URAU_CODE.", call. = FALSE)
}
if (anyNA(meta[, .(CNTR_CODE, region, lon, lat)])) {
  stop("City metadata contains missing country, region or coordinate fields.", call. = FALSE)
}

obs <- coefs[agegroup %in% age_groups]
obs <- merge(obs, meta, by = "URAU_CODE", all.x = TRUE, sort = FALSE)
if (nrow(obs) != uniqueN(obs$URAU_CODE) * length(age_groups)) {
  stop("Expected exactly one coefficient row per city and requested age group.", call. = FALSE)
}
if (anyNA(obs[, c("CNTR_CODE", "region", coef_cols), with = FALSE])) {
  stop("Coefficient/metadata join contains missing values.", call. = FALSE)
}
if (!all(is.finite(as.matrix(obs[, ..coef_cols])))) {
  stop("Observed coefficient vectors contain non-finite values.", call. = FALSE)
}
setnames(obs, coef_cols, paste0("obs_", coef_cols))

cities <- sort(unique(obs$URAU_CODE))
countries <- sort(unique(obs$CNTR_CODE))
message(sprintf(
  "Benchmark domain: %d cities, %d countries, %d age groups.",
  length(cities), length(countries), length(age_groups)
))

# ERA5 provides climate features for donor selection and the target-local
# percentile grid used to evaluate transferred coefficient vectors. The ERF
# basis matches pipeline/00_pkg_params.R and pipeline/00_prep_data.R.
predper <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))
varper <- c(10, 75, 90)
eval_keep <- predper >= 1 & predper <= 99
mmt_keep <- predper >= 25 & predper <= 99

message("Reading ERA5 and computing target-local climate features...")
era5 <- as.data.table(read_parquet(
  era5_file,
  col_select = c("URAU_CODE", "era5landtmean")
))
setnames(era5, "era5landtmean", "tmean")
era5 <- era5[URAU_CODE %in% cities]
if (!all(is.finite(era5$tmean))) {
  stop("ERA5 input contains non-finite temperatures.", call. = FALSE)
}

climate_features <- era5[, .(
  climate_mean = mean(tmean),
  climate_sd = sd(tmean),
  climate_p05 = as.numeric(quantile(tmean, 0.05, names = FALSE)),
  climate_p50 = as.numeric(quantile(tmean, 0.50, names = FALSE)),
  climate_p95 = as.numeric(quantile(tmean, 0.95, names = FALSE))
), by = URAU_CODE]
temp_grids <- era5[, .(
  temperature_grid = list(as.numeric(quantile(
    tmean,
    probs = predper / 100,
    names = FALSE,
    type = 7
  )))
), by = URAU_CODE]
rm(era5)
invisible(gc())

if (!setequal(cities, climate_features$URAU_CODE) ||
    !setequal(cities, temp_grids$URAU_CODE)) {
  stop("ERA5 does not cover every benchmark city.", call. = FALSE)
}
feature_cols <- setdiff(names(climate_features), "URAU_CODE")
feature_matrix <- as.matrix(climate_features[, ..feature_cols])
feature_scaled <- scale(feature_matrix)
if (any(!is.finite(feature_scaled))) {
  stop("Climate-feature standardisation produced non-finite values.", call. = FALSE)
}
climate_features[, paste0("z_", feature_cols) := as.data.table(feature_scaled)]
z_cols <- paste0("z_", feature_cols)
climate_features <- merge(
  climate_features,
  meta[, .(URAU_CODE, CNTR_CODE, region, lon, lat)],
  by = "URAU_CODE",
  all.x = TRUE,
  sort = FALSE
)

haversine_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  dlon <- (lon2 - lon1) * rad
  dlat <- (lat2 - lat1) * rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371.0088 * 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

nearest_rows <- rbindlist(lapply(seq_len(nrow(climate_features)), function(i) {
  target <- climate_features[i]
  candidates <- climate_features[CNTR_CODE != target$CNTR_CODE]
  if (!nrow(candidates)) {
    stop("No country-excluded donor candidates for ", target$URAU_CODE, call. = FALSE)
  }
  geo_distance <- haversine_km(
    target$lon, target$lat, candidates$lon, candidates$lat
  )
  climate_difference <- sweep(
    as.matrix(candidates[, ..z_cols]),
    2L,
    as.numeric(target[, ..z_cols]),
    "-"
  )
  climate_distance <- sqrt(rowSums(climate_difference^2))
  geo_i <- which.min(geo_distance)
  climate_i <- which.min(climate_distance)
  data.table(
    URAU_CODE = target$URAU_CODE,
    target_country = target$CNTR_CODE,
    geographic_donor = candidates$URAU_CODE[geo_i],
    geographic_donor_country = candidates$CNTR_CODE[geo_i],
    geographic_distance_km = geo_distance[geo_i],
    climate_donor = candidates$URAU_CODE[climate_i],
    climate_donor_country = candidates$CNTR_CODE[climate_i],
    climate_distance_z = climate_distance[climate_i]
  )
}))

observed_coef_cols <- paste0("obs_b", 1:5)
pred_coef_cols <- paste0("pred_b", 1:5)

# Country-excluded global means.
global_total <- obs[, c(
  list(total_n = .N),
  lapply(.SD, sum)
), by = agegroup, .SDcols = observed_coef_cols]
setnames(global_total, observed_coef_cols, paste0("total_", observed_coef_cols))
country_total <- obs[, c(
  list(country_n = .N),
  lapply(.SD, sum)
), by = .(CNTR_CODE, agegroup), .SDcols = observed_coef_cols]
setnames(country_total, observed_coef_cols, paste0("country_", observed_coef_cols))
global_pred <- merge(
  obs,
  merge(global_total, country_total, by = "agegroup", allow.cartesian = TRUE),
  by = c("CNTR_CODE", "agegroup"),
  all.x = TRUE,
  sort = FALSE
)
global_pred <- copy(global_pred)
for (j in seq_along(observed_coef_cols)) {
  global_pred[[pred_coef_cols[j]]] <-
    (global_pred[[paste0("total_", observed_coef_cols[j])]] -
       global_pred[[paste0("country_", observed_coef_cols[j])]]) /
    (global_pred$total_n - global_pred$country_n)
}
set(global_pred, j = "method", value = "country-excluded global mean")
set(global_pred, j = "donor_id", value = NA_character_)
set(global_pred, j = "donor_country", value = NA_character_)
set(global_pred, j = "donor_n", value = global_pred$total_n - global_pred$country_n)
set(global_pred, j = "donor_distance", value = NA_real_)

# Country-excluded mean within the target city's European region.
region_total <- obs[, c(
  list(total_n = .N),
  lapply(.SD, sum)
), by = .(region, agegroup), .SDcols = observed_coef_cols]
setnames(region_total, observed_coef_cols, paste0("total_", observed_coef_cols))
country_region_total <- obs[, c(
  list(country_n = .N),
  lapply(.SD, sum)
), by = .(CNTR_CODE, region, agegroup), .SDcols = observed_coef_cols]
setnames(country_region_total, observed_coef_cols, paste0("country_", observed_coef_cols))
region_pred <- merge(
  obs,
  merge(
    region_total,
    country_region_total,
    by = c("region", "agegroup"),
    allow.cartesian = TRUE
  ),
  by = c("CNTR_CODE", "region", "agegroup"),
  all.x = TRUE,
  sort = FALSE
)
region_pred <- copy(region_pred)
for (j in seq_along(observed_coef_cols)) {
  region_pred[[pred_coef_cols[j]]] <-
    (region_pred[[paste0("total_", observed_coef_cols[j])]] -
       region_pred[[paste0("country_", observed_coef_cols[j])]]) /
    (region_pred$total_n - region_pred$country_n)
}
set(region_pred, j = "method", value = "country-excluded regional mean")
set(region_pred, j = "donor_id", value = NA_character_)
set(region_pred, j = "donor_country", value = NA_character_)
set(region_pred, j = "donor_n", value = region_pred$total_n - region_pred$country_n)
set(region_pred, j = "donor_distance", value = NA_real_)

make_nearest_prediction <- function(map, donor_col, country_col, distance_col, method_name) {
  target <- merge(obs, map, by = "URAU_CODE", all.x = TRUE, sort = FALSE)
  donor <- copy(obs[, c(
    "URAU_CODE", "agegroup", "CNTR_CODE", observed_coef_cols
  ), with = FALSE])
  setnames(
    donor,
    c("URAU_CODE", "CNTR_CODE", observed_coef_cols),
    c("donor_id", "donor_country_from_coef", pred_coef_cols)
  )
  target[, donor_id := get(donor_col)]
  target[, donor_country := get(country_col)]
  target[, donor_distance := as.numeric(get(distance_col))]
  target <- merge(
    target,
    donor,
    by = c("donor_id", "agegroup"),
    all.x = TRUE,
    sort = FALSE
  )
  target[, `:=`(method = method_name, donor_n = 1L)]
  target
}

geo_pred <- make_nearest_prediction(
  nearest_rows,
  "geographic_donor",
  "geographic_donor_country",
  "geographic_distance_km",
  "nearest geographic city"
)
climate_pred <- make_nearest_prediction(
  nearest_rows,
  "climate_donor",
  "climate_donor_country",
  "climate_distance_z",
  "nearest ERA5 climate analogue"
)

keep_cols <- c(
  "URAU_CODE", "LABEL", "CNTR_CODE", "cntr_name", "region", "agegroup",
  observed_coef_cols, pred_coef_cols, "method", "donor_id", "donor_country",
  "donor_n", "donor_distance"
)
predictions <- rbindlist(list(
  global_pred[, ..keep_cols],
  region_pred[, ..keep_cols],
  geo_pred[, ..keep_cols],
  climate_pred[, ..keep_cols]
), use.names = TRUE)

setorder(predictions, CNTR_CODE, URAU_CODE, agegroup, method)

grid_lookup <- setNames(temp_grids$temperature_grid, temp_grids$URAU_CODE)
basis_lookup <- lapply(grid_lookup, function(tper) {
  onebasis(
    tper,
    fun = "bs",
    degree = 2,
    knots = tper[match(varper, predper)],
    Bound = range(tper)
  )
})
coef_sd <- obs[, lapply(.SD, sd), by = agegroup, .SDcols = observed_coef_cols]
coef_sd_lookup <- setNames(
  lapply(seq_len(nrow(coef_sd)), function(i) as.numeric(coef_sd[i, ..observed_coef_cols])),
  coef_sd$agegroup
)

evaluate_row <- function(i) {
  row <- predictions[i]
  city <- row$URAU_CODE
  age <- row$agegroup
  tper <- grid_lookup[[city]]
  basis <- basis_lookup[[city]]
  observed_coef <- as.numeric(row[, ..observed_coef_cols])
  predicted_coef <- as.numeric(row[, ..pred_coef_cols])
  observed_eta <- drop(basis %*% observed_coef)
  predicted_eta <- drop(basis %*% predicted_coef)
  observed_mmt_i <- which(mmt_keep)[which.min(observed_eta[mmt_keep])]
  predicted_mmt_i <- which(mmt_keep)[which.min(predicted_eta[mmt_keep])]
  observed_log_rr <- observed_eta - observed_eta[observed_mmt_i]
  predicted_log_rr <- predicted_eta - predicted_eta[predicted_mmt_i]
  curve_observed <- observed_log_rr[eval_keep]
  curve_predicted <- predicted_log_rr[eval_keep]
  curve_correlation <- if (
    sd(curve_observed) > 0 && sd(curve_predicted) > 0
  ) cor(curve_observed, curve_predicted) else NA_real_
  i01 <- which.min(abs(predper - 1))
  i99 <- which.min(abs(predper - 99))
  sd_vec <- coef_sd_lookup[[age]]
  data.table(
    coefficient_rmse = sqrt(mean((predicted_coef - observed_coef)^2)),
    scaled_coefficient_rmse = sqrt(mean(((predicted_coef - observed_coef) / sd_vec)^2)),
    log_rr_curve_rmse = sqrt(mean((curve_predicted - curve_observed)^2)),
    log_rr_curve_correlation = curve_correlation,
    observed_mmt_c = tper[observed_mmt_i],
    predicted_mmt_c = tper[predicted_mmt_i],
    mmt_abs_error_c = abs(tper[predicted_mmt_i] - tper[observed_mmt_i]),
    mmt_abs_error_percentile = abs(predper[predicted_mmt_i] - predper[observed_mmt_i]),
    cold_tail_log_rr_abs_error = abs(predicted_log_rr[i01] - observed_log_rr[i01]),
    heat_tail_log_rr_abs_error = abs(predicted_log_rr[i99] - observed_log_rr[i99])
  )
}

message(sprintf("Evaluating %s held-out city-age-method predictions...", nrow(predictions)))
metrics <- rbindlist(lapply(seq_len(nrow(predictions)), evaluate_row))
predictions <- cbind(predictions, metrics)

method_order <- c(
  "country-excluded regional mean",
  "nearest ERA5 climate analogue",
  "nearest geographic city",
  "country-excluded global mean"
)
predictions[, method := factor(method, levels = method_order)]

metric_summary <- function(x) {
  list(
    mean = mean(x),
    median = median(x),
    p25 = as.numeric(quantile(x, 0.25)),
    p75 = as.numeric(quantile(x, 0.75)),
    p90 = as.numeric(quantile(x, 0.90))
  )
}

method_age_summary <- predictions[, {
  curve <- metric_summary(log_rr_curve_rmse)
  mmt <- metric_summary(mmt_abs_error_c)
  coef <- metric_summary(scaled_coefficient_rmse)
  list(
    n = .N,
    log_rr_rmse_mean = curve$mean,
    log_rr_rmse_median = curve$median,
    log_rr_rmse_p25 = curve$p25,
    log_rr_rmse_p75 = curve$p75,
    log_rr_rmse_p90 = curve$p90,
    mmt_mae_c = mmt$mean,
    mmt_median_ae_c = mmt$median,
    mmt_ae_p25_c = mmt$p25,
    mmt_ae_p75_c = mmt$p75,
    mmt_ae_p90_c = mmt$p90,
    scaled_coef_rmse_mean = coef$mean,
    curve_correlation_median = median(log_rr_curve_correlation, na.rm = TRUE),
    cold_tail_log_rr_mae = mean(cold_tail_log_rr_abs_error),
    heat_tail_log_rr_mae = mean(heat_tail_log_rr_abs_error)
  )
}, by = .(agegroup, method)]

method_overall_summary <- predictions[, {
  curve <- metric_summary(log_rr_curve_rmse)
  mmt <- metric_summary(mmt_abs_error_c)
  coef <- metric_summary(scaled_coefficient_rmse)
  list(
    n = .N,
    log_rr_rmse_mean = curve$mean,
    log_rr_rmse_median = curve$median,
    log_rr_rmse_p25 = curve$p25,
    log_rr_rmse_p75 = curve$p75,
    log_rr_rmse_p90 = curve$p90,
    mmt_mae_c = mmt$mean,
    mmt_median_ae_c = mmt$median,
    mmt_ae_p25_c = mmt$p25,
    mmt_ae_p75_c = mmt$p75,
    mmt_ae_p90_c = mmt$p90,
    scaled_coef_rmse_mean = coef$mean,
    curve_correlation_median = median(log_rr_curve_correlation, na.rm = TRUE),
    cold_tail_log_rr_mae = mean(cold_tail_log_rr_abs_error),
    heat_tail_log_rr_mae = mean(heat_tail_log_rr_abs_error)
  )
}, by = method]
method_overall_summary[, agegroup := "All 65+"]
setcolorder(method_overall_summary, names(method_age_summary))
all_method_summary <- rbindlist(list(method_age_summary, method_overall_summary))

country_method_summary <- predictions[, .(
  n = .N,
  log_rr_rmse_mean = mean(log_rr_curve_rmse),
  mmt_mae_c = mean(mmt_abs_error_c),
  scaled_coef_rmse_mean = mean(scaled_coefficient_rmse)
), by = .(CNTR_CODE, cntr_name, method)]
country_method_summary[, relative_curve_error :=
  log_rr_rmse_mean / min(log_rr_rmse_mean), by = CNTR_CODE]
country_method_summary[, method_rank :=
  frank(log_rr_rmse_mean, ties.method = "average"), by = CNTR_CODE]

rank_stability <- copy(method_age_summary)
rank_stability[, `:=`(
  curve_rank = frank(log_rr_rmse_mean, ties.method = "average"),
  mmt_rank = frank(mmt_mae_c, ties.method = "average")
), by = agegroup]
overall_ranks <- copy(method_overall_summary)
overall_ranks[, `:=`(
  curve_rank = frank(log_rr_rmse_mean, ties.method = "average"),
  mmt_rank = frank(mmt_mae_c, ties.method = "average")
)]
rank_stability <- rbindlist(list(rank_stability, overall_ranks), use.names = TRUE)
rank_range <- rank_stability[agegroup != "All 65+", .(
  curve_rank_min = min(curve_rank),
  curve_rank_max = max(curve_rank),
  curve_rank_range = max(curve_rank) - min(curve_rank),
  mmt_rank_min = min(mmt_rank),
  mmt_rank_max = max(mmt_rank),
  mmt_rank_range = max(mmt_rank) - min(mmt_rank)
), by = method]
rank_stability <- merge(rank_stability, rank_range, by = "method", all.x = TRUE)

age_pairs <- combn(age_groups, 2L, simplify = FALSE)
rank_correlations <- rbindlist(lapply(age_pairs, function(pair) {
  a <- rank_stability[agegroup == pair[1], .(method, curve_rank, mmt_rank)]
  b <- rank_stability[agegroup == pair[2], .(method, curve_rank, mmt_rank)]
  x <- merge(a, b, by = "method", suffixes = c("_a", "_b"))
  data.table(
    agegroup_a = pair[1],
    agegroup_b = pair[2],
    curve_rank_spearman = cor(x$curve_rank_a, x$curve_rank_b, method = "spearman"),
    mmt_rank_spearman = cor(x$mmt_rank_a, x$mmt_rank_b, method = "spearman")
  )
}))

expected_rows <- length(cities) * length(age_groups) * length(method_order)
nearest_only <- predictions[method %in% c(
  "nearest geographic city", "nearest ERA5 climate analogue"
)]
validation_checks <- data.table(
  check_name = c(
    "expected_city_age_method_rows",
    "one_row_per_city_age_method",
    "all_854_cities_present",
    "all_countries_present",
    "all_age_groups_present",
    "all_methods_present",
    "era5_complete",
    "five_coefficient_dimensions",
    "regional_pools_nonempty",
    "nearest_donor_country_excluded",
    "nearest_donor_identity_excluded",
    "nearest_donor_country_matches_coefficients",
    "predicted_coefficients_finite",
    "evaluation_metrics_finite",
    "climate_selection_uses_no_erf_outcome"
  ),
  status = c(
    if (nrow(predictions) == expected_rows) "PASS" else "FAIL",
    if (predictions[, anyDuplicated(paste(URAU_CODE, agegroup, method))] == 0L) "PASS" else "FAIL",
    if (uniqueN(predictions$URAU_CODE) == 854L) "PASS" else "FAIL",
    if (uniqueN(predictions$CNTR_CODE) == length(countries)) "PASS" else "FAIL",
    if (setequal(unique(predictions$agegroup), age_groups)) "PASS" else "FAIL",
    if (setequal(as.character(unique(predictions$method)), method_order)) "PASS" else "FAIL",
    if (nrow(climate_features) == length(cities)) "PASS" else "FAIL",
    if (length(observed_coef_cols) == 5L && length(pred_coef_cols) == 5L) "PASS" else "FAIL",
    if (all(predictions[method == "country-excluded regional mean"]$donor_n > 0L)) "PASS" else "FAIL",
    if (all(nearest_only$donor_country != nearest_only$CNTR_CODE)) "PASS" else "FAIL",
    if (all(nearest_only$donor_id != nearest_only$URAU_CODE)) "PASS" else "FAIL",
    if (
      all(geo_pred$donor_country == geo_pred$donor_country_from_coef) &&
      all(climate_pred$donor_country == climate_pred$donor_country_from_coef)
    ) "PASS" else "FAIL",
    if (all(is.finite(as.matrix(predictions[, ..pred_coef_cols])))) "PASS" else "FAIL",
    if (all(is.finite(as.matrix(predictions[, .(
      coefficient_rmse, scaled_coefficient_rmse, log_rr_curve_rmse,
      log_rr_curve_correlation,
      observed_mmt_c, predicted_mmt_c, mmt_abs_error_c,
      mmt_abs_error_percentile, cold_tail_log_rr_abs_error,
      heat_tail_log_rr_abs_error
    )])))) "PASS" else "FAIL",
    "PASS"
  ),
  detail = c(
    sprintf("observed=%d expected=%d", nrow(predictions), expected_rows),
    "unique target city x age group x method",
    sprintf("observed=%d", uniqueN(predictions$URAU_CODE)),
    sprintf("observed=%d", uniqueN(predictions$CNTR_CODE)),
    paste(sort(unique(predictions$agegroup)), collapse = ", "),
    paste(as.character(unique(predictions$method)), collapse = ", "),
    sprintf("ERA5 features for %d cities", nrow(climate_features)),
    "b1-b5 for observed and transferred vectors",
    sprintf("minimum donor city-age rows=%d", min(predictions[
      method == "country-excluded regional mean"
    ]$donor_n)),
    "all nearest donors belong to another country",
    "no city can donate to itself",
    "donor-map country agrees with coefficient-table country",
    "all transferred b1-b5 values finite",
    "all primary numerical metrics finite",
    "donors selected only from region, coordinates or ERA5 distribution features"
  )
)
if (any(validation_checks$status != "PASS")) {
  fwrite(validation_checks, file.path(out_dir, "validation_checks.csv"))
  stop("Benchmark failed one or more validation checks.", call. = FALSE)
}

predictions[, method := as.character(method)]
all_method_summary[, method := as.character(method)]
country_method_summary[, method := as.character(method)]
rank_stability[, method := as.character(method)]

fwrite(predictions, file.path(out_dir, "city_age_predictions.csv"))
fwrite(all_method_summary, file.path(out_dir, "method_age_summary.csv"))
fwrite(country_method_summary, file.path(out_dir, "country_method_summary.csv"))
fwrite(rank_stability, file.path(out_dir, "method_rank_stability.csv"))
fwrite(rank_correlations, file.path(out_dir, "rank_correlations.csv"))
fwrite(validation_checks, file.path(out_dir, "validation_checks.csv"))
fwrite(nearest_rows, file.path(out_dir, "nearest_donor_map.csv"))

method_labels <- c(
  "country-excluded regional mean" = "Regional mean",
  "nearest ERA5 climate analogue" = "Nearest climate analogue",
  "nearest geographic city" = "Nearest geographic city",
  "country-excluded global mean" = "Global mean"
)
age_colors <- c("65-74" = "#4B8DC1", "75-84" = "#D97A3A", "85+" = "#7A5BA7")

plot_data <- rbindlist(list(
  method_age_summary[, .(
    agegroup, method,
    metric = "Log-RR curve RMSE",
    median = log_rr_rmse_median,
    p25 = log_rr_rmse_p25,
    p75 = log_rr_rmse_p75
  )],
  method_age_summary[, .(
    agegroup, method,
    metric = "MMT absolute error (°C)",
    median = mmt_median_ae_c,
    p25 = mmt_ae_p25_c,
    p75 = mmt_ae_p75_c
  )]
))
plot_data[, method_label := factor(
  method_labels[as.character(method)],
  levels = rev(unname(method_labels))
)]
p1 <- ggplot(plot_data, aes(x = median, y = method_label, color = agegroup)) +
  geom_errorbar(
    aes(xmin = p25, xmax = p75),
    width = 0.22,
    linewidth = 0.7,
    orientation = "y"
  ) +
  geom_point(size = 2.7) +
  facet_wrap(~metric, scales = "free_x", ncol = 2) +
  scale_color_manual(values = age_colors) +
  labs(
    title = "Leave-one-country-out transfer performance",
    subtitle = "Median and interquartile range across held-out cities; lower is better",
    x = NULL, y = NULL, color = "Age group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 18),
    legend.position = "bottom"
  )
ggsave(
  file.path(out_dir, "01_transfer_performance_by_age.png"),
  p1, width = 12, height = 6.8, dpi = 200, bg = "white"
)

heat <- copy(country_method_summary)
heat[, method_label := method_labels[as.character(method)]]
country_order <- heat[, .(best_error = min(log_rr_rmse_mean)), by = CNTR_CODE][
  order(best_error), CNTR_CODE
]
heat[, country_axis := factor(CNTR_CODE, levels = rev(country_order))]
p2 <- ggplot(heat, aes(x = method_label, y = country_axis, fill = relative_curve_error)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(
    low = "#EDF8E9", high = "#B2182B",
    name = "Error / country best",
    labels = function(x) sprintf("%.2fx", x)
  ) +
  labs(
    title = "Transfer-rule performance varies across held-out countries",
    subtitle = "Mean log-RR curve RMSE across ages 65+; 1.00x is the best rule within each country",
    x = NULL, y = "Held-out country"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold", size = 17),
    legend.position = "right"
  )
ggsave(
  file.path(out_dir, "02_country_method_heatmap.png"),
  p2, width = 10.5, height = 9.5, dpi = 200, bg = "white"
)

git_commit <- tryCatch(
  system("git rev-parse HEAD", intern = TRUE),
  error = function(e) "unavailable"
)
git_branch <- tryCatch(
  system("git branch --show-current", intern = TRUE),
  error = function(e) "unavailable"
)
host <- Sys.info()[["nodename"]]
input_manifest <- data.table(
  role = c("central coefficients", "city metadata", "ERA5 daily temperature"),
  path = normalizePath(c(coefs_file, city_file, era5_file)),
  bytes = file.info(c(coefs_file, city_file, era5_file))$size,
  md5 = unname(tools::md5sum(c(coefs_file, city_file, era5_file)))
)
fwrite(input_manifest, file.path(out_dir, "input_manifest.csv"))

overall_print <- method_overall_summary[order(log_rr_rmse_mean), .(
  method = as.character(method),
  log_rr_rmse_mean,
  mmt_mae_c,
  curve_correlation_median
)]
top_method <- overall_print$method[1]
top_curve <- overall_print$log_rr_rmse_mean[1]
top_mmt_method <- method_overall_summary[which.min(mmt_mae_c), as.character(method)]
top_mmt <- min(method_overall_summary$mmt_mae_c)
rank_stable <- rank_range[curve_rank_range == 0, as.character(method)]
rank_stable_text <- if (length(rank_stable)) {
  paste(rank_stable, collapse = ", ")
} else {
  "No method held exactly the same curve-error rank in all three age groups."
}

command_text <- paste0(
  "ERA5_FILE=", normalizePath(era5_file), " ",
  "OUT_DIR=", out_dir, " ",
  "Rscript analysis/non_eu_transfer_benchmark.R"
)
readme <- c(
  "# Non-EU transfer precursor: leave-one-country-out benchmark",
  "",
  sprintf("Generated on `%s` at `%s`.", host, format(Sys.time(), tz = "UTC")),
  sprintf("Repository: `%s`", normalizePath(getwd())),
  sprintf("Branch: `%s`", git_branch),
  sprintf("Commit: `%s`", git_commit),
  "",
  "## Command",
  "",
  paste0("`", command_text, "`"),
  "",
  "## Design",
  "",
  paste0(
    "All cities in one country are held out together. Four transfer rules predict each held-out ",
    "city's five central ERF coefficients separately for ages 65-74, 75-84 and 85+: ",
    "country-excluded regional mean, country-excluded global mean, nearest geographic city, ",
    "and nearest ERA5 climate analogue. Nearest donors always come from another country."
  ),
  "",
  paste0(
    "Transferred vectors are evaluated on the held-out city's local percentile spline basis. ",
    "Primary metrics are log-relative-risk curve RMSE across temperature percentiles 1-99 and ",
    "absolute MMT error. Coefficient RMSE and cold/heat tail errors are secondary diagnostics."
  ),
  "",
  "## Headline result",
  "",
  sprintf(
    "The lowest overall mean log-RR curve RMSE was `%s` (%.6f).",
    top_method, top_curve
  ),
  sprintf(
    "The lowest overall mean MMT absolute error was `%s` (%.3f °C).",
    top_mmt_method, top_mmt
  ),
  sprintf("Age-rank stability: %s", rank_stable_text),
  "",
  "See `method_age_summary.csv`, `method_rank_stability.csv` and the two figures for the full result.",
  "",
  "## Validation",
  "",
  sprintf("All %d validation checks passed.", nrow(validation_checks)),
  sprintf(
    "Coverage: %d cities × %d age groups × %d methods = %d evaluated predictions.",
    length(cities), length(age_groups), length(method_order), nrow(predictions)
  ),
  "",
  "## Outputs",
  "",
  "- `city_age_predictions.csv`: prediction and error for every city-age-method row.",
  "- `method_age_summary.csv`: method performance overall and by age group.",
  "- `country_method_summary.csv`: country-level performance and within-country ranks.",
  "- `method_rank_stability.csv`: method ranks across age groups.",
  "- `rank_correlations.csv`: Spearman rank correlations between age groups.",
  "- `nearest_donor_map.csv`: country-excluded geographic and climate donors.",
  "- `validation_checks.csv`: leakage, coverage, dimension and finiteness checks.",
  "- `input_manifest.csv`: exact input paths, sizes and MD5 hashes.",
  "- `01_transfer_performance_by_age.png`: primary performance figure.",
  "- `02_country_method_heatmap.png`: country heterogeneity diagnostic.",
  "",
  "## Limitations",
  "",
  paste0(
    "This tests internal European transportability, not Türkiye. The benchmark transfers central ",
    "coefficient vectors and does not propagate coefficient covariance. ERA5 climate similarity ",
    "uses distribution summaries only and omits income, vulnerability, healthcare, air conditioning ",
    "and mortality-system covariates. Regional labels are European and may not map naturally outside ",
    "the study domain. Performance against held-out central estimates does not establish external ",
    "validity or justify narrow uncertainty intervals."
  )
)
writeLines(readme, file.path(out_dir, "README.md"), useBytes = TRUE)

message("Completed. Outputs: ", normalizePath(out_dir))
print(overall_print)
print(validation_checks)
