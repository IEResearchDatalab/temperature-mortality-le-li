#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(arrow); library(dplyr)})

out_root <- Sys.getenv("UNCERTAINTY_ROOT", "results/europe/uncertainty_20261104")
timing_root <- Sys.getenv(
  "TIMING500_ROOT",
  file.path(out_root, "timing500", "part01", "ssp3", "GFDL_ESM4")
)
dir.create(file.path(out_root, "checks"), recursive = TRUE, showWarnings = FALSE)

draws <- as.data.table(read_parquet("references/2025-masselot-data/coef_simu.gz.parquet"))
draws <- draws[agegroup %chin% c("65-74", "75-84", "85+")]
coef_cols <- paste0("b", 1:5)
available <- draws[, .(available_draws = uniqueN(sim), min_available = min(sim),
                       max_available = max(sim)), by = .(URAU_CODE, agegroup)]
draws <- draws[sim %between% c(1L, 500L)]
coverage <- draws[, .(
  rows = .N,
  n_draws = uniqueN(sim),
  min_draw = min(sim),
  max_draw = max(sim),
  duplicate_draws = .N - uniqueN(sim),
  finite_coefficients = all(vapply(.SD, function(x) all(is.finite(x)), logical(1)))
), by = .(URAU_CODE, agegroup), .SDcols = coef_cols]
coverage <- merge(coverage, available, by = c("URAU_CODE", "agegroup"))
fwrite(coverage, file.path(out_root, "checks", "coefficient_draw_coverage.csv"))
if (nrow(coverage) != 854L * 3L || any(coverage$n_draws != 500L) ||
    any(coverage$rows != 500L) || any(coverage$available_draws < 500L) ||
    any(coverage$duplicate_draws != 0L) ||
    any(!coverage$finite_coefficients)) {
  stop("Coefficient draw coverage failed", call. = FALSE)
}

central <- fread("data/coefs.csv")[agegroup %chin% c("65-74", "75-84", "85+")]
medians <- draws[, lapply(.SD, median), by = .(URAU_CODE, agegroup), .SDcols = coef_cols]
coef_compare <- merge(central, medians, by = c("URAU_CODE", "agegroup"), suffixes = c("_central", "_draw_median"))
for (b in coef_cols) {
  coef_compare[, paste0(b, "_difference") :=
                 get(paste0(b, "_draw_median")) - get(paste0(b, "_central"))]
}
fwrite(coef_compare, file.path(out_root, "checks", "coefficient_medians_vs_central.csv"))

pooled <- as.data.table(read_parquet(file.path(timing_root, "pooled_counts.parquet")))
pooled[, `:=`(draw = as.integer(draw), year = as.integer(year))]
median_counts <- pooled[geography == "Europe", .(draw_median_an = median(an)),
                        by = .(branch, year, agegroup, range)]
canonical <- fread(
  "results/europe/heatsplit_20261013/ssp3/ES001C/GFDL_ESM4/01_attribution_grouped.csv"
)[, .(branch, year, agegroup, range, central_an = an)]
count_compare <- merge(canonical, median_counts, by = c("branch", "year", "agegroup", "range"))
count_compare[, `:=`(
  difference = draw_median_an - central_an,
  relative_error = fifelse(abs(central_an) > 1e-8,
                           abs(draw_median_an / central_an - 1), NA_real_)
)]
fwrite(count_compare, file.path(out_root, "checks", "draw_median_vs_central_madrid.csv"))

aggregate_draw <- pooled[geography == "Europe", .(an = sum(an)),
                        by = .(draw, branch, agegroup)]
aggregate_draw <- aggregate_draw[, .(draw_median_an = median(an)),
                                 by = .(branch, agegroup)]
aggregate_central <- canonical[, .(central_an = sum(central_an)),
                               by = .(branch, agegroup)]
aggregate_compare <- merge(aggregate_central, aggregate_draw,
                           by = c("branch", "agegroup"))
aggregate_compare[, `:=`(
  difference = draw_median_an - central_an,
  relative_error = abs(draw_median_an / central_an - 1)
)]
fwrite(aggregate_compare,
       file.path(out_root, "checks", "draw_median_vs_central_madrid_aggregate.csv"))

summary <- data.table(
  check = c(
    "city_age_groups", "draws_per_city_age", "coefficient_basis_columns",
    "all_coefficients_finite", "median_count_relative_error_median",
    "median_count_relative_error_p95", "median_count_relative_error_max_nonzero",
    "aggregate_count_relative_error_max"
  ),
  value = c(
    nrow(coverage), unique(coverage$n_draws), length(coef_cols),
    all(coverage$finite_coefficients),
    median(count_compare$relative_error, na.rm = TRUE),
    quantile(count_compare$relative_error, 0.95, na.rm = TRUE),
    max(count_compare$relative_error, na.rm = TRUE),
    max(aggregate_compare$relative_error)
  ),
  verdict = c(
    rep("match", 4),
    rep("reported; cellwise errors are unstable when central counts are small", 3),
    if (max(aggregate_compare$relative_error) < 0.05) "match" else "not explained"
  )
)
fwrite(summary, file.path(out_root, "checks", "coefficient_draw_validation_summary.csv"))
print(summary)
