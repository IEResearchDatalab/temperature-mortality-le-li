#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(parallel)})

canonical <- Sys.getenv("CANONICAL_ROOT", "results/europe/heatsplit_20261013")
sensitivity_root <- Sys.getenv(
  "THRESHOLD_ROOT", "results/europe/threshold_sensitivities_20261104"
)
out_root <- Sys.getenv("CENTRAL_EXTRAS_ROOT", "results/europe/central_extras_20261104")
ncores <- as.integer(Sys.getenv("THRESHOLD_SUMMARY_CORES", "16"))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

meta <- unique(fread("data/city_results.csv", select = c("URAU_CODE", "region")))
cities <- sort(meta$URAU_CODE)
gcms <- c(
  "ACCESS_CM2", "ACCESS_ESM1_5", "BCC_CSM2_MR", "CanESM5", "CMCC_ESM2",
  "EC_Earth3", "EC_Earth3_Veg_LR", "GFDL_ESM4", "IITM_ESM", "INM_CM4_8",
  "INM_CM5_0", "IPSL_CM6A_LR", "KACE_1_0_G", "MIROC6", "MPI_ESM1_2_HR",
  "MPI_ESM1_2_LR", "MRI_ESM2_0", "NorESM2_LM", "NorESM2_MM"
)
period_years <- list(
  `2020-2024` = 2020:2024,
  `2055-2059` = 2055:2059,
  `2095-2099` = 2095:2099
)

read_city <- function(root, label, ssp, city) {
  files <- if (label == "canonical") {
    file.path(root, paste0("ssp", ssp), city, gcms, "01_attribution_grouped.csv")
  } else {
    file.path(root, label, paste0("ssp", ssp), city, gcms, "01_attribution_grouped.csv")
  }
  if (any(!file.exists(files))) stop("Missing Part 01 files for ", label, " SSP", ssp, " ", city)
  x <- rbindlist(lapply(seq_along(files), function(i) {
    z <- fread(files[i], select = c("year", "agegroup", "range", "branch", "an"))
    z <- z[year %in% unlist(period_years, use.names = FALSE)]
    z[, gcm := gcms[i]]
    z
  }))
  x <- x[, .(an = mean(an)), by = .(year, agegroup, range, branch)]
  x[, `:=`(city = city, ssp = ssp, variant = label)]
  x
}

variants <- data.table(
  variant = c("canonical", "trailing30_p975", "trailing15_p95"),
  root = c(canonical, sensitivity_root, sensitivity_root)
)
cases <- CJ(v = seq_len(nrow(variants)), ssp = 1:3, city = cities)
parts <- mclapply(seq_len(nrow(cases)), function(i) {
  v <- variants[cases$v[i]]
  read_city(v$root, v$variant, cases$ssp[i], cases$city[i])
}, mc.cores = ncores, mc.preschedule = FALSE)
if (any(vapply(parts, inherits, logical(1), what = "try-error"))) {
  stop("One or more threshold summary inputs failed", call. = FALSE)
}
x <- rbindlist(parts)
x <- merge(x, meta, by.x = "city", by.y = "URAU_CODE")
x[, period := fifelse(year %in% period_years[[1]], names(period_years)[1],
               fifelse(year %in% period_years[[2]], names(period_years)[2], names(period_years)[3]))]

pool <- rbindlist(list(
  x[, .(an = sum(an)), by = .(variant, ssp, geography = region, period, branch, range)],
  x[, .(an = sum(an)), by = .(variant, ssp, period, branch, range)][, geography := "Europe"]
), use.names = TRUE)

metrics <- pool[, {
  with <- .SD[branch == "with_cc"]
  without <- .SD[branch == "without_cc"]
  with_total <- with[, sum(an)]
  with_heat <- with[range %chin% c("ModHeat", "NormExtrHeat", "ExtrHeat"), sum(an)]
  without_heat <- without[range %chin% c("ModHeat", "NormExtrHeat", "ExtrHeat"), sum(an)]
  added_heat <- with_heat - without_heat
  normal_extreme_added <- with[range == "NormExtrHeat", sum(an)] -
    without[range == "NormExtrHeat", sum(an)]
  list(
    temperature_deaths = with_total,
    added_heat_deaths = added_heat,
    normal_extreme_added_heat_deaths = normal_extreme_added,
    normal_extreme_share_of_added_heat_pct = 100 * normal_extreme_added / added_heat,
    sustained_tier_share_pct = 100 * with[range == "NormExtrHeat", sum(an)] / with_total
  )
}, by = .(variant, ssp, geography, period)]

wide <- dcast(
  metrics,
  ssp + geography + period ~ variant,
  value.var = c("normal_extreme_share_of_added_heat_pct", "sustained_tier_share_pct")
)
wide[, `:=`(
  normal_extreme_30yr_minus_canonical_pp =
    normal_extreme_share_of_added_heat_pct_trailing30_p975 -
    normal_extreme_share_of_added_heat_pct_canonical,
  normal_extreme_p95_minus_canonical_pp =
    normal_extreme_share_of_added_heat_pct_trailing15_p95 -
    normal_extreme_share_of_added_heat_pct_canonical,
  sustained_30yr_minus_canonical_pp = sustained_tier_share_pct_trailing30_p975 -
    sustained_tier_share_pct_canonical,
  sustained_p95_minus_canonical_pp = sustained_tier_share_pct_trailing15_p95 -
    sustained_tier_share_pct_canonical
)]
fwrite(metrics, file.path(out_root, "threshold_sensitivity_metrics.csv"))
fwrite(wide, file.path(out_root, "threshold_sensitivity_comparison.csv"))
message("Wrote threshold sensitivity summaries to ", out_root)
