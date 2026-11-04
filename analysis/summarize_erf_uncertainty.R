#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(arrow)})

root <- Sys.getenv("UNCERTAINTY_ROOT", "results/europe/uncertainty_20261104")
ssp <- as.integer(Sys.getenv("SSP", "3"))
table_dir <- file.path(root, "tables", paste0("ssp", ssp))
check_dir <- file.path(root, "checks")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

gcms <- scan(file.path(check_dir, "gcms.txt"), what = character(), quiet = TRUE)
metric_files <- file.path(root, "part01", paste0("ssp", ssp), gcms,
                          "experiment_metrics.parquet")
risk_files <- file.path(root, "part01", paste0("ssp", ssp), gcms,
                        "annual_risk.parquet")
if (length(gcms) != 19L || any(!file.exists(metric_files)) || any(!file.exists(risk_files))) {
  stop("Expected 19 complete GCM outputs", call. = FALSE)
}

metrics <- rbindlist(lapply(metric_files, function(f) as.data.table(read_parquet(f))))
annual <- rbindlist(lapply(risk_files, function(f) as.data.table(read_parquet(f))))
write_parquet(metrics, file.path(table_dir, "experiment_level.parquet"), compression = "zstd")
write_parquet(annual, file.path(table_dir, "annual_risk_experiment_level.parquet"), compression = "zstd")

id_cols <- c("ssp", "gcm", "draw", "geography")
estimand_cols <- c(
  "lifetime_risk_all_pct", "lifetime_risk_heat_pct", "persistent_share_pct",
  "general_tier_pct", "sustained_tier_pct", "episodic_tier_pct",
  "temperature_death_rate_65plus_end", "climate_death_rate_65plus_end",
  "heat_erosion_ratio_pct", "all_erosion_ratio_pct",
  "delta_LE65_end", "delta_LI65_end"
)
long <- melt(metrics, id.vars = id_cols, measure.vars = estimand_cols,
             variable.name = "estimand", value.name = "value")

combined <- long[, .(
  experiment_median = median(value),
  combined_low = quantile(value, 0.025),
  combined_high = quantile(value, 0.975),
  n_experiments = .N
), by = .(ssp, geography, estimand)]
erf_by_draw <- long[, .(value = mean(value)),
                    by = .(ssp, geography, estimand, draw)]
erf <- erf_by_draw[, .(
  erf_only_low = quantile(value, 0.025),
  erf_only_high = quantile(value, 0.975),
  erf_only_median = median(value)
), by = .(ssp, geography, estimand)]
gcm_by_model <- long[, .(value = mean(value)),
                     by = .(ssp, geography, estimand, gcm)]
gcm <- gcm_by_model[, .(
  gcm_low = quantile(value, 0.025),
  gcm_high = quantile(value, 0.975),
  gcm_min = min(value), gcm_max = max(value), gcm_median = median(value)
), by = .(ssp, geography, estimand)]
summary <- Reduce(function(x, y) merge(x, y, by = c("ssp", "geography", "estimand")),
                  list(combined, erf, gcm))

# Published central values for Table 1 metrics.
table1 <- fread("results/figures_simon_20oct/table1/table1_full.csv")[
  ssp == as.integer(Sys.getenv("SSP", "3")) & level %chin% c("europe", "region")
]
central <- rbindlist(list(
  table1[, .(ssp, geography, estimand = "lifetime_risk_all_pct", central = lifetime_risk_pct)],
  table1[, .(ssp, geography, estimand = "lifetime_risk_heat_pct", central = lifetime_risk_heat_pct)],
  table1[, .(ssp, geography, estimand = "persistent_share_pct", central = persistent_share_pct)],
  table1[, .(ssp, geography, estimand = "general_tier_pct", central = general_pct)],
  table1[, .(ssp, geography, estimand = "sustained_tier_pct", central = sustained_pct)],
  table1[, .(ssp, geography, estimand = "episodic_tier_pct", central = episodic_pct)],
  table1[, .(ssp, geography, estimand = "heat_erosion_ratio_pct", central = heat_ratio_median_pct)],
  table1[, .(ssp, geography, estimand = "all_erosion_ratio_pct", central = all_ratio_pct)]
))

geo_map <- data.table(
  geography = c("Europe", "Eastern Europe", "Northern Europe", "Southern Europe", "Western Europe"),
  level = c("europe", rep("region", 4)),
  id = c("EUROPE", "Eastern", "Northern", "Southern", "Western")
)
levels <- rbindlist(lapply(seq_len(nrow(geo_map)), function(i) {
  f <- file.path("results/europe/heatsplit_20261013/geographies", paste0("ssp", ssp),
                 geo_map$level[i], geo_map$id[i], "04_le_li_levels.csv")
  x <- fread(f)[year %between% c(2095L, 2099L),
                .(LE65 = mean(LE65), LI65 = mean(LI65)), by = branch]
  data.table(
    ssp = ssp, geography = geo_map$geography[i],
    delta_LE65_end = x[branch == "with_cc", LE65] - x[branch == "without_cc", LE65],
    delta_LI65_end = x[branch == "with_cc", LI65] - x[branch == "without_cc", LI65]
  )
}))
central <- rbind(
  central,
  melt(levels, id.vars = c("ssp", "geography"), variable.name = "estimand",
       value.name = "central"),
  fill = TRUE
)
summary <- merge(summary, central, by = c("ssp", "geography", "estimand"), all.x = TRUE)
summary[is.na(central), central := experiment_median]
summary[, central_source := fifelse(
  estimand %chin% c("temperature_death_rate_65plus_end", "climate_death_rate_65plus_end"),
  "experiment median (no matching published central table column)",
  "canonical central pipeline"
)]
signs <- merge(
  long,
  summary[, .(ssp, geography, estimand, central)],
  by = c("ssp", "geography", "estimand")
)[, .(sign_consistency = mean(sign(value) == sign(central))),
  by = .(ssp, geography, estimand)]
summary <- merge(summary, signs, by = c("ssp", "geography", "estimand"))
fwrite(summary, file.path(table_dir, "summary.csv"))

conv <- rbindlist(lapply(c(50L, 100L, 250L, 500L), function(n) {
  long[draw <= n, .(
    low = quantile(value, 0.025), high = quantile(value, 0.975)
  ), by = .(ssp, geography, estimand)][, draws := n]
}))
final <- conv[draws == 500L, .(ssp, geography, estimand, low500 = low, high500 = high)]
conv <- merge(conv, final, by = c("ssp", "geography", "estimand"))
conv[, `:=`(low_difference_from_500 = low - low500,
            high_difference_from_500 = high - high500)]
fwrite(conv, file.path(table_dir, "convergence.csv"))

annual_summary <- annual[, .(
  median = median(risk_per_100000),
  low = quantile(risk_per_100000, 0.025),
  high = quantile(risk_per_100000, 0.975),
  sign_consistency = mean(sign(risk_per_100000) == sign(median(risk_per_100000)))
), by = .(ssp, geography, year, age_scope, estimand)]
write_parquet(annual_summary, file.path(table_dir, "annual_risk_summary.parquet"),
              compression = "zstd")

checks <- data.table(
  check = c("gcms", "draws", "geographies", "metric_experiments",
            "annual_rows", "finite_metrics", "finite_annual"),
  value = c(uniqueN(metrics$gcm), uniqueN(metrics$draw), uniqueN(metrics$geography),
            nrow(metrics), nrow(annual),
            all(vapply(metrics, function(x) !is.numeric(x) || all(is.finite(x)), logical(1))),
            all(is.finite(annual$risk_per_100000))),
  expected = c("19", "500", "5", "47500", "15200000", "TRUE", "TRUE")
)
fwrite(checks, file.path(table_dir, "checks.csv"))
message("Wrote uncertainty summaries to ", table_dir)
