#!/usr/bin/env Rscript

################################################################################
# Summarise ERF-draw and GCM uncertainty for one SSP.
# Requires: UNCERTAINTY_ROOT/checks/gcms.txt and
# UNCERTAINTY_ROOT/part01/ssp<SSP>/<GCM>/{experiment_metrics,annual_risk}.parquet
# from analysis/run_erf_uncertainty.sh; Table 1 table1_full.csv and pooled
# <batch>/geographies/ssp<SSP>/<level>/<id>/04_le_li_levels.csv.
# Produces: UNCERTAINTY_OUTPUT_ROOT/tables/ssp<SSP>/ summaries and checks.
# How to run: UNCERTAINTY_ROOT=results/europe/uncertainty_20261110 \
#   UNCERTAINTY_OUTPUT_ROOT=results/paper_output_new/uncertainty \
#   BATCH_ROOT=results/europe/heatsplit_20261013 \
#   TABLE1_FULL=results/paper_output_new/table1/table1_full.csv SSP=3 \
#   Rscript analysis/summarize_erf_uncertainty.R
# Notes: SSP1 replacement and SSP1/SSP2 shape-check bypasses are exploratory;
#   see the caveats in README.md.
################################################################################

suppressPackageStartupMessages({library(data.table); library(arrow)})

root <- Sys.getenv("UNCERTAINTY_ROOT", "results/europe/uncertainty_20261104")
output_root <- Sys.getenv("UNCERTAINTY_OUTPUT_ROOT", root)
batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
table1_file <- Sys.getenv("TABLE1_FULL", "results/figures_simon_20oct/table1/table1_full.csv")
ssp <- as.integer(Sys.getenv("SSP", "3"))
if (length(ssp) != 1L || is.na(ssp) || !ssp %in% c(2L, 3L)) {
  stop("SSP must be 2 or 3 for the available ERF uncertainty runs.", call. = FALSE)
}
table_dir <- file.path(output_root, "tables", paste0("ssp", ssp))
check_dir <- file.path(root, "checks")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

gcm_override <- Sys.getenv("GCM_IDS", "")
gcms <- if (nzchar(gcm_override)) {
  strsplit(gcm_override, ",", fixed = TRUE)[[1]]
} else {
  scan(file.path(check_dir, "gcms.txt"), what = character(), quiet = TRUE)
}
metric_files <- file.path(root, "part01", paste0("ssp", ssp), gcms,
                          "experiment_metrics.parquet")
risk_files <- file.path(root, "part01", paste0("ssp", ssp), gcms,
                        "annual_risk.parquet")
if (!length(gcms) || any(!file.exists(metric_files)) || any(!file.exists(risk_files))) {
  stop("Expected complete GCM outputs", call. = FALSE)
}

metrics <- rbindlist(lapply(metric_files, function(f) as.data.table(read_parquet(f))))
annual <- rbindlist(lapply(risk_files, function(f) as.data.table(read_parquet(f))))
age85_draw_metrics <- c(
  "temperature_risk_change_85_89_with_cc_pct",
  "temperature_risk_change_85_89_without_cc_pct",
  "cold_risk_change_85_89_with_cc_pct",
  "cold_risk_change_85_89_without_cc_pct",
  "heat_risk_change_85_89_with_cc_pct",
  "heat_risk_change_85_89_without_cc_pct"
)
li65_draw_metrics <- c(
  "LI65plus_change_with_cc_years", "LI65plus_change_without_cc_years",
  "LI65plus_warming_effect_on_change_months"
)
updated_draw_metrics <- c(age85_draw_metrics, li65_draw_metrics)
n_updated_draw_metrics <- sum(updated_draw_metrics %in% names(metrics))
if (n_updated_draw_metrics != 0L && n_updated_draw_metrics != length(updated_draw_metrics)) {
  stop("Uncertainty draws contain only a subset of the updated metrics.", call. = FALSE)
}
has_updated_draw_metrics <- n_updated_draw_metrics == length(updated_draw_metrics)
if (ssp == 3L && !has_updated_draw_metrics) {
  stop(
    "Uncertainty draws lack the updated LI65+ and age-85-89 metrics; ",
    "regenerate them with analysis/run_erf_uncertainty.sh.",
    call. = FALSE
  )
}
metrics[, `:=`(
  table1_all_loss_months = 12 * all_erosion_years,
  table1_heat_loss_months = 12 * heat_erosion_years,
  table1_cold_loss_months = 12 * (all_erosion_years - heat_erosion_years),
  table1_cold_erosion_ratio_pct =
    100 * (all_erosion_years - heat_erosion_years) / g0_years
)]
scenario_status <- if (ssp == 3L) "production" else "exploratory"
scenario_caveat <- if (ssp == 2L) {
  paste(
    "Exploratory SSP2: includes 12 Swiss cities with the mortality-shape",
    "check bypassed; not a production run.",
    if (has_updated_draw_metrics) {
      "Updated age-85-89 and LI65+ draw metrics are included."
    } else {
      paste(
        "Archived draw inputs predate the added age-85-89 and LI65+ metrics;",
        "those estimands are unavailable for SSP2."
      )
    }
  )
} else {
  paste(
    "SSP3 production scenario. Intervals combine 19-GCM spread and ERF",
    "coefficient draws; projected mortality is deterministic."
  )
}
metrics[, `:=`(scenario_status = scenario_status, scenario_caveat = scenario_caveat)]
annual[, `:=`(scenario_status = scenario_status, scenario_caveat = scenario_caveat)]
write_parquet(metrics, file.path(table_dir, "experiment_level.parquet"), compression = "zstd")
write_parquet(annual, file.path(table_dir, "annual_risk_experiment_level.parquet"), compression = "zstd")

id_cols <- c("ssp", "gcm", "draw", "geography")
estimand_cols <- c(
  "lifetime_risk_all_pct", "lifetime_risk_heat_pct", "persistent_share_pct",
  "general_tier_pct", "sustained_tier_pct", "episodic_tier_pct",
  "temperature_death_rate_65plus_end", "climate_death_rate_65plus_end",
  "heat_erosion_ratio_pct", "all_erosion_ratio_pct",
  "g0_years", "g1_years", "heat_erosion_years", "all_erosion_years",
  "table1_all_loss_months", "table1_heat_loss_months",
  "table1_cold_loss_months", "table1_cold_erosion_ratio_pct",
  "LI65plus_change_with_cc_years", "LI65plus_change_without_cc_years",
  "LI65plus_warming_effect_on_change_months",
  "temperature_risk_change_85_89_with_cc_pct",
  "temperature_risk_change_85_89_without_cc_pct",
  "cold_risk_change_85_89_with_cc_pct",
  "cold_risk_change_85_89_without_cc_pct",
  "heat_risk_change_85_89_with_cc_pct",
  "heat_risk_change_85_89_without_cc_pct",
  "delta_LE65_end", "delta_LI65_end"
)
available_estimand_cols <- intersect(estimand_cols, names(metrics))
if (!length(available_estimand_cols)) {
  stop("No supported uncertainty estimands are available.", call. = FALSE)
}
long <- melt(metrics, id.vars = id_cols, measure.vars = available_estimand_cols,
             variable.name = "estimand", value.name = "value")

combined <- long[, .(
  experiment_median = median(value),
  combined_low = quantile(value, 0.025),
  combined_high = quantile(value, 0.975),
  proportion_above_zero = mean(value > 0),
  interval_crosses_zero = quantile(value, 0.025) <= 0 &&
    quantile(value, 0.975) >= 0,
  combined_proportion_above_zero = mean(value > 0),
  combined_interval_crosses_zero = quantile(value, 0.025) <= 0 &&
    quantile(value, 0.975) >= 0,
  n_experiments = .N
), by = .(ssp, geography, estimand)]
erf_by_draw <- long[, .(value = mean(value)),
                    by = .(ssp, geography, estimand, draw)]
erf <- erf_by_draw[, .(
  erf_only_low = quantile(value, 0.025),
  erf_only_high = quantile(value, 0.975),
  erf_only_median = median(value),
  erf_only_proportion_above_zero = mean(value > 0),
  erf_only_n_draws = .N,
  erf_only_interval_crosses_zero = quantile(value, 0.025) <= 0 &&
    quantile(value, 0.975) >= 0
), by = .(ssp, geography, estimand)]
gcm_by_model <- long[, .(value = mean(value)),
                     by = .(ssp, geography, estimand, gcm)]
gcm <- gcm_by_model[, .(
  gcm_low = quantile(value, 0.025),
  gcm_high = quantile(value, 0.975),
  gcm_min = min(value), gcm_max = max(value), gcm_median = median(value),
  gcm_only_low = quantile(value, 0.025),
  gcm_only_high = quantile(value, 0.975),
  gcm_only_proportion_above_zero = mean(value > 0),
  gcm_only_n_models = .N,
  gcm_only_interval_crosses_zero = quantile(value, 0.025) <= 0 &&
    quantile(value, 0.975) >= 0
), by = .(ssp, geography, estimand)]
summary <- Reduce(function(x, y) merge(x, y, by = c("ssp", "geography", "estimand")),
                  list(combined, erf, gcm))
summary[, `:=`(scenario_status = scenario_status, scenario_caveat = scenario_caveat)]

# Published central values for Table 1 metrics.
table1 <- fread(table1_file)[
  ssp == as.integer(Sys.getenv("SSP", "3")) & level %chin% c("europe", "region")
]
central <- rbindlist(list(
  table1[, .(ssp, geography, estimand = "lifetime_risk_all_pct", central = lifetime_risk_pct)],
  table1[, .(ssp, geography, estimand = "lifetime_risk_heat_pct", central = lifetime_risk_heat_pct)],
  table1[, .(ssp, geography, estimand = "persistent_share_pct", central = persistent_share_pct)],
  table1[, .(ssp, geography, estimand = "general_tier_pct", central = general_pct)],
  table1[, .(ssp, geography, estimand = "sustained_tier_pct", central = sustained_pct)],
  table1[, .(ssp, geography, estimand = "episodic_tier_pct", central = episodic_pct)],
  table1[, .(ssp, geography, estimand = "table1_all_loss_months", central = all_loss_months)],
  table1[, .(ssp, geography, estimand = "table1_heat_loss_months",
             central = heat_loss_ensemble_months)],
  table1[, .(ssp, geography, estimand = "table1_cold_loss_months",
             central = cold_loss_ensemble_months)],
  table1[, .(ssp, geography, estimand = "table1_cold_erosion_ratio_pct",
             central = 100 * cold_loss_ensemble_months / g0_months)],
  table1[, .(ssp, geography, estimand = "heat_erosion_ratio_pct", central = heat_ratio_median_pct)],
  table1[, .(ssp, geography, estimand = "all_erosion_ratio_pct", central = all_ratio_pct)]
))

geo_map <- data.table(
  geography = c("Europe", "Eastern Europe", "Northern Europe", "Southern Europe", "Western Europe"),
  level = c("europe", rep("region", 4)),
  id = c("EUROPE", "Eastern", "Northern", "Southern", "Western")
)
levels <- rbindlist(lapply(seq_len(nrow(geo_map)), function(i) {
  f <- file.path(batch_root, "geographies", paste0("ssp", ssp),
                 geo_map$level[i], geo_map$id[i], "04_le_li_levels.csv")
  x <- fread(f)[
    gcm == "ENSEMBLE" & branch %chin% c("with_cc", "without_cc") &
      year %in% c(2020:2024, 2095:2099)
  ]
  x[, period := fifelse(year <= 2024L, "start", "end")]
  x <- x[, .(LE65 = mean(LE65), LI65 = mean(LI65)), by = .(branch, period)]
  get_level <- function(branch_name, period_name, measure) {
    x[branch == branch_name & period == period_name, get(measure)]
  }
  le_change_with <- get_level("with_cc", "end", "LE65") -
    get_level("with_cc", "start", "LE65")
  le_change_without <- get_level("without_cc", "end", "LE65") -
    get_level("without_cc", "start", "LE65")
  li_change_with <- get_level("with_cc", "end", "LI65") -
    get_level("with_cc", "start", "LI65")
  li_change_without <- get_level("without_cc", "end", "LI65") -
    get_level("without_cc", "start", "LI65")
  data.table(
    ssp = ssp, geography = geo_map$geography[i],
    delta_LE65_end = get_level("with_cc", "end", "LE65") -
      get_level("without_cc", "end", "LE65"),
    delta_LI65_end = get_level("with_cc", "end", "LI65") -
      get_level("without_cc", "end", "LI65"),
    LI65plus_change_with_cc_years = li_change_with,
    LI65plus_change_without_cc_years = li_change_without,
    LI65plus_warming_effect_on_change_months =
      12 * (li_change_with - li_change_without),
    LE65_change_with_cc_annual_levels_years = le_change_with,
    LE65_change_without_cc_annual_levels_years = le_change_without
  )
}))
central <- rbind(
  central,
  melt(levels, id.vars = c("ssp", "geography"), variable.name = "estimand",
       value.name = "central"),
  fill = TRUE
)
summary <- merge(summary, central, by = c("ssp", "geography", "estimand"), all.x = TRUE)
summary[, central_source := fifelse(
  is.na(central),
  "experiment median (no matching published central table column)",
  "canonical central pipeline"
)]
summary[is.na(central), central := experiment_median]
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
    low = quantile(value, 0.025), high = quantile(value, 0.975),
    proportion_above_zero = mean(value > 0),
    interval_crosses_zero = quantile(value, 0.025) <= 0 &&
      quantile(value, 0.975) >= 0
  ), by = .(ssp, geography, estimand)][, draws := n]
}))
final <- conv[draws == 500L, .(ssp, geography, estimand, low500 = low, high500 = high)]
conv <- merge(conv, final, by = c("ssp", "geography", "estimand"))
conv[, `:=`(low_difference_from_500 = low - low500,
            high_difference_from_500 = high - high500,
            combined_proportion_above_zero = proportion_above_zero,
            combined_interval_crosses_zero = interval_crosses_zero,
            scenario_status = scenario_status,
            scenario_caveat = scenario_caveat)]
fwrite(conv, file.path(table_dir, "convergence.csv"))

annual_summary <- annual[, .(
  median = median(risk_per_100000),
  low = quantile(risk_per_100000, 0.025),
  high = quantile(risk_per_100000, 0.975),
  proportion_above_zero = mean(risk_per_100000 > 0),
  interval_crosses_zero = quantile(risk_per_100000, 0.025) <= 0 &&
    quantile(risk_per_100000, 0.975) >= 0,
  combined_proportion_above_zero = mean(risk_per_100000 > 0),
  combined_interval_crosses_zero = quantile(risk_per_100000, 0.025) <= 0 &&
    quantile(risk_per_100000, 0.975) >= 0,
  sign_consistency = mean(sign(risk_per_100000) == sign(median(risk_per_100000)))
), by = .(ssp, geography, year, age_scope, estimand)]
annual_summary[, `:=`(scenario_status = scenario_status, scenario_caveat = scenario_caveat)]
write_parquet(annual_summary, file.path(table_dir, "annual_risk_summary.parquet"),
              compression = "zstd")

checks <- data.table(
  check = c("gcms", "draws", "geographies", "metric_experiments",
            "annual_rows", "headline_estimands", "separate_interval_flags",
            "interval_proportions", "age85_risk_estimands", "LI65_change_estimands",
            "finite_metrics", "finite_annual", "scenario_caveat"),
  value = c(
    as.character(uniqueN(metrics$gcm)), as.character(uniqueN(metrics$draw)),
    as.character(uniqueN(metrics$geography)), as.character(nrow(metrics)),
    as.character(nrow(annual)),
    as.character(uniqueN(summary$estimand)),
    as.character(all(vapply(
      summary[, .(
        gcm_only_interval_crosses_zero, erf_only_interval_crosses_zero,
        combined_interval_crosses_zero
      )],
      function(x) !anyNA(x), logical(1)
    ))),
    as.character(all(vapply(
      summary[, .(
        gcm_only_proportion_above_zero, erf_only_proportion_above_zero,
        combined_proportion_above_zero
      )],
      function(x) all(is.finite(x) & x >= 0 & x <= 1), logical(1)
    ))),
    as.character(all(age85_draw_metrics %in% summary$estimand)),
    as.character(all(li65_draw_metrics %in% summary$estimand)),
    if (all(vapply(metrics, function(x) !is.numeric(x) || all(is.finite(x)), logical(1))))
      "TRUE" else "FALSE",
    if (all(is.finite(annual$risk_per_100000))) "TRUE" else "FALSE",
    scenario_status
  ),
  expected = c(
    as.character(length(gcms)), "500", as.character(uniqueN(metrics$geography)),
    as.character(length(gcms) * 500L * uniqueN(metrics$geography)),
    as.character(length(gcms) * 500L * uniqueN(metrics$geography) * 80L * 2L * 2L),
    as.character(length(available_estimand_cols)), "TRUE", "TRUE",
    as.character(has_updated_draw_metrics), as.character(has_updated_draw_metrics),
    "TRUE", "TRUE", scenario_status
  )
)
checks[, scenario_caveat := scenario_caveat]
fwrite(checks, file.path(table_dir, "checks.csv"))
if (any(checks$value != checks$expected)) {
  stop("Uncertainty summary checks failed; see checks.csv", call. = FALSE)
}
message("Wrote uncertainty summaries to ", table_dir)
