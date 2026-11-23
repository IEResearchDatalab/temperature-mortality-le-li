#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
source("analysis/decade_decomposition.R")
source("analysis/output_names.R")

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
figure_root <- Sys.getenv("FIGURE_ROOT", "results/figures_simon_20261123")
table1_file <- Sys.getenv(
  "TABLE1_SSP3",
  "results/paper_figures/20261016_ssp1_ssp2_ssp3_ordered/table1/table1_ssp3_geographies.csv"
)
uncertainty_summary_file <- Sys.getenv(
  "UNCERTAINTY_SUMMARY",
  "results/europe/uncertainty_20261110/tables/ssp3/summary.csv"
)
output_date <- Sys.getenv("TABLE_DATE", "20261123")
if (!grepl("^[0-9]{8}$", output_date)) {
  stop("TABLE_DATE must use YYYYMMDD.", call. = FALSE)
}

data_dir <- file.path(figure_root, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

geographies <- data.table(
  level = c("europe", rep("region", 4L)),
  id = c("EUROPE", "Eastern", "Northern", "Southern", "Western"),
  geography = c("Europe", "Eastern Europe", "Northern Europe", "Southern Europe", "Western Europe")
)
branch_names <- c(with_cc = "with climate change", without_cc = "without climate change")
start_years <- 2020:2024
end_years <- 2095:2099
start_period <- "2020-2024"
end_period <- "2095-2099"
cold_causes <- c("ExtrCold", "ModCold")
heat_causes <- c("ModHeat", "NormExtrHeat", "ExtrHeat")
temperature_causes <- c(cold_causes, heat_causes)
cause_levels <- decade_cause_levels

check_rows <- list()
add_check <- function(name, geography, ssp, observed = NA_real_, expected = NA_real_,
                      tolerance = NA_real_, status, note = "") {
  check_rows[[length(check_rows) + 1L]] <<- data.table(
    check = name, geography = geography, ssp = ssp,
    observed = observed, expected = expected, difference = observed - expected,
    tolerance = tolerance, status = status, note = note
  )
}

check_numeric <- function(name, geography, ssp, observed, expected, tolerance,
                          note = "") {
  observed_value <- if (length(observed) == 1L) as.numeric(observed) else NA_real_
  expected_value <- if (length(expected) == 1L) as.numeric(expected) else NA_real_
  status <- if (is.finite(observed_value) && is.finite(expected_value) &&
                abs(observed_value - expected_value) <= tolerance) "PASS" else "FAIL"
  add_check(name, geography, ssp, observed_value, expected_value, tolerance, status, note)
}

period_mean_schedule <- function(ssp, level, id) {
  path <- file.path(
    batch_root, "geographies", paste0("ssp", ssp), level, id, "03_master_table.csv"
  )
  if (!file.exists(path)) stop("Missing pooled master table: ", path, call. = FALSE)
  x <- fread(path, select = c(
    "gcm", "branch", "year", "age", "range", "deaths_component", "pop"
  ))
  x <- x[
    gcm == "ENSEMBLE" & age %in% decade_age_levels &
      year %in% c(start_years, end_years) & range %chin% cause_levels
  ]
  if (!nrow(x) || any(!is.finite(x$deaths_component)) ||
      any(!is.finite(x$pop)) || any(x$pop <= 0)) {
    stop("Invalid pooled ENSEMBLE master rows in ", path, call. = FALSE)
  }
  x[, period := fifelse(year %in% start_years, start_period, end_period)]
  x[, mx_cause := deaths_component / pop]
  schedule <- x[, .(mx_cause = mean(mx_cause)),
                by = .(branch, period, age, cause = range)]
  expected <- 2L * 2L * length(decade_age_levels) * length(cause_levels)
  if (nrow(schedule) != expected ||
      anyDuplicated(schedule, by = c("branch", "period", "age", "cause")) ||
      any(!is.finite(schedule$mx_cause)) || any(schedule$mx_cause < 0)) {
    stop("Period-mean schedule is incomplete or invalid for SSP", ssp, " ",
         level, "/", id, ".", call. = FALSE)
  }
  schedule
}

period_le <- function(schedule) {
  schedule[, .(
    LE65 = decade_life_expectancy(decade_cause_vector(.SD))
  ), by = .(branch, period)]
}

life_expectancy_without_heat <- function(schedule) {
  schedule[cause %chin% c(cold_causes, "rest"),
           .(mx = sum(mx_cause)), by = age][order(age), decade_life_expectancy(mx)]
}

read_level_periods <- function(ssp, level, id) {
  path <- file.path(
    batch_root, "geographies", paste0("ssp", ssp), level, id, "04_le_li_levels.csv"
  )
  if (!file.exists(path)) stop("Missing pooled LE65 levels: ", path, call. = FALSE)
  x <- fread(path, select = c("gcm", "branch", "year", "LE65"))
  x <- x[
    gcm == "ENSEMBLE" & branch %chin% names(branch_names) &
      year %in% c(start_years, end_years)
  ]
  x[, period := fifelse(year %in% start_years, start_period, end_period)]
  result <- x[, .(
    LE65 = mean(LE65),
    years = uniqueN(year)
  ), by = .(branch, period)]
  if (nrow(result) != 4L || any(result$years != 5L) ||
      any(!is.finite(result$LE65))) {
    stop("Pooled ENSEMBLE LE65 levels are incomplete in ", path, call. = FALSE)
  }
  result
}

read_temperature_decomposition <- function(ssp, level, id) {
  path <- file.path(
    batch_root, "geographies", paste0("ssp", ssp), level, id,
    "04_within_branch_period_decomposition.csv"
  )
  if (!file.exists(path)) stop("Missing within-branch decomposition: ", path, call. = FALSE)
  x <- fread(path)
  x <- x[gcm == "ENSEMBLE" & branch == "with_cc"]
  x[, `:=`(
    from_start = as.integer(substr(period_from, 1L, 4L)),
    to_start = as.integer(substr(period_to, 1L, 4L))
  )]
  x <- x[from_start >= 2020L & to_start <= 2095L]
  expected <- 15L * length(decade_age_levels) * length(cause_levels)
  if (nrow(x) != expected ||
      anyDuplicated(x, by = c("period_from", "period_to", "age", "cause")) ||
      any(!is.finite(x$le_contribution))) {
    stop("Within-branch ENSEMBLE decomposition does not cover all 5-year steps in SSP",
         ssp, " ", level, "/", id, ".", call. = FALSE)
  }
  x
}

figure_panel_path <- function(name, id) {
  figure_file(figure_root, paste0("figure02_le65_contributions_", name, "_data.csv"), id)
}

calculate_one <- function(ssp, geo) {
  level <- geo$level
  id <- geo$id
  label <- geo$geography
  schedule <- period_mean_schedule(ssp, level, id)
  period_levels <- period_le(schedule)
  get_period_le <- function(branch, period) {
    branch_value <- branch
    period_value <- period
    period_levels$LE65[
      period_levels$branch == branch_value & period_levels$period == period_value
    ]
  }
  g0 <- get_period_le("without_cc", end_period) -
    get_period_le("without_cc", start_period)
  g1 <- get_period_le("with_cc", end_period) -
    get_period_le("with_cc", start_period)
  if (length(g0) != 1L || length(g1) != 1L || g0 == 0) {
    stop("Invalid period-mean LE65 gains for SSP", ssp, " ", label, ".", call. = FALSE)
  }

  level_data <- read_level_periods(ssp, level, id)
  level_le <- function(branch, period) {
    branch_value <- branch
    period_value <- period
    level_data$LE65[
      level_data$branch == branch_value & level_data$period == period_value
    ]
  }
  gain_level_g0 <- level_le("without_cc", end_period) -
    level_le("without_cc", start_period)
  gain_level_g1 <- level_le("with_cc", end_period) -
    level_le("with_cc", start_period)

  noheat_gain <- schedule[, .(
    gain = life_expectancy_without_heat(.SD[period == end_period]) -
      life_expectancy_without_heat(.SD[period == start_period])
  ), by = branch]
  g0_noheat <- noheat_gain[branch == "without_cc", gain]
  g1_noheat <- noheat_gain[branch == "with_cc", gain]
  all_loss <- g0 - g1
  cold_only_loss <- g0_noheat - g1_noheat
  heat_loss <- all_loss - cold_only_loss

  decomposition <- read_temperature_decomposition(ssp, level, id)
  cause_months <- decomposition[, .(
    months = 12 * sum(le_contribution)
  ), by = cause]
  cold_months <- cause_months[cause %chin% cold_causes, sum(months)]
  heat_months <- cause_months[cause %chin% heat_causes, sum(months)]
  together_months <- cause_months[cause %chin% temperature_causes, sum(months)]
  rest_months <- cause_months[cause == "rest", sum(months)]
  decomp_gain <- sum(decomposition$le_contribution)

  add_check(
    "period-mean decomposition closes to G1", label, ssp,
    decomp_gain, g1, 1e-6,
    if (abs(decomp_gain - g1) <= 1e-6) "PASS" else "FAIL",
    "Sum of all causes over the 15 consecutive 5-year steps."
  )
  check_numeric(
    "cold + heat equals five-range diamond", label, ssp,
    cold_months + heat_months, together_months, 1e-8,
    "Cold is extreme + moderate; heat is moderate + normal-extreme + extreme heat."
  )
  check_numeric(
    "temperature + other-causes contribution equals G1",
    label, ssp, (together_months + rest_months) / 12, g1, 1e-6,
    "Other causes are the Part 04 rest contribution."
  )
  add_check(
    "annual-level gain versus period-mean gain", label, ssp,
    gain_level_g1, g1, NA_real_,
    "INFO",
    sprintf("04_le_li_levels mean-of-annual-LE gain differs from period-mean-rate G1 by %.9f years.",
            gain_level_g1 - g1)
  )
  add_check(
    "annual-level G0 versus period-mean G0", label, ssp,
    gain_level_g0, g0, NA_real_,
    "INFO",
    sprintf("04_le_li_levels mean-of-annual-LE gain differs from period-mean-rate G0 by %.9f years.",
            gain_level_g0 - g0)
  )

  list(
    long = data.table(
      ssp = ssp, scenario = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")[ssp],
      geography = label, metric = c(
        "LE65_2020_24_with_climate_change",
        "LE65_2095_99_with_climate_change",
        "LE65_2020_24_without_climate_change",
        "LE65_2095_99_without_climate_change",
        "gain_with_climate_change",
        "gain_without_climate_change",
        "gain_with_climate_change_from_annual_levels",
        "gain_without_climate_change_from_annual_levels",
        "cold_contribution", "heat_contribution", "cold_and_heat_together_contribution",
        "cold_contribution_pct_of_gain", "heat_contribution_pct_of_gain",
        "cold_and_heat_together_pct_of_gain",
        "reduction_due_to_warming_all_temperatures",
        "reduction_due_to_warming_through_heat_alone"
      ),
      value = c(
        level_le("with_cc", start_period), level_le("with_cc", end_period),
        level_le("without_cc", start_period), level_le("without_cc", end_period),
        g1, g0, gain_level_g1, gain_level_g0,
        cold_months, heat_months, together_months,
        100 * cold_months / (12 * g1),
        100 * heat_months / (12 * g1),
        100 * together_months / (12 * g1),
        100 * all_loss / g0,
        100 * heat_loss / g0
      ),
      units = c(rep("years", 8L), rep("months", 3L), rep("%", 5L)),
      definition = c(
        "Mean annual LE65 in 2020-24 from pooled 04_le_li_levels.csv; ENSEMBLE.",
        "Mean annual LE65 in 2095-99 from pooled 04_le_li_levels.csv; ENSEMBLE.",
        "Mean annual LE65 in 2020-24 from pooled 04_le_li_levels.csv; ENSEMBLE.",
        "Mean annual LE65 in 2095-99 from pooled 04_le_li_levels.csv; ENSEMBLE.",
        "Period-mean-rate LE65(2095-99) minus LE65(2020-24), with climate change.",
        "Period-mean-rate LE65(2095-99) minus LE65(2020-24), without climate change; G0.",
        "Difference of the period means of annual LE65 values read from 04_le_li_levels.csv.",
        "Difference of the period means of annual LE65 values read from 04_le_li_levels.csv.",
        "Five-year-step Horiuchi contributions, extreme + moderate cold, summed over ages.",
        "Five-year-step Horiuchi contributions, moderate + normal-extreme + extreme heat, summed over ages.",
        "Five-year-step Horiuchi contributions across all five temperature ranges, summed over ages.",
        "Cold contribution / period-mean-rate G1.",
        "Heat contribution / period-mean-rate G1.",
        "All five temperature contributions / period-mean-rate G1.",
        "100 * (G0 - G1) / G0, using the existing period-mean-rate Table 1 baseline.",
        "100 * heat loss / G0; heat loss follows the heat-deleted versus full schedule definition in erf_uncertainty_pilot.R."
      )
    ),
    wide = data.table(
      ssp = ssp, scenario = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")[ssp],
      geography = label, gain_in_years = g1,
      cold_pct = 100 * cold_months / (12 * g1),
      heat_pct = 100 * heat_months / (12 * g1),
      cold_and_heat_together_pct = 100 * together_months / (12 * g1),
      reduction_due_to_warming_pct = 100 * all_loss / g0,
      note = "SSP1 uses exploratory 0.9995 replacement; SSP2 includes 12 Swiss cities with mortality-shape check bypassed. Neither is a production run."
    ),
    sx = data.table(
      ssp = ssp, scenario = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")[ssp],
      geography = label, gain_in_years = g1,
      cold_months = cold_months, heat_months = heat_months,
      cold_and_heat_together_months = together_months,
      reduction_due_to_warming_all_temperatures_pct = 100 * all_loss / g0,
      reduction_due_to_warming_through_heat_alone_pct = 100 * heat_loss / g0,
      note = "SSP1 uses exploratory 0.9995 replacement; SSP2 includes 12 Swiss cities with mortality-shape check bypassed. Neither is a production run."
    ),
    comparison = list(
      ssp = ssp, level = level, id = id, geography = label,
      g0 = g0, g1 = g1, g0_level = gain_level_g0, g1_level = gain_level_g1,
      cold_months = cold_months, heat_months = heat_months,
      together_months = together_months, rest_months = rest_months,
      all_ratio = 100 * all_loss / g0, heat_ratio = 100 * heat_loss / g0
    ),
    cause_contributions = cause_months[, .(
      ssp = ssp, id = id, geography = label, cause, months
    )]
  )
}

results <- list()
for (ssp in 1:3) {
  for (i in seq_len(nrow(geographies))) {
    results[[length(results) + 1L]] <- calculate_one(ssp, geographies[i])
  }
}
long <- rbindlist(lapply(results, `[[`, "long"), use.names = TRUE)
wide <- rbindlist(lapply(results, `[[`, "wide"), use.names = TRUE)
sx <- rbindlist(lapply(results, `[[`, "sx"), use.names = TRUE)
comparisons <- rbindlist(lapply(results, `[[`, "comparison"), use.names = TRUE)
cause_contributions <- rbindlist(
  lapply(results, `[[`, "cause_contributions"), use.names = TRUE
)

setorder(long, ssp, geography, metric)
wide[, geography_order := match(geography, geographies$geography)]
sx[, geography_order := match(geography, geographies$geography)]
setorder(wide, ssp, geography_order)
setorder(sx, ssp, geography_order)
wide[, geography_order := NULL]
sx[, geography_order := NULL]

add_check(
  "long table row count and key uniqueness", "all geographies", NA_integer_,
  observed = nrow(long), expected = 240, tolerance = 0,
  status = if (nrow(long) == 240L &&
               uniqueN(long, by = c("ssp", "geography", "metric")) == 240L) {
    "PASS"
  } else {
    "FAIL"
  },
  note = "16 metrics for each of 3 scenarios and 5 geographies."
)
add_check(
  "wide table row count and key uniqueness", "all geographies", NA_integer_,
  observed = nrow(wide), expected = 15, tolerance = 0,
  status = if (nrow(wide) == 15L &&
               uniqueN(wide, by = c("ssp", "geography")) == 15L) "PASS" else "FAIL",
  note = "One row for each scenario and geography."
)
add_check(
  "Table Sx row count and key uniqueness", "all geographies", NA_integer_,
  observed = nrow(sx), expected = 15, tolerance = 0,
  status = if (nrow(sx) == 15L &&
               uniqueN(sx, by = c("ssp", "geography")) == 15L) "PASS" else "FAIL",
  note = "One row for each scenario and geography."
)
add_check(
  "finite long-table values", "all geographies", NA_integer_,
  observed = sum(!is.finite(long$value)), expected = 0, tolerance = 0,
  status = if (all(is.finite(long$value))) "PASS" else "FAIL",
  note = "All reported metric values must be finite."
)

draft_compare <- function(geography, ssp, metric, draft, tolerance) {
  geography_value <- geography
  ssp_value <- ssp
  row <- wide[wide$geography == geography_value & wide$ssp == ssp_value]
  observed <- switch(
    metric,
    gain = row$gain_in_years,
    cold = row$cold_pct,
    heat = row$heat_pct,
    together = row$cold_and_heat_together_pct,
    all_ratio = row$reduction_due_to_warming_pct
  )
  status <- if (length(observed) == 1L && abs(observed - draft) <= tolerance) {
    "MATCH_ROUNDED"
  } else {
    "DIFF"
  }
  add_check(
    paste0("Simon draft: ", metric), geography, ssp,
    if (length(observed)) observed else NA_real_, draft, tolerance, status,
    "Provisional draft figure; difference is observed minus draft."
  )
}

for (ssp in 1:3) {
  for (metric_entry in list(
    list("gain", c(18, 10.5, 4.5), c(0.55, 0.05, 0.05)),
    list("cold", c(7, 7, 9), rep(0.55, 3L)),
    list("heat", c(2, 1, -5), rep(0.55, 3L)),
    list("together", c(9, 8, 4), rep(0.55, 3L)),
    list("all_ratio", c(0.1, 0.4, 4.3), rep(0.05, 3L))
  )) {
    draft_compare("Europe", ssp, metric_entry[[1]], metric_entry[[2]][ssp],
                  metric_entry[[3]][ssp])
  }
  for (metric_entry in list(
    list("gain", c(18.5, 10.5, 4), c(0.05, 0.05, 0.55)),
    list("cold", c(6, 6, 8), rep(0.55, 3L)),
    list("heat", c(4, 2, -7), rep(0.55, 3L)),
    list("together", c(9, 8, 1), rep(0.55, 3L)),
    list("all_ratio", c(0.3, 1.2, 8.9), rep(0.05, 3L))
  )) {
    draft_compare("Southern Europe", ssp, metric_entry[[1]], metric_entry[[2]][ssp],
                  metric_entry[[3]][ssp])
  }
  for (region in c("Northern Europe", "Western Europe", "Eastern Europe")) {
    reduction <- switch(
      region,
      "Northern Europe" = c(-0.1, -0.3, -0.1),
      "Western Europe" = c(0.0, 0.3, 3.4),
      "Eastern Europe" = c(0.0, 0.4, 4.8)
    )
    draft_compare(region, ssp, "all_ratio", reduction[ssp], 0.05)
  }
}

europe_ssp3 <- wide[geography == "Europe" & ssp == 3L]
check_numeric(
  "Simon draft SSP3 Europe cold months", "Europe", 3L,
  sx[geography == "Europe" & ssp == 3L, cold_months], 4.6, 0.05,
  "Draft results-text value."
)
check_numeric(
  "Simon draft SSP3 Europe heat months", "Europe", 3L,
  sx[geography == "Europe" & ssp == 3L, heat_months], -2.4, 0.05,
  "Draft results-text value."
)

if (!file.exists(uncertainty_summary_file)) {
  add_check(
    "Earlier uncertainty-run comparison", "Europe", 3L,
    status = "FAIL", note = paste("Missing file:", uncertainty_summary_file)
  )
} else {
  uncertainty_summary <- fread(uncertainty_summary_file)
  uncertainty_values <- uncertainty_summary[
    geography == "Europe" &
      estimand %chin% c("all_erosion_ratio_pct", "heat_erosion_ratio_pct"),
    .(
      metric = estimand,
      run_median_pct = experiment_median,
      canonical_central_pct = central,
      central_source
    )
  ]
  current_values <- data.table(
    metric = c("all_erosion_ratio_pct", "heat_erosion_ratio_pct"),
    central_pct = c(
      europe_ssp3$reduction_due_to_warming_pct,
      sx[geography == "Europe" & ssp == 3L,
         reduction_due_to_warming_through_heat_alone_pct]
    )
  )
  uncertainty_comparison <- merge(uncertainty_values, current_values, by = "metric")
  for (i in seq_len(nrow(uncertainty_comparison))) {
    row <- uncertainty_comparison[i]
    add_check(
      paste0("Earlier uncertainty-run median: ", row$metric),
      "Europe", 3L, row$central_pct, row$run_median_pct, NA_real_, "INFO",
      sprintf("Central estimate minus 19-GCM x 500-draw experiment median = %.4f percentage points; source: %s.",
              row$central_pct - row$run_median_pct, uncertainty_summary_file)
    )
    add_check(
      paste0("Earlier canonical central comparison: ", row$metric),
      "Europe", 3L, row$central_pct, row$canonical_central_pct,
      NA_real_, "INFO",
      paste("The summary's central value uses:", row$central_source,
            "Source:", uncertainty_summary_file)
    )
  }
}

if (!file.exists(table1_file)) {
  add_check(
    "SSP3 Table 1 cross-check", "Europe and four regions", 3L,
    status = "FAIL", note = paste("Missing file:", table1_file)
  )
} else {
  table1 <- fread(table1_file)
  expected_labels <- c("Europe", "Eastern Europe", "Northern Europe",
                       "Southern Europe", "Western Europe")
  table1 <- table1[geography %chin% expected_labels]
  for (label in expected_labels) {
    current <- sx[geography == label & ssp == 3L]
    prior <- table1[geography == label]
    if (nrow(current) != 1L || nrow(prior) != 1L) {
      add_check(
        "SSP3 Table 1 row coverage", label, 3L,
        observed = nrow(prior), expected = 1, tolerance = 0,
        status = "FAIL", note = table1_file
      )
    } else {
      check_numeric(
        "SSP3 Table 1 G0 cross-check", label, 3L,
        comparisons[
          geography == label & ssp == 3L, g0
        ] * 12,
        prior$le65_gain_without_cc_months, 1e-5, table1_file
      )
      check_numeric(
        "SSP3 Table 1 all-temperature erosion cross-check", label, 3L,
        current$reduction_due_to_warming_all_temperatures_pct,
        100 * prior$all_temperature_erosion_ratio, 0.005, table1_file
      )
      add_check(
        "Archived Table 1 heat-erosion comparison", label, 3L,
        current$reduction_due_to_warming_through_heat_alone_pct,
        100 * prior$heat_only_erosion_ratio, NA_real_, "INFO",
        paste(
          "The archived Table 1 column uses its earlier heat-erosion calculation;",
          "Table Sx applies the heat-deleted pooled central-schedule definition.",
          table1_file
        )
      )
    }
  }
}

for (i in seq_len(nrow(geographies))) {
  geo <- geographies[i]
  geo_id <- geo$id
  for (ssp in 1:3) {
    ssp_value <- ssp
    comparison <- comparisons[
      comparisons$ssp == ssp_value & comparisons$id == geo_id
    ]
    expected_by_cause <- cause_contributions[
      cause_contributions$ssp == ssp_value & cause_contributions$id == geo_id &
        cause_contributions$cause %chin% temperature_causes
    ]
    figure_values <- NULL
    panel_missing <- FALSE
    panel_a_path <- figure_panel_path("panel_a", geo$id)
    decade_path <- file.path(
      data_dir, sprintf("decade_decomposition_%s_%s.csv", geo$level, geo$id)
    )
    if (!file.exists(panel_a_path) || !file.exists(decade_path)) {
      add_check(
        "Figure 2 panel a decade-data cross-check", geo$geography, ssp,
        status = "FAIL",
        note = paste("Missing:", c(panel_a_path, decade_path)[
          !file.exists(c(panel_a_path, decade_path))
        ], collapse = "; ")
      )
    } else {
      panel_a <- fread(panel_a_path)[
        ssp_label == c("SSP1–2.6", "SSP2–4.5", "SSP3–7.0")[ssp]
      ]
      panel_a[, cause := as.character(cause)]
      panel_a[, decade_label := as.character(decade_label)]
      panel_a <- panel_a[, .(value = sum(value)), by = .(decade_label, cause)]
      decade_source <- fread(decade_path)
      decade_values <- decade_source[
        decade_source$ssp == ssp_value &
          decade_source$cause %chin% temperature_causes,
        .(
          value = sum(le_contribution_months),
          decade_label = sprintf("%d–%02d", decade_to, (decade_to + 9L) %% 100L)
        ),
        by = .(decade_to, cause)
      ][, .(value = sum(value)), by = .(decade_label, cause)]
      decade_contributions <- decade_source[
        decade_source$ssp == ssp_value &
          decade_source$cause %chin% temperature_causes,
        .(months = sum(le_contribution_months)), by = cause
      ]
      panel_a_check <- merge(
        panel_a, decade_values, by = c("decade_label", "cause"),
        suffixes = c("_figure", "_decomposition"), all = TRUE
      )
      panel_a_gap <- max(abs(
        panel_a_check$value_figure - panel_a_check$value_decomposition
      ))
      check_numeric(
        "Figure 2 panel a decade-data cross-check", geo$geography, ssp,
        panel_a_gap, 0, 1e-12,
        "Panel a bars equal the age-summed decade Horiuchi contributions."
      )
      decade_cold <- decade_contributions[cause %chin% cold_causes, sum(months)]
      decade_heat <- decade_contributions[cause %chin% heat_causes, sum(months)]
      decade_together <- decade_cold + decade_heat
      for (metric_entry in list(
        list("cold", comparison$cold_months, decade_cold),
        list("heat", comparison$heat_months, decade_heat),
        list("cold and heat together", comparison$together_months, decade_together)
      )) {
        add_check(
          paste0("Five-year versus decade contribution: ", metric_entry[[1L]]),
          geo$geography, ssp, metric_entry[[2L]], metric_entry[[3L]],
          NA_real_, "INFO",
          "The table uses all consecutive 5-year transitions; panel a recomputes contributions over 10-year schedules."
        )
      }
    }
    for (panel in c("b", "c")) {
      path <- figure_panel_path(paste0("panel_", panel), geo$id)
      if (!file.exists(path)) {
        add_check(
          paste0("Figure 2 panel ", panel, " cross-check"), geo$geography,
          ssp, status = "FAIL", note = paste("Missing file:", path)
        )
        panel_missing <- TRUE
        next
      }
      plotted <- fread(path)[      ssp_label == c("SSP1–2.6", "SSP2–4.5", "SSP3–7.0")[ssp]]
      plotted[, cause := as.character(cause)]
      factor <- if (panel == "b") 3.5 else 4
      plotted <- plotted[, .(months = sum(value) * factor), by = cause]
      if (panel == "b") {
        figure_values <- plotted
      } else {
        figure_values <- rbindlist(list(figure_values, plotted))[
          , .(months = sum(months)), by = cause
        ]
      }
    }
    if (!panel_missing) {
      for (cause in temperature_causes) {
        cause_value <- cause
        observed <- figure_values$months[figure_values$cause == cause_value]
        expected <- expected_by_cause$months[expected_by_cause$cause == cause_value]
        check_numeric(
          paste0("Figure 2 panels b+c ", cause, " cross-check"),
          geo$geography, ssp, observed, expected, 1e-7,
          "Panel bars rescaled by their 35- or 40-year span to total months."
        )
      }
    }
  }
}

check_rows[[length(check_rows) + 1L]] <- data.table(
  check = "SSP1 sensitivity provenance", geography = NA_character_, ssp = 1L,
  observed = NA_real_, expected = NA_real_, difference = NA_real_, tolerance = NA_real_,
  status = "INFO",
  note = paste(
    "The 0.9995 survival-ratio replacement is an exploratory sensitivity, not a production run.",
    "See SSP1 run manifests under", batch_root
  )
)
check_rows[[length(check_rows) + 1L]] <- data.table(
  check = "SSP2 sensitivity provenance", geography = NA_character_, ssp = 2L,
  observed = NA_real_, expected = NA_real_, difference = NA_real_, tolerance = NA_real_,
  status = "INFO",
  note = paste(
    "Includes 12 Swiss cities with the mortality-shape check bypassed; not a production run.",
    "See SSP2 run manifests under", batch_root
  )
)
for (ssp in 1:2) {
  manifest <- file.path(batch_root, paste0("ssp", ssp), "run_manifest.txt")
  if (!file.exists(manifest)) {
    add_check(
      "Sensitivity run manifest", NA_character_, ssp,
      status = "FAIL", note = paste("Missing file:", manifest)
    )
    next
  }
  settings <- readLines(manifest, warn = FALSE)
  replacement <- sub("^assr_one_replacement=", "", grep(
    "^assr_one_replacement=", settings, value = TRUE
  ))
  bypass <- sub("^skip_mx_plausibility_check=", "", grep(
    "^skip_mx_plausibility_check=", settings, value = TRUE
  ))
  expected_replacement <- if (ssp == 1L) "0.9995" else "none"
  provenance_ok <- length(replacement) == 1L &&
    identical(replacement, expected_replacement) &&
    length(bypass) == 1L && identical(bypass, "1")
  add_check(
    "Sensitivity run manifest", NA_character_, ssp,
    observed = as.numeric(provenance_ok), expected = 1, tolerance = 0,
    status = if (provenance_ok) "PASS" else "FAIL",
    note = paste(manifest, "replacement=", paste(replacement, collapse = ","),
                 "skip_mx_plausibility_check=", paste(bypass, collapse = ","))
  )
}

checks <- rbindlist(check_rows, use.names = TRUE, fill = TRUE)
fwrite(long, file.path(data_dir, paste0("table02_le65_gain_long_", output_date, ".csv")))
fwrite(wide, file.path(data_dir, "table02_le65_gain_wide.csv"))
fwrite(sx, file.path(data_dir, "tableSx_le65_gain_months.csv"))
fwrite(checks, file.path(data_dir, "table02_checks.csv"))

commit <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)
commit <- if (length(commit) && !is.na(commit[1L])) commit[1L] else "unavailable"
uncertainty_note <- if (exists("uncertainty_comparison")) {
  all <- uncertainty_comparison[metric == "all_erosion_ratio_pct"]
  heat <- uncertainty_comparison[metric == "heat_erosion_ratio_pct"]
  sprintf(
    paste(
      "The 2026-11-10 uncertainty summary reports experiment medians of",
      "%.4f%% (all temperatures) and %.4f%% (heat). These are 19-GCM by",
      "500-draw medians; the pooled central estimates here are %.4f%% and",
      "%.4f%%, respectively."
    ),
    all$run_median_pct, heat$run_median_pct, all$central_pct, heat$central_pct
  )
} else {
  "No prior uncertainty summary was available for comparison."
}
readme <- c(
  "# Table 2 and Table Sx",
  "",
  paste("Analysis commit:", commit),
  paste("Input batch:", batch_root),
  "",
  paste(
    "G1 is the with-climate-change LE65 gain from the 2020-24 to 2095-99",
    "period-mean mortality schedules; G0 is the same gain without climate",
    "change. The 2020-24 baseline for both is the pooled ENSEMBLE, age 65-100",
    "cause-specific mortality-rate schedule averaged over those five years,",
    "matching `compute_period_metrics()` in `analysis/paper_table1.R`."
  ),
  "",
  paste(
    "The four endpoint LE65 levels are the means of annual ENSEMBLE values in",
    "`04_le_li_levels.csv`. Their differences are also included in the long",
    "table, separately labelled; they need not equal period-mean-rate gains.",
    "Table columns and reductions use the period-mean-rate G0/G1 so they remain",
    "consistent with Table 1 and the five-year-step decomposition."
  ),
  "",
  paste(
    "Temperature contributions are sums of the with-CC, ENSEMBLE Horiuchi",
    "contributions across ages 65-100 and all 15 consecutive five-year",
    "transitions. Cold is extreme plus moderate cold; heat is moderate,",
    "normal-extreme and extreme heat; together is the sum of all five",
    "temperature ranges. The diamond check verifies that cold plus heat equals",
    "this five-range total. Other causes are `rest`."
  ),
  "",
  paste(
    "All-temperature reduction is `100 * (G0 - G1) / G0`, matching the",
    "erosion-ratio definition used in methods section 2.8 and",
    "`analysis/paper_table1.R`. Heat-only reduction follows the heat-deleted",
    "versus full schedule calculation in `analysis/erf_uncertainty_pilot.R`:",
    "subtract the no-heat cold-only gain difference from the all-temperature",
    "gain loss, then divide by the same G0."
  ),
  "",
  uncertainty_note,
  "",
  paste(
    "The SSP1 0.9995 survival-ratio replacement and SSP2 inclusion of 12",
    "Swiss cities with the mortality-shape check bypassed are exploratory",
    "sensitivities, not production runs. They are retained as caveats in this",
    "table output."
  )
)
writeLines(readme, file.path(data_dir, "table02_README.md"), useBytes = TRUE)

failed <- checks[status == "FAIL"]
message("Wrote Table 2 outputs under ", data_dir)
message("Checks: ", checks[, sum(status == "PASS")], " PASS; ",
        checks[, sum(status == "FAIL")], " FAIL; ",
        checks[, sum(status %chin% c("INFO", "DIFF", "MATCH_ROUNDED"))], " informational/comparison rows.")
if (nrow(failed)) {
  stop(nrow(failed), " Table 2 checks failed; see table02_checks.csv.", call. = FALSE)
}
