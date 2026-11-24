#!/usr/bin/env Rscript

################################################################################
# Table 2 and Table Sx: gains in LE65 and temperature contributions.
# Requires: Parts 00-04 for all three SSPs, pooled Europe and four regions;
# each <batch>/geographies/ssp<SSP>/<level>/<id>/{03_master_table,
# 04_le_li_levels,04_within_branch_period_decomposition}.csv; run manifests
# <batch>/ssp<SSP>/run_manifest.txt; Figure 2 panel and decade CSVs in
# FIGURE_ROOT/data/; TABLE1_SSP3 and UNCERTAINTY_SUMMARY comparison CSVs;
# analysis/decade_decomposition.R and analysis/output_names.R.
# Produces: FIGURE_ROOT/data/table02_* and tableSx_le65_gain_months.csv,
# plus LE65/LI65+ level, change and LI65+ contribution exports.
# How to run: BATCH_ROOT=results/europe/heatsplit_20261013 \
#   FIGURE_ROOT=results/paper_output_new \
#   TABLE1_SSP3=results/paper_output_new/table1/table1_ssp3_geographies.csv \
#   UNCERTAINTY_SUMMARY=results/paper_output_new/uncertainty/tables/ssp3/summary.csv \
#   Rscript analysis/paper_table2.R
# Notes: SSP1 replacement and SSP1/SSP2 shape-check bypasses are exploratory;
#   see the caveats in README.md.
################################################################################

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
output_date <- Sys.getenv("TABLE_DATE", "20261124")
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
middle_years <- 2055:2059
end_years <- 2095:2099
start_period <- "2020-2024"
middle_period <- "2055-2059"
end_period <- "2095-2099"
level_periods <- data.table(
  period = c(start_period, middle_period, end_period),
  period_start = c(2020L, 2055L, 2095L),
  period_end = c(2024L, 2059L, 2099L)
)
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
  x <- fread(path, select = c("gcm", "branch", "year", "LE65", "LI65"))
  x <- x[
    gcm == "ENSEMBLE" & branch %chin% names(branch_names) &
      year %in% c(start_years, middle_years, end_years)
  ]
  x[, period := fifelse(
    year %in% start_years, start_period,
    fifelse(year %in% middle_years, middle_period, end_period)
  )]
  result <- x[, .(
    LE65 = mean(LE65),
    LI65plus = mean(LI65),
    years = uniqueN(year)
  ), by = .(branch, period)]
  if (nrow(result) != 6L || any(result$years != 5L) ||
      any(!is.finite(result$LE65)) || any(!is.finite(result$LI65plus))) {
    stop("Pooled ENSEMBLE LE65/LI65+ levels are incomplete in ",
         path, call. = FALSE)
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
      any(!is.finite(x$le_contribution)) ||
      any(!is.finite(x$li_contribution))) {
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
  level_li <- function(branch, period) {
    branch_value <- branch
    period_value <- period
    level_data$LI65plus[
      level_data$branch == branch_value & level_data$period == period_value
    ]
  }
  levels_long <- melt(
    copy(level_data),
    id.vars = c("branch", "period", "years"),
    measure.vars = c("LE65", "LI65plus"),
    variable.name = "measure",
    value.name = "level_years"
  )
  levels_long <- merge(levels_long, level_periods, by = "period", all.x = TRUE)
  levels_long[, `:=`(
    ssp = ssp,
    scenario = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")[ssp],
    exploratory = ssp < 3L,
    level = level,
    id = id,
    geography = label,
    n_years = years,
    units = "years",
    definition = paste(
      "Mean annual ENSEMBLE value over the five years from pooled",
      "04_le_li_levels.csv."
    )
  )]
  levels_long[, years := NULL]
  level_changes_long <- rbindlist(lapply(c("LE65", "LI65plus"), function(measure) {
    start_with <- level_data[
      branch == "with_cc" & period == start_period, get(measure)
    ]
    end_with <- level_data[
      branch == "with_cc" & period == end_period, get(measure)
    ]
    start_without <- level_data[
      branch == "without_cc" & period == start_period, get(measure)
    ]
    end_without <- level_data[
      branch == "without_cc" & period == end_period, get(measure)
    ]
    if (any(lengths(list(start_with, end_with, start_without, end_without)) != 1L)) {
      stop("Incomplete start/end levels for ", label, " ", measure, ".",
           call. = FALSE)
    }
    change_with <- end_with - start_with
    change_without <- end_without - start_without
    data.table(
      ssp = ssp,
      scenario = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")[ssp],
      exploratory = ssp < 3L,
      level = level,
      id = id,
      geography = label,
      measure = measure,
      estimand = c(
        "change_with_cc_2020-24_to_2095-99",
        "change_without_cc_2020-24_to_2095-99",
        "warming_effect_on_change_2020-24_to_2095-99"
      ),
      start_level_years = c(start_with, start_without, NA_real_),
      end_level_years = c(end_with, end_without, NA_real_),
      change_years = c(change_with, change_without, change_with - change_without)
    )
  }))

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
  li_cause_months <- decomposition[, .(
    months = 12 * sum(li_contribution)
  ), by = cause]
  li_cold_months <- li_cause_months[cause %chin% cold_causes, sum(months)]
  li_heat_months <- li_cause_months[cause %chin% heat_causes, sum(months)]
  li_together_months <- li_cold_months + li_heat_months
  li_temperature_months <-
    li_cause_months[cause %chin% temperature_causes, sum(months)]
  li_total_months <- sum(li_cause_months$months)
  li_change_with <- level_li("with_cc", end_period) -
    level_li("with_cc", start_period)
  li_change_without <- level_li("without_cc", end_period) -
    level_li("without_cc", start_period)
  li_change_with_months <- 12 * li_change_with
  li_change_without_months <- 12 * li_change_without
  if (!is.finite(li_change_with_months) ||
      abs(li_change_with_months) < 1e-9 ||
      any(!is.finite(c(li_cold_months, li_heat_months, li_together_months)))) {
    stop("Invalid LI65+ changes or contributions for SSP", ssp, " ",
         label, ".", call. = FALSE)
  }
  check_numeric(
    "LI65+ level change closes to within-branch decomposition",
    label, ssp, li_total_months, li_change_with_months, 1e-3,
    "Sum of all Part 04 li_contribution values over 2020-24 to 2095-99."
  )
  check_numeric(
    "LI65+ cold + heat equals five-range contribution",
    label, ssp, li_together_months, li_temperature_months, 1e-8,
    "Cold is extreme + moderate; heat is moderate + normal-extreme + extreme heat."
  )
  direction_expected <- if (ssp == 1L && id == "Western") 1 else -1
  direction_ok <- if (direction_expected > 0) li_together_months > 0 else
    li_together_months < 0
  add_check(
    "Simon draft direction: cold and heat LI65+ contribution",
    label, ssp, li_together_months, direction_expected, NA_real_,
    if (direction_ok) "PASS" else "FAIL",
    "Cold and heat together lower LI65+ except Western Europe SSP1."
  )

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
      exploratory = ssp < 3L,
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
      exploratory = ssp < 3L,
      geography = label, gain_in_years = g1,
      cold_pct = 100 * cold_months / (12 * g1),
      heat_pct = 100 * heat_months / (12 * g1),
      cold_and_heat_together_pct = 100 * together_months / (12 * g1),
      reduction_due_to_warming_pct = 100 * all_loss / g0,
      note = "SSP1 uses exploratory 0.9995 replacement; SSP2 includes 12 Swiss cities with mortality-shape check bypassed. Neither is a production run."
    ),
    sx = data.table(
      ssp = ssp, scenario = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")[ssp],
      exploratory = ssp < 3L,
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
    )],
    levels_long = levels_long,
    level_changes_long = level_changes_long,
    li_decomposition = data.table(
      ssp = ssp,
      scenario = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")[ssp],
      exploratory = ssp < 3L,
      level = level,
      id = id,
      geography = label,
      LI65plus_start_with_cc_years = level_li("with_cc", start_period),
      LI65plus_end_with_cc_years = level_li("with_cc", end_period),
      LI65plus_change_with_cc_years = li_change_with,
      LI65plus_change_without_cc_years = li_change_without,
      LI65plus_warming_effect_on_change_months =
        12 * (li_change_with - li_change_without),
      cold_contribution_months = li_cold_months,
      heat_contribution_months = li_heat_months,
      cold_and_heat_together_contribution_months = li_together_months,
      cold_contribution_pct_of_LI_change =
        100 * li_cold_months / li_change_with_months,
      heat_contribution_pct_of_LI_change =
        100 * li_heat_months / li_change_with_months,
      cold_and_heat_together_pct_of_LI_change =
        100 * li_together_months / li_change_with_months,
      total_decomposition_change_months = li_total_months,
      level_change_months = li_change_with_months
    )
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
levels_long <- rbindlist(
  lapply(results, `[[`, "levels_long"), use.names = TRUE
)
level_changes_long <- rbindlist(
  lapply(results, `[[`, "level_changes_long"), use.names = TRUE
)
li_decomposition <- rbindlist(
  lapply(results, `[[`, "li_decomposition"), use.names = TRUE
)

setorder(long, ssp, geography, metric)
setorder(levels_long, ssp, geography, period, branch, measure)
levels_wide <- dcast(
  levels_long,
  ssp + scenario + exploratory + level + id + geography + period +
    period_start + period_end + n_years ~ measure + branch,
  value.var = "level_years"
)
levels_wide[, `:=`(
  units = "years",
  definition = "Mean annual ENSEMBLE LE65 or LI65+ over the stated five-year period."
)]
setorder(level_changes_long, ssp, geography, measure, estimand)
level_changes_wide <- dcast(
  level_changes_long,
  ssp + scenario + exploratory + level + id + geography ~ measure + estimand,
  value.var = "change_years"
)
setorder(li_decomposition, ssp, geography)
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
  "LE65/LI65+ levels row count and key uniqueness", "all geographies", NA_integer_,
  observed = nrow(levels_long), expected = 180, tolerance = 0,
  status = if (nrow(levels_long) == 180L &&
               uniqueN(levels_long, by = c(
                 "ssp", "geography", "period", "branch", "measure"
               )) == 180L) "PASS" else "FAIL",
  note = "Two measures, two branches, three periods, three SSPs, and five geographies."
)
add_check(
  "LE65/LI65+ wide levels row count and key uniqueness", "all geographies",
  NA_integer_, observed = nrow(levels_wide), expected = 45, tolerance = 0,
  status = if (nrow(levels_wide) == 45L &&
               uniqueN(levels_wide, by = c("ssp", "geography", "period")) == 45L) {
    "PASS"
  } else {
    "FAIL"
  },
  note = "One row for each scenario, geography, and period."
)
add_check(
  "LE65/LI65+ change row count and key uniqueness", "all geographies", NA_integer_,
  observed = nrow(level_changes_long), expected = 90, tolerance = 0,
  status = if (nrow(level_changes_long) == 90L &&
               uniqueN(level_changes_long, by = c(
                 "ssp", "geography", "measure", "estimand"
               )) == 90L) "PASS" else "FAIL",
  note = "Two measures and three endpoint-change estimands for each scenario/geography."
)
add_check(
  "LI65+ decomposition row count and key uniqueness", "all geographies",
  NA_integer_, observed = nrow(li_decomposition), expected = 15, tolerance = 0,
  status = if (nrow(li_decomposition) == 15L &&
               uniqueN(li_decomposition, by = c("ssp", "geography")) == 15L) {
    "PASS"
  } else {
    "FAIL"
  },
  note = "Cold, heat, combined contribution, percent of change, and warming-alone effect."
)
level_value_columns <- c(
  "LE65_with_cc", "LE65_without_cc",
  "LI65plus_with_cc", "LI65plus_without_cc"
)
level_values_ok <- all(level_value_columns %in% names(levels_wide)) &&
  all(vapply(levels_wide[, ..level_value_columns], function(x) all(is.finite(x)),
             logical(1)))
add_check(
  "LE65/LI65+ levels are finite for all branches and periods",
  "all geographies", NA_integer_, observed = as.numeric(level_values_ok),
  expected = 1, tolerance = 0,
  status = if (level_values_ok) "PASS" else "FAIL",
  note = "Annual means of the ENSEMBLE rows in 04_le_li_levels.csv."
)
for (metric_pair in list(
  c("change_with_cc_2020-24_to_2095-99",
    "gain_with_climate_change_from_annual_levels"),
  c("change_without_cc_2020-24_to_2095-99",
    "gain_without_climate_change_from_annual_levels")
)) {
  level_change <- level_changes_long[
    measure == "LE65" & estimand == metric_pair[[1L]],
    .(ssp, geography, level_change_years = change_years)
  ]
  table2_gain <- long[
    metric == metric_pair[[2L]],
    .(ssp, geography, table2_gain_years = value)
  ]
  compared <- merge(level_change, table2_gain, by = c("ssp", "geography"))
  max_gap <- if (nrow(compared) == 15L) {
    max(abs(compared$level_change_years - compared$table2_gain_years))
  } else {
    Inf
  }
  check_numeric(
    "LE65 level changes agree with Table 2 annual-level gains",
    "all geographies", NA_integer_, max_gap, 0, 1e-12, metric_pair[[2L]]
  )
}
flags_match <- function(x) {
  nrow(x) > 0L &&
    all(x[ssp %in% 1:2, exploratory]) &&
    all(!x[ssp == 3L, exploratory])
}
exploratory_flags_ok <- all(vapply(
  list(
    long, wide, sx, levels_long, levels_wide, level_changes_long,
    level_changes_wide, li_decomposition
  ), flags_match, logical(1)
))
add_check(
  "exploratory flags mark SSP1 and SSP2 only", "all geographies", NA_integer_,
  observed = as.numeric(exploratory_flags_ok), expected = 1, tolerance = 0,
  status = if (exploratory_flags_ok) "PASS" else "FAIL",
  note = "All Table 2 and Table Sx rows carry an explicit logical flag."
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
fwrite(levels_long, file.path(
  data_dir, paste0("table02_le_li_levels_long_", output_date, ".csv")
))
fwrite(levels_wide, file.path(data_dir, "table02_le_li_levels_wide.csv"))
fwrite(level_changes_long, file.path(
  data_dir, paste0("table02_le_li_level_changes_long_", output_date, ".csv")
))
fwrite(level_changes_wide, file.path(
  data_dir, "table02_le_li_level_changes_wide.csv"
))
fwrite(li_decomposition, file.path(
  data_dir, "table02_li65_change_decomposition.csv"
))
fwrite(checks, file.path(data_dir, "table02_checks.csv"))

commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)
if (!is.null(attr(commit, "status")) || length(commit) != 1L || !nzchar(commit)) {
  stop("Could not determine the analysis Git commit.", call. = FALSE)
}
git_status <- system2(
  "git", c("status", "--porcelain", "--untracked-files=normal"),
  stdout = TRUE, stderr = FALSE
)
if (!is.null(attr(git_status, "status"))) {
  stop("Could not determine the analysis Git worktree status.", call. = FALSE)
}
worktree <- if (length(git_status)) "modified" else "clean"
uncertainty_note <- if (exists("uncertainty_comparison")) {
  all <- uncertainty_comparison[metric == "all_erosion_ratio_pct"]
  heat <- uncertainty_comparison[metric == "heat_erosion_ratio_pct"]
  sprintf(
    paste(
      "Table 2 uses pooled 19-GCM ENSEMBLE central schedules: all-temperature",
      "warming reduction %.4f%% (4.3%% rounded). The heat-only central value",
      "in Table Sx is %.4f%%. The 2026-11-10 summary's 19-GCM by 500-draw",
      "experiment medians are %.4f%% for all temperatures (3.9%% rounded)",
      "and %.4f%% for heat (7.0%% rounded). These are different summaries,",
      "not conflicting estimates: the paper tables use the ENSEMBLE central",
      "schedule, not the median over GCM-by-ERF-draw experiments."
    ),
    all$central_pct, heat$central_pct, all$run_median_pct, heat$run_median_pct
  )
} else {
  "No prior uncertainty summary was available for comparison."
}
readme <- c(
  "# Table 2 and Table Sx",
  "",
  paste("Analysis commit:", commit),
  paste("Analysis worktree:", worktree),
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
    "The LE65 and LI65+ level exports give the means of annual ENSEMBLE values",
    "in 2020-24, 2055-59, and 2095-99 from `04_le_li_levels.csv`, with and",
    "without climate change. The endpoint LE65 changes match the annual-level",
    "gain metrics in the Table 2 long file. These annual-level changes need not",
    "equal period-mean-rate gains; Table 2 columns and reductions use the",
    "period-mean-rate G0/G1 to remain consistent with Table 1."
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
    "The LI65+ decomposition export gives cold, heat, and their combined",
    "contributions in months (12 times the within-branch Part 04 contributions),",
    "and as percentages of the with-CC LI65+ endpoint change. The warming-alone",
    "effect is the with-CC LI65+ change minus the without-CC change, in months.",
    "The closure check compares all Part 04 contributions with the annual-level",
    "LI65+ change."
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
    "table output. The Table 2, Table Sx, level, level-change, and LI65+",
    "decomposition CSVs carry an `exploratory` flag (`TRUE` for SSP1 and SSP2;",
    "`FALSE` for SSP3)."
  ),
  paste(
    "Table Sx reports the with-CC ENSEMBLE LE65 gain, cold and heat",
    "contributions in months, their total, and the all-temperature and",
    "heat-only reductions due to warming; all five geographies are populated."
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
