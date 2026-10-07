#!/usr/bin/env Rscript

################################################################################
# Madrid SSP3 coefficient-draw and adaptation pilot
#
# This analysis is deliberately outside the five-step production pipeline. It
# reuses the production calibration, attributable-number and Lloyd life-table
# definitions without changing canonical outputs. The experiment propagates
# every available Masselot ERF coefficient draw for Madrid through one validated
# GCM (GFDL_ESM4), then applies the published heat-only excess-RR attenuation
# rule at 0%, 10%, 50% and 90%.
################################################################################

source("pipeline/00_pkg_params.R")

source_root <- Sys.getenv("SOURCE_ROOT", "/home/SHARED/temperature-mortality-le-li")
canonical_root <- Sys.getenv("CANONICAL_ROOT", file.path(source_root, "results/europe"))
output_dir <- Sys.getenv(
  "OUTPUT_DIR",
  "results/europe/auxiliary/vig_phase2_20261007_madrid_mc_adaptation"
)
n_draws_request <- as.integer(Sys.getenv("N_DRAWS", "0"))
chunk_size <- as.integer(Sys.getenv("DRAW_CHUNK_SIZE", "100"))
if (!is.finite(n_draws_request) || n_draws_request < 0L) {
  stop("N_DRAWS must be zero (all) or a positive integer.", call. = FALSE)
}
if (!is.finite(chunk_size) || chunk_size < 1L) {
  stop("DRAW_CHUNK_SIZE must be a positive integer.", call. = FALSE)
}
if (city_id != "ES001C" || ssp_name != "3" || gcm_name != "GFDL_ESM4") {
  stop("This validated pilot is fixed to ES001C / SSP3 / GFDL_ESM4.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
started_at <- Sys.time()
adaptation_pct <- c(0, 10, 50, 90)
baseline_years <- 2020:2024
late_years <- 2095:2099
analysis_years <- c(baseline_years, late_years)
annualization_rule <- "sum(daily AN) / 365 (29 Feb removed)"

input_paths <- list(
  prep = file.path(source_root, "data/prep_data.RData"),
  temperature = file.path(source_root, "data/tmeanproj.gz.parquet"),
  draws = file.path(source_root, "data/coef_simu.csv"),
  central = "data/coefs.csv",
  city_results = "data/city_results.csv",
  demography_grouped = file.path(canonical_root, "ssp3", city_id, "demography", "00_demography_grouped.csv"),
  demography_single = file.path(canonical_root, "ssp3", city_id, "demography", "00_demography_single_age.csv"),
  production_levels = file.path(canonical_root, "ssp3", city_id, gcm_name, "04_le_li_levels.csv")
)
missing_inputs <- names(input_paths)[!vapply(input_paths, file.exists, logical(1))]
if (length(missing_inputs)) {
  stop("Missing inputs: ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}

# Load the production-prepared observed series and fixed thresholds.
prep_env <- new.env(parent = emptyenv())
load(input_paths$prep, envir = prep_env)
obs_data <- as.data.table(prep_env$obs_data)
thresholds <- as.data.table(prep_env$thresholds)
rm(prep_env)

city_meta <- unique(fread(input_paths$city_results)[
  URAU_CODE == city_id & agegroup %in% agelabs
])
city_meta <- city_meta[match(agelabs, agegroup)]
if (nrow(city_meta) != length(agelabs) || anyNA(city_meta)) {
  stop("Madrid city metadata are incomplete.", call. = FALSE)
}

city_thresholds <- thresholds[URAU_CODE == city_id & agegroup %in% agelabs]
city_thresholds <- city_thresholds[match(agelabs, agegroup)]
if (nrow(city_thresholds) != length(agelabs) || anyNA(city_thresholds)) {
  stop("Madrid threshold rows are incomplete.", call. = FALSE)
}

dem_grouped <- fread(input_paths$demography_grouped)[
  geo_id == city_id & ssp == as.integer(ssp_name) & agegroup %in% agelabs
]
dem_single <- fread(input_paths$demography_single)[
  geo_id == city_id & ssp == as.integer(ssp_name) & age %in% age_levels
]
if (nrow(dem_grouped) != length(future_years) * length(agelabs)) {
  stop("Grouped demography does not have the expected year-by-age-group grid.", call. = FALSE)
}
if (nrow(dem_single) != length(future_years) * length(age_levels)) {
  stop("Single-age demography does not have the expected year-by-age grid.", call. = FALSE)
}

# Coefficient draws. The source has 1,000 draws; N_DRAWS=0 uses all of them.
coef_draws <- fread(input_paths$draws)[
  URAU_CODE == city_id & agegroup %in% agelabs
]
available_draws <- sort(unique(coef_draws$sim))
if (n_draws_request > 0L) {
  available_draws <- head(available_draws, n_draws_request)
  coef_draws <- coef_draws[sim %in% available_draws]
}
n_draws <- length(available_draws)
expected_draw_rows <- n_draws * length(agelabs)
draw_grid_complete <- nrow(coef_draws) == expected_draw_rows &&
  !anyDuplicated(coef_draws, by = c("agegroup", "sim")) &&
  setequal(unique(coef_draws$agegroup), agelabs)
coef_cols <- paste0("b", 1:5)
draw_values_finite <- all(is.finite(as.matrix(coef_draws[, ..coef_cols])))

coef_central <- fread(input_paths$central)[
  URAU_CODE == city_id & agegroup %in% agelabs
]
coef_central <- coef_central[match(agelabs, agegroup)]
if (nrow(coef_central) != length(agelabs) || anyNA(coef_central[, ..coef_cols])) {
  stop("Central Madrid coefficients are incomplete.", call. = FALSE)
}

coef_means <- coef_draws[, lapply(.SD, mean), by = agegroup, .SDcols = coef_cols]
coef_means <- coef_means[match(agelabs, agegroup)]
coef_validation <- rbindlist(lapply(coef_cols, function(cc) {
  data.table(
    agegroup = agelabs,
    coefficient = cc,
    central = coef_central[[cc]],
    draw_mean = coef_means[[cc]],
    difference = coef_means[[cc]] - coef_central[[cc]]
  )
}))
coef_validation[, abs_difference := abs(difference)]
coef_validation[, draw_mean_correlation := cor(central, draw_mean), by = coefficient]
fwrite(coef_validation, file.path(output_dir, "coefficient_validation.csv"))

# Production ISIMIP3BASD function, copied without scientific changes from
# pipeline/01_attribution.R (Masselot & Gasparrini, 2025).
isimip3 <- function(obshist, simhist, simfut, yearobshist, yearsimhist, yearsimfut, detrend = TRUE) {
  if (detrend) {
    obstrend <- lm(obshist ~ yearobshist, na.action = na.exclude) |> predict() |> scale(scale = FALSE)
    simhisttrend <- lm(simhist ~ yearsimhist, na.action = na.exclude) |> predict() |> scale(scale = FALSE)
    simfuttrend <- lm(simfut ~ yearsimfut, na.action = na.exclude) |> predict() |> scale(scale = FALSE)
    obshist <- obshist - obstrend
    simhist <- simhist - simhisttrend
    simfut <- simfut - simfuttrend
  }
  ecdfobs <- ecdf(obshist)(obshist)
  deltaadd <- quantile(simfut, ecdfobs, na.rm = TRUE) - quantile(simhist, ecdfobs, na.rm = TRUE)
  obsfut <- deltaadd + obshist
  simfutcdf <- pnorm(simfut, mean(simfut, na.rm = TRUE), sd(simfut, na.rm = TRUE))
  calsimfut <- qnorm(p = simfutcdf, mean = mean(obsfut, na.rm = TRUE), sd = sd(obsfut, na.rm = TRUE))
  if (detrend) calsimfut <- calsimfut + simfuttrend
  calsimfut
}

obs_city <- obs_data[URAU_CODE == city_id]
obs_city[, `:=`(
  year = as.integer(format(date, "%Y")),
  month = as.integer(format(date, "%m")),
  month_day = format(date, "%m-%d")
)]
obs_hist <- obs_city[year %between% histrange & month_day != "02-29"]

proj_ds <- open_dataset(input_paths$temperature)
tmean_all <- proj_ds |>
  filter(URAU_CODE == city_id, ssp %in% c("hist", ssp_name)) |>
  collect() |>
  as.data.table()
gcm_col <- paste0("tas_", gcm_name)
if (!gcm_col %in% names(tmean_all)) stop("GCM temperature column is missing.", call. = FALSE)
tmean_all <- tmean_all[, .(date, ssp, tmean = get(gcm_col))]
tmean_all[, `:=`(
  year = as.integer(format(date, "%Y")),
  month = as.integer(format(date, "%m")),
  month_day = format(date, "%m-%d")
)]
tmean_all <- tmean_all[month_day != "02-29"]
tmean_future <- tmean_all[ssp == ssp_name & year %in% future_years]
hist_sim <- tmean_all[ssp == "hist" & year %between% histrange]

cal_src <- rbind(
  tmean_all[ssp == "hist" & year %between% histrange],
  tmean_all[ssp == ssp_name & year >= min(projrange)]
)
cal_src[, calperiod := cut(year, c(histrange[1], projrange), right = FALSE)]
if (anyNA(cal_src$calperiod)) stop("Calibration-period assignment failed.", call. = FALSE)
cal_src[, full := {
  m <- .BY$month
  obs_m <- obs_hist[month == m]
  sim_h_m <- hist_sim[month == m]
  isimip3(
    obshist = obs_m$tmean_obs, simhist = sim_h_m$tmean, simfut = tmean,
    yearobshist = obs_m$year, yearsimhist = sim_h_m$year, yearsimfut = year
  )
}, by = .(month, calperiod)]
cal_src[, year5 := (year %/% 5L) * 5L]
cal_src[, demo := {
  m <- .BY$month
  ref <- cal_src[month == m & year5 == counterfactual_ref_year5]
  isimip3(
    obshist = ref$full, simhist = full, simfut = full,
    yearobshist = ref$year, yearsimhist = year, yearsimfut = year
  )
}, by = .(month, year5)]

temp_variants <- list(
  with_cc = merge(
    copy(tmean_future), cal_src[, .(date, tmean_variant = full)],
    by = "date", all.x = TRUE, sort = FALSE
  ),
  without_cc = merge(
    copy(tmean_future), cal_src[, .(date, tmean_variant = demo)],
    by = "date", all.x = TRUE, sort = FALSE
  )
)
if (any(vapply(temp_variants, function(x) anyNA(x$tmean_variant), logical(1)))) {
  stop("A calibrated temperature branch contains missing values.", call. = FALSE)
}
if (any(vapply(temp_variants, function(x) {
  any(x[, .N, by = year]$N != 365L)
}, logical(1)))) {
  stop("Daily temperature coverage is not 365 days per year.", call. = FALSE)
}

tper <- quantile(obs_city$tmean_obs, predper / 100, na.rm = TRUE)
knots <- tper[paste0(varper, ".0%")]
bound <- range(tper)

# Draw 0 is the published central coefficient vector; positive IDs are the
# simulation rows. Keeping it in the same array makes validation direct.
experiment_ids <- c(0L, available_draws)
n_experiments <- length(experiment_ids)
an_total <- array(
  NA_real_,
  dim = c(n_experiments, length(adaptation_pct), length(branch_levels),
          length(analysis_years), length(agelabs)),
  dimnames = list(
    draw = as.character(experiment_ids),
    adaptation = as.character(adaptation_pct),
    branch = branch_levels,
    year = as.character(analysis_years),
    agegroup = agelabs
  )
)

for (agegrp in agelabs) {
  message("Computing ", agegrp, "...")
  age_row <- city_thresholds[agegroup == agegrp]
  mmt <- age_row$mmt
  death_annual <- dem_grouped[agegroup == agegrp, .(year, death)]

  draw_age <- coef_draws[agegroup == agegrp][match(available_draws, sim)]
  coef_matrix <- rbind(
    as.matrix(coef_central[agegroup == agegrp, ..coef_cols]),
    as.matrix(draw_age[, ..coef_cols])
  )
  if (nrow(coef_matrix) != n_experiments) stop("Coefficient matrix size mismatch.", call. = FALSE)

  for (branch in branch_levels) {
    age_work <- merge(copy(temp_variants[[branch]]), death_annual, by = "year", all.x = TRUE, sort = FALSE)
    if (anyNA(age_work$death)) stop("Annual deaths are missing after the temperature join.", call. = FALSE)
    b_fut <- onebasis(age_work$tmean_variant, fun = varfun, degree = vardegree, knots = knots, Bound = bound)
    b_mmt <- onebasis(mmt, fun = varfun, degree = vardegree, knots = knots, Bound = bound)
    b_centered <- scale(b_fut, center = b_mmt, scale = FALSE)
    heat_row <- age_work$tmean_variant >= mmt
    target_row <- age_work$year %in% analysis_years

    for (start in seq.int(1L, n_experiments, by = chunk_size)) {
      stop_idx <- min(start + chunk_size - 1L, n_experiments)
      idx <- start:stop_idx
      rr_unadapted <- exp(b_centered %*% t(coef_matrix[idx, , drop = FALSE]))
      rr_unadapted[rr_unadapted < 1] <- 1

      for (adapt in adaptation_pct) {
        rr <- rr_unadapted
        if (adapt > 0) {
          rr[heat_row, ] <- 1 + (rr[heat_row, , drop = FALSE] - 1) * (1 - adapt / 100)
        }
        af <- 1 - 1 / rr
        an_daily <- af * age_work$death
        annual <- rowsum(an_daily[target_row, , drop = FALSE], age_work$year[target_row], reorder = FALSE) / 365
        annual <- annual[match(analysis_years, as.integer(rownames(annual))), , drop = FALSE]
        if (anyNA(annual)) stop("Annual AN aggregation produced missing values.", call. = FALSE)
        an_total[idx, as.character(adapt), branch, , agegrp] <- t(annual)
      }
    }
  }
}

an_values_finite_nonnegative <- all(is.finite(an_total)) && min(an_total) >= -1e-12

# Lloyd/Aburto functions copied without scientific changes from Part 04.
life_expectancy_from_mx_65plus <- function(mx, x, nx = c(rep(1, 100 - 65), Inf), age = 0) {
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  ex[age + 1]
}

sd_from_mx_fun_65_plus <- function(mx, x, nx = c(rep(1, 100 - 65), Inf), age = 0) {
  x_conditional <- x - 65
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  sqrt(sum(dx * (x_conditional + 0.5 - ex[age + 1])^2))
}

dem_target <- dem_single[year %in% analysis_years]
setorder(dem_target, year, age)
pop_matrix <- matrix(dem_target$pop, nrow = length(analysis_years), byrow = TRUE)
death_matrix <- matrix(dem_target$death, nrow = length(analysis_years), byrow = TRUE)
if (any(!is.finite(pop_matrix)) || any(pop_matrix <= 0) ||
    any(!is.finite(death_matrix)) || any(death_matrix < 0)) {
  stop("Single-age demography is invalid.", call. = FALSE)
}

weight_matrices <- lapply(agelabs, function(grp) {
  ages <- age_slices[[grp]]
  d <- dem_target[age %in% ages]
  setorder(d, year, age)
  w <- matrix(d$death, nrow = length(analysis_years), byrow = TRUE)
  w / rowSums(w)
})
names(weight_matrices) <- agelabs
weight_sums_ok <- all(vapply(weight_matrices, function(w) {
  max(abs(rowSums(w) - 1)) <= 1e-12
}, logical(1)))

base_levels <- rbindlist(lapply(seq_along(analysis_years), function(yi) {
  mx <- death_matrix[yi, ] / pop_matrix[yi, ]
  data.table(
    year = analysis_years[yi],
    LE65_without = life_expectancy_from_mx_65plus(mx, age_levels, nx),
    LI65_without = sd_from_mx_fun_65_plus(mx, age_levels, nx)
  )
}))

estimand_rows <- vector("list", n_experiments * length(adaptation_pct))
central_annual_rows <- list()
row_counter <- 0L
minimum_adjusted_death <- Inf

for (ei in seq_len(n_experiments)) {
  for (adapt in adaptation_pct) {
    an_age <- list()
    for (branch in branch_levels) {
      mat <- matrix(0, nrow = length(analysis_years), ncol = length(age_levels))
      for (grp in agelabs) {
        cols <- match(age_slices[[grp]], age_levels)
        group_an <- an_total[ei, as.character(adapt), branch, , grp]
        mat[, cols] <- weight_matrices[[grp]] * as.numeric(group_an)
      }
      an_age[[branch]] <- mat
    }
    adjusted_with <- death_matrix + an_age$with_cc - an_age$without_cc
    minimum_adjusted_death <- min(minimum_adjusted_death, adjusted_with)
    if (any(adjusted_with <= 0)) stop("A climate-adjusted death schedule is non-positive.", call. = FALSE)

    annual_levels <- rbindlist(lapply(seq_along(analysis_years), function(yi) {
      mx_with <- adjusted_with[yi, ] / pop_matrix[yi, ]
      data.table(
        year = analysis_years[yi],
        LE65_with = life_expectancy_from_mx_65plus(mx_with, age_levels, nx),
        LI65_with = sd_from_mx_fun_65_plus(mx_with, age_levels, nx)
      )
    }))
    annual_levels <- merge(base_levels, annual_levels, by = "year")

    if (experiment_ids[ei] == 0L) {
      central_annual_rows[[as.character(adapt)]] <- copy(annual_levels)[, `:=`(
        draw = 0L,
        adaptation_pct = adapt
      )]
    }

    b <- annual_levels[year %in% baseline_years]
    l <- annual_levels[year %in% late_years]
    row_counter <- row_counter + 1L
    estimand_rows[[row_counter]] <- data.table(
      city = city_id,
      city_name = city_meta$LABEL[1],
      ssp = as.integer(ssp_name),
      gcm = gcm_name,
      draw = experiment_ids[ei],
      draw_type = if (experiment_ids[ei] == 0L) "central" else "simulation",
      adaptation_pct = adapt,
      LE65_with_baseline = mean(b$LE65_with),
      LE65_without_baseline = mean(b$LE65_without),
      LE65_gap_baseline = mean(b$LE65_with - b$LE65_without),
      LE65_with_late = mean(l$LE65_with),
      LE65_without_late = mean(l$LE65_without),
      LE65_gap_late = mean(l$LE65_with - l$LE65_without),
      LE65_gain_change = mean(l$LE65_with - l$LE65_without) - mean(b$LE65_with - b$LE65_without),
      LI65_with_baseline = mean(b$LI65_with),
      LI65_without_baseline = mean(b$LI65_without),
      LI65_gap_baseline = mean(b$LI65_with - b$LI65_without),
      LI65_with_late = mean(l$LI65_with),
      LI65_without_late = mean(l$LI65_without),
      LI65_gap_late = mean(l$LI65_with - l$LI65_without),
      LI65_change = mean(l$LI65_with - l$LI65_without) - mean(b$LI65_with - b$LI65_without)
    )
  }
}

estimands <- rbindlist(estimand_rows)
central_annual <- rbindlist(central_annual_rows)
setorder(estimands, draw, adaptation_pct)
setorder(central_annual, adaptation_pct, year)
fwrite(estimands, file.path(output_dir, "pilot_draw_estimands.csv"))
fwrite(central_annual, file.path(output_dir, "central_annual_levels.csv"))

metrics <- c("LE65_gap_late", "LE65_gain_change", "LI65_gap_late", "LI65_change")
summary_long <- melt(
  estimands,
  id.vars = c("draw", "draw_type", "adaptation_pct"),
  measure.vars = metrics,
  variable.name = "estimand",
  value.name = "value"
)
pilot_summary <- summary_long[draw_type == "simulation", .(
  n_draws = .N,
  mean = mean(value),
  sd = sd(value),
  q025 = quantile(value, 0.025),
  median = median(value),
  q975 = quantile(value, 0.975),
  minimum = min(value),
  maximum = max(value),
  proportion_negative = mean(value < 0),
  proportion_positive = mean(value > 0)
), by = .(adaptation_pct, estimand)]
central_values <- summary_long[draw_type == "central", .(central = value), by = .(adaptation_pct, estimand)]
pilot_summary <- merge(pilot_summary, central_values, by = c("adaptation_pct", "estimand"))
pilot_summary[, `:=`(
  central_minus_draw_mean = central - mean,
  central_inside_empirical_95 = central >= q025 & central <= q975
)]
setorder(pilot_summary, estimand, adaptation_pct)
fwrite(pilot_summary, file.path(output_dir, "pilot_summary.csv"))

# Independent reproduction check against the production GFDL central levels.
production_levels <- fread(input_paths$production_levels)[
  branch %in% branch_levels & year %in% analysis_years
]
production_wide <- dcast(production_levels, year ~ branch, value.var = c("LE65", "LI65"))
central_noadapt <- central_annual[adaptation_pct == 0]
repro <- merge(central_noadapt, production_wide, by = "year")
repro[, `:=`(
  LE65_with_abs_error = abs(LE65_with - LE65_with_cc),
  LE65_without_abs_error = abs(LE65_without - LE65_without_cc),
  LI65_with_abs_error = abs(LI65_with - LI65_with_cc),
  LI65_without_abs_error = abs(LI65_without - LI65_without_cc)
)]
fwrite(repro, file.path(output_dir, "production_reproduction.csv"))
max_reproduction_error <- max(repro[, c(
  "LE65_with_abs_error", "LE65_without_abs_error",
  "LI65_with_abs_error", "LI65_without_abs_error"
)], na.rm = TRUE)

# Monotonicity is evaluated for the headline LE65 gain-change estimand. Greater
# attenuation should make the heat-driven climate penalty no more negative.
mono <- estimands[order(draw, adaptation_pct), .(
  monotonic = all(diff(LE65_gain_change) >= -1e-10)
), by = .(draw, draw_type)]
central_monotonic <- mono[draw_type == "central", all(monotonic)]
simulation_monotonic_share <- mono[draw_type == "simulation", mean(monotonic)]

levels_finite_plausible <- all(is.finite(as.matrix(estimands[, setdiff(names(estimands), c(
  "city", "city_name", "gcm", "draw_type"
)), with = FALSE]))) &&
  estimands[, all(LE65_with_baseline > 0 & LE65_with_baseline < 50 &
                  LE65_with_late > 0 & LE65_with_late < 50 &
                  LI65_with_baseline > 0 & LI65_with_baseline < 30 &
                  LI65_with_late > 0 & LI65_with_late < 30)]

validation <- data.table(
  check_name = c(
    "coefficient_draw_grid_complete",
    "coefficient_values_finite",
    "draw_mean_coefficients_reproduce_central",
    "daily_to_annual_an_finite_nonnegative",
    "single_age_death_weights_sum_to_one",
    "central_noadapt_reproduces_production_levels",
    "life_table_outputs_finite_plausible",
    "central_adaptation_response_monotonic",
    "simulation_adaptation_response_monotonic",
    "central_estimands_inside_draw_intervals"
  ),
  status = c(
    if (draw_grid_complete) "PASS" else "FAIL",
    if (draw_values_finite) "PASS" else "FAIL",
    if (max(coef_validation$abs_difference) <= 0.03 &&
        min(coef_validation$draw_mean_correlation) >= 0.99) "PASS" else "FAIL",
    if (an_values_finite_nonnegative) "PASS" else "FAIL",
    if (weight_sums_ok) "PASS" else "FAIL",
    if (max_reproduction_error <= 1e-9) "PASS" else "FAIL",
    if (levels_finite_plausible && minimum_adjusted_death > 0) "PASS" else "FAIL",
    if (central_monotonic) "PASS" else "FAIL",
    if (simulation_monotonic_share >= 0.99) "PASS" else "WARNING",
    if (all(pilot_summary$central_inside_empirical_95)) "PASS" else "WARNING"
  ),
  value = c(
    sprintf("%d/%d rows; %d draws; duplicate keys=%d", nrow(coef_draws), expected_draw_rows, n_draws,
            anyDuplicated(coef_draws, by = c("agegroup", "sim"))),
    sprintf("nonfinite coefficient values=%d", sum(!is.finite(as.matrix(coef_draws[, ..coef_cols])))),
    sprintf("max |draw mean - central|=%.6f; min correlation=%.6f",
            max(coef_validation$abs_difference), min(coef_validation$draw_mean_correlation)),
    sprintf("minimum annual AN=%.9f", min(an_total)),
    sprintf("maximum |sum(weights)-1|=%.3e", max(vapply(weight_matrices, function(w) max(abs(rowSums(w) - 1)), numeric(1)))),
    sprintf("max absolute LE/LI error=%.3e", max_reproduction_error),
    sprintf("minimum adjusted death=%.9f", minimum_adjusted_death),
    sprintf("central monotonic=%s", central_monotonic),
    sprintf("monotonic simulation draws=%d/%d (%.3f)", sum(mono[draw_type == "simulation"]$monotonic), n_draws, simulation_monotonic_share),
    sprintf("central inside 95%% interval=%d/%d", sum(pilot_summary$central_inside_empirical_95), nrow(pilot_summary))
  ),
  threshold = c(
    "exact agegroup x draw grid",
    "all finite",
    "max absolute coefficient difference <= 0.03 and correlation >= 0.99",
    "all finite and >= 0",
    "within 1e-12",
    "<= 1e-9",
    "finite, LE65 in (0,50), LI65 in (0,30), adjusted deaths > 0",
    "all successive changes >= -1e-10",
    ">= 99% of simulation draws",
    "all four adaptations x four estimands"
  )
)
fwrite(validation, file.path(output_dir, "validation_checks.csv"))

# Meeting-ready figures. The band is ERF coefficient-draw uncertainty for one
# GCM only; it is not combined uncertainty and not GCM spread.
plot_data <- pilot_summary[estimand %in% c("LE65_gain_change", "LI65_change")]
adaptation_labels <- sprintf("%d%%", adaptation_pct)
plot_data[, adaptation_label := factor(
  sprintf("%d%%", adaptation_pct),
  levels = adaptation_labels
)]

le_plot <- ggplot(plot_data[estimand == "LE65_gain_change"], aes(adaptation_label, mean * 12)) +
  geom_hline(yintercept = 0, colour = "#8A8A8A", linewidth = 0.4) +
  geom_errorbar(aes(ymin = q025 * 12, ymax = q975 * 12), width = 0.12, colour = "#005F9E", linewidth = 0.8) +
  geom_point(size = 3, colour = "#005F9E") +
  geom_point(aes(y = central * 12), shape = 21, size = 3.2, fill = "white", colour = "#D7191C", stroke = 0.9) +
  labs(
    title = "Madrid: ERF uncertainty and heat-risk attenuation",
    subtitle = "SSP3-7.0, GFDL-ESM4; change in climate effect from 2020-24 to 2095-99",
    x = "Attenuation of heat excess relative risk",
    y = "Change in LE65 climate effect (months)",
    caption = "Blue: mean and empirical 2.5th-97.5th percentiles across ERF coefficient draws. Red ring: central coefficients. One GCM only."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(output_dir, "01_le65_mc_adaptation.png"), le_plot, width = 9, height = 5.4, dpi = 200, bg = "white")

li_plot <- ggplot(plot_data[estimand == "LI65_change"], aes(adaptation_label, mean * 12)) +
  geom_hline(yintercept = 0, colour = "#8A8A8A", linewidth = 0.4) +
  geom_errorbar(aes(ymin = q025 * 12, ymax = q975 * 12), width = 0.12, colour = "#005F9E", linewidth = 0.8) +
  geom_point(size = 3, colour = "#005F9E") +
  geom_point(aes(y = central * 12), shape = 21, size = 3.2, fill = "white", colour = "#D7191C", stroke = 0.9) +
  labs(
    title = "Madrid: lifespan-inequality sensitivity",
    subtitle = "SSP3-7.0, GFDL-ESM4; change in climate effect from 2020-24 to 2095-99",
    x = "Attenuation of heat excess relative risk",
    y = "Change in LI65 climate effect (SD-months)",
    caption = "Blue: mean and empirical 2.5th-97.5th percentiles across ERF coefficient draws. Red ring: central coefficients. One GCM only."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(output_dir, "02_li65_mc_adaptation.png"), li_plot, width = 9, height = 5.4, dpi = 200, bg = "white")

repo_commit <- system("git rev-parse HEAD", intern = TRUE)
repo_branch <- system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
command_text <- sprintf(
  "CITY_ID=%s SSP=%s GCM=%s SOURCE_ROOT=%s CANONICAL_ROOT=%s OUTPUT_DIR=%s N_DRAWS=%d DRAW_CHUNK_SIZE=%d Rscript analysis/madrid_mc_adaptation_pilot.R",
  city_id, ssp_name, gcm_name, source_root, canonical_root, output_dir, n_draws, chunk_size
)
manifest <- c(
  "server=vangelis",
  sprintf("repository=%s", normalizePath(".")),
  sprintf("branch=%s", repo_branch),
  sprintf("commit=%s", repo_commit),
  "script=analysis/madrid_mc_adaptation_pilot.R",
  sprintf("started_at=%s", format(started_at, tz = "UTC", usetz = TRUE)),
  sprintf("finished_at=%s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("city=%s", city_id),
  sprintf("ssp=%s", ssp_name),
  sprintf("gcm=%s", gcm_name),
  sprintf("draws=%d", n_draws),
  sprintf("adaptation_pct=%s", paste(adaptation_pct, collapse = ",")),
  sprintf("baseline_years=%s", paste(range(baseline_years), collapse = "-")),
  sprintf("late_years=%s", paste(range(late_years), collapse = "-")),
  sprintf("annualization=%s", annualization_rule),
  sprintf("inputs=%s", paste(unlist(input_paths), collapse = ";")),
  sprintf("outputs=%s", paste(sort(list.files(output_dir)), collapse = ";")),
  sprintf("command=%s", command_text),
  "interpretation=Empirical ERF coefficient-draw distribution conditional on SSP3-7.0, Madrid demography, GFDL-ESM4 and fixed heat excess-RR attenuation.",
  "limitations=Single city and single GCM; not combined climate-model uncertainty; fixed attenuation is a sensitivity, not a forecast; PCLM structural uncertainty excluded."
)
writeLines(manifest, file.path(output_dir, "provenance.txt"))

if (any(validation$status == "FAIL")) {
  stop("Pilot failed validation; see validation_checks.csv.", call. = FALSE)
}

message("Pilot complete: ", output_dir)
print(pilot_summary)
print(validation)
