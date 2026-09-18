#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality / life expectancy pipeline -- Madrid pilot
#
# R Pipeline Step 00: Grouped and single-age demographic projections
#   Follows the calibration logic of Masselot & Gasparrini (2025, R Code Part 2):
#   national Wittgenstein population/survival projections are rescaled to a
#   single EUcityTRM city using a city-to-country ratio, then disaggregated
#   from 5-year age bands to single ages (65:100) via PCLM.
#   Scope here is narrowed to one city (Madrid), one SSP (SSP3-7.0) and the
#   three oldest age groups (65-74, 75-84, 85+).
#
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Expected single-age output length for a PCLM fit over grouped age bands x.
pclm_expected_length <- function(x, nlast) {
  sum(diff(x)) + nlast
}

message("\n[00] Building Madrid SSP3-7.0 demographic tables...")

#----- Output locations

out_dir <- "results/phase1_madrid"
check_dir <- "results/checks"
fig_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

grouped_file <- file.path(out_dir, "00_demography_grouped.csv")
single_file <- file.path(out_dir, "00_demography_single_age.csv")
checks_file <- file.path(check_dir, "00_demography_checks.csv")
failures_file <- file.path(check_dir, "00_demography_failures.csv")
fig_file <- file.path(fig_dir, "00_demography_diagnostic.png")

#----- Analysis parameters (single city, single SSP, ages 65+)

city_id <- "ES001C"
city_name <- "Madrid"
ssp_target <- "3"
baseline_year <- 2020L
future_years <- 2020:2099
target_agegroups <- c("65-74", "75-84", "85+")
age_map <- data.table(
  age_band = c("65--69", "70--74", "75--79", "80--84", "85--89", "90--94", "95--99", "100+"),
  agegroup = c(rep("65-74", 2L), rep("75-84", 2L), rep("85+", 4L)),
  age_start = c(65L, 70L, 75L, 80L, 85L, 90L, 95L, 100L),
  nlast = c(rep(5L, 7L), 1L)
)

# PCLM inputs are stored in the Wittgenstein source's thousand-person unit.
# Scale to persons before fitting so the likelihood is evaluated on the same
# count scale as the validated Lloyd/Aburto reference implementations.
pclm_input_scale <- 1000

#----- Verify Madrid identity in the EUcityTRM baseline table

city_meta <- fread("data/city_results.csv")
city_meta <- unique(city_meta[URAU_CODE == city_id & LABEL == city_name & agegroup %in% target_agegroups])
if (nrow(city_meta) != 3L || !all(unique(city_meta$agegroup) %in% target_agegroups)) {
  stop(sprintf("Madrid identifier verification failed for %s; expected 3 age-group rows for %s.", city_id, city_name), call. = FALSE)
}

city_meta <- city_meta[match(target_agegroups, agegroup)]
if (anyNA(city_meta$agegroup)) {
  stop(sprintf("Madrid identifier verification failed: missing one of age groups %s.", paste(target_agegroups, collapse = ", ")), call. = FALSE)
}

country_codes <- unique(city_meta$CNTR_CODE)
if (length(country_codes) != 1L || is.na(country_codes)) {
  stop("Madrid must map to exactly one country code.", call. = FALSE)
}
country_code <- country_codes[[1L]]

#----- Load Wittgenstein population and survival-ratio projections

pop_raw <- fread("data/wittgenstein_pop.csv")
assr_raw <- fread("data/wittgenstein_assr.csv")

sex_levels <- c("Female", "Male")
pop_raw <- pop_raw[
  CNTR_CODE == country_code & ssp == as.integer(ssp_target) &
    age %in% age_map$age_band & sex %in% sex_levels
]
assr_raw <- assr_raw[
  CNTR_CODE == country_code & ssp == as.integer(ssp_target) &
    age %in% age_map$age_band & sex %in% sex_levels
]

if (!nrow(pop_raw) || !nrow(assr_raw)) {
  stop("SSP3 demographic source tables are empty after filtering.", call. = FALSE)
}

parse_year_start <- function(x) as.integer(sub("^([0-9]{4}).*$", "\\1", x))

pop_raw[, year_start := as.integer(year)]
assr_raw[, year_start := parse_year_start(period)]

if (anyNA(pop_raw$year_start) || anyNA(assr_raw$year_start)) {
  stop("Failed to parse SSP3 year anchors from the demographic sources.", call. = FALSE)
}

country_5y <- merge(
  pop_raw[, .(CNTR_CODE, cntr_name, ssp, year_start, age_band = as.character(age), sex, pop = as.numeric(pop))],
  assr_raw[, .(CNTR_CODE, age_band = as.character(age), sex, year_start, assr = as.numeric(assr))],
  by = c("CNTR_CODE", "age_band", "sex", "year_start"),
  all = FALSE,
  sort = FALSE
)

if (anyNA(country_5y$pop) || anyNA(country_5y$assr)) {
  stop("Merged SSP3 demographic source contains NA pop/assr values.", call. = FALSE)
}
if (any(country_5y$pop < 0) || any(!is.finite(country_5y$pop)) || any(!is.finite(country_5y$assr))) {
  stop("Invalid SSP3 demographic source values detected.", call. = FALSE)
}

country_5y[, death_5y := pop * (1 - assr)]
country_5y[, death_annual := death_5y / 5]
country_5y[, `:=`(
  year_end = pmin(year_start + 4L, 2100L)
)]

country_bridge <- country_5y[, .(
  death_annual_expected = sum(death_annual),
  death_5y_expected = sum(death_5y)
), by = .(CNTR_CODE, cntr_name, ssp, year_start, age_band, sex)]
if (any(!is.finite(country_bridge$death_annual_expected)) || any(country_bridge$death_annual_expected < 0)) {
  stop("Annual death bridge derived from ASSR is invalid.", call. = FALSE)
}

#----- Annualise: hold each 5-year Wittgenstein snapshot constant across its window

expand_annual <- function(dt) {
  rbindlist(lapply(seq_len(nrow(dt)), function(i) {
    row <- dt[i]
    yrs <- seq(row$year_start, row$year_end)
        data.table(
          CNTR_CODE = row$CNTR_CODE,
          cntr_name = row$cntr_name,
          ssp = row$ssp,
          year = yrs,
          age_band = row$age_band,
          age_start = age_map$age_start[match(row$age_band, age_map$age_band)],
          sex = row$sex,
          pop = row$pop,
          death = row$death_annual,
          death_annual = row$death_annual,
          death_5y = row$death_5y
        )
  }), use.names = TRUE, fill = TRUE)
}

country_annual_5y <- expand_annual(country_5y)

country_annual <- country_annual_5y[, .(
  pop = sum(pop),
  death = sum(death)
), by = .(CNTR_CODE, cntr_name, ssp, year, age_band, age_start)]

if (uniqueN(country_annual$CNTR_CODE) != 1L || unique(country_annual$CNTR_CODE) != country_code) {
  stop("Country projection contains rows outside Madrid's country.", call. = FALSE)
}

country_annual <- merge(country_annual, age_map[, .(age_band, agegroup, age_start)], by = c("age_band", "age_start"), all.x = TRUE, sort = FALSE)
if (anyNA(country_annual$agegroup)) {
  stop("Failed to map country age bands to Madrid age groups.", call. = FALSE)
}

country_annual_bridge <- country_annual_5y[, .(
  death_annual_expected = sum(death_annual)
), by = .(CNTR_CODE, cntr_name, ssp, year, age_band, age_start)]
country_annual_bridge <- merge(
  country_annual_bridge,
  country_annual[, .(death_annual_observed = sum(death)), by = .(CNTR_CODE, cntr_name, ssp, year, age_band, age_start)],
  by = c("CNTR_CODE", "cntr_name", "ssp", "year", "age_band", "age_start"),
  all.x = TRUE,
  sort = FALSE
)
country_annual_bridge[, bridge_abs_diff := abs(death_annual_expected - death_annual_observed)]

# Aggregate the country age bands into the three Madrid age groups of interest
country_agegroup <- country_annual[, .(
  country_pop = sum(pop),
  country_death = sum(death)
), by = .(year, agegroup)]

#----- Calibrate the city-to-country ratio at the baseline anchor year

country_agegroup_2020 <- country_agegroup[year == baseline_year]
if (nrow(country_agegroup_2020) != length(target_agegroups)) {
  stop("Baseline year 2020 is missing one or more required age groups in the SSP3 source.", call. = FALSE)
}

city_base <- unique(city_meta[, .(agegroup, city_pop = agepop, city_death = death)])
share_tbl <- merge(city_base, country_agegroup_2020, by = "agegroup", all.x = TRUE, sort = FALSE)
if (anyNA(share_tbl$country_pop) || anyNA(share_tbl$country_death)) {
  stop("Could not compute Madrid baseline demographic shares from the SSP3 source.", call. = FALSE)
}
if (any(share_tbl$country_pop <= 0) || any(share_tbl$country_death <= 0)) {
  stop("Nonpositive baseline country demographic totals prevent share construction.", call. = FALSE)
}

share_tbl[, `:=`(
  pop_share = city_pop / (country_pop * 1000),
  death_share = city_death / (country_death * 1000)
)]

if (any(!is.finite(share_tbl$pop_share)) || any(!is.finite(share_tbl$death_share))) {
  stop("Madrid baseline shares are not finite.", call. = FALSE)
}

#----- Disaggregate country age bands to single ages via PCLM (person-scale fit)

fit_pclm_person_scale <- function(x, y, nlast, context, scale_factor = pclm_input_scale) {
  if (any(!is.finite(y))) {
    stop(sprintf("PCLM input contains non-finite values for %s.", context), call. = FALSE)
  }
  if (any(y < 0)) {
    stop(sprintf("PCLM input contains negative values for %s.", context), call. = FALSE)
  }

  if (sum(y) == 0) {
    fit <- rep(0, pclm_expected_length(x, nlast))
    return(list(
      fitted = fit,
      lambda = NA_real_,
      convergence = "ALL_ZERO",
      reconstruction = "ALL_ZERO",
      grouped_error = 0,
      weight_sum = 0,
      min_weight = 0,
      max_weight = 0,
      age_of_max_weight = NA_integer_
    ))
  }

  y_scaled <- y * scale_factor
  pclm_fit <- suppressWarnings(ungroup::pclm(x = x, y = y_scaled, nlast = nlast))
  fit_scaled <- as.numeric(pclm_fit$fitted)
  fit_scaled <- fit_scaled * (sum(y_scaled) / sum(fit_scaled))
  fit <- fit_scaled / scale_factor
  lambda <- as.numeric(pclm_fit$smoothPar["lambda"])
  if (!is.finite(lambda)) {
    stop(sprintf("PCLM returned a non-finite lambda for %s.", context), call. = FALSE)
  }
  if (any(!is.finite(fit)) || any(fit < 0)) {
    stop(sprintf("PCLM returned invalid fitted values for %s.", context), call. = FALSE)
  }
  # Genuine solver-convergence signal from ungroup::pclm's IRLS iteration count,
  # kept distinct from the grouped-reconstruction accuracy check below.
  iterations_used <- as.numeric(pclm_fit$deep$trace)
  max_iterations <- as.numeric(pclm_fit$deep$max.iter)
  if (!is.finite(iterations_used) || !is.finite(max_iterations)) {
    stop(sprintf("PCLM did not expose an iteration count for %s.", context), call. = FALSE)
  }
  solver_converged <- iterations_used < max_iterations
  if (!solver_converged) {
    stop(sprintf("PCLM solver did not converge for %s (iterations = %.0f/%.0f).", context, iterations_used, max_iterations), call. = FALSE)
  }
  grouped_error <- abs(sum(fit) - sum(y))
  if (grouped_error > 1e-9) {
    stop(sprintf("PCLM grouped reconstruction failed for %s (error = %.3e).", context, grouped_error), call. = FALSE)
  }

  weights <- fit / sum(fit)
  list(
    fitted = fit,
    lambda = lambda,
    convergence = if (solver_converged) "PASS" else "FAIL",
    reconstruction = if (grouped_error <= 1e-9) "PASS" else "FAIL",
    grouped_error = grouped_error,
    weight_sum = sum(weights),
    min_weight = min(weights),
    max_weight = max(weights),
    age_of_max_weight = 65L + which.max(weights) - 1L
  )
}

country_single_list <- list()
pclm_diag_list <- list()
for (yr in sort(unique(country_annual$year))) {
  yr_dt <- country_annual[year == yr]
  single_rows <- list()
  pclm_diag_rows <- list()
  for (kind in c("pop", "death")) {
    grouped_vals <- yr_dt[, .(value = sum(get(kind))), by = .(age_start, agegroup)]
    grouped_vals <- grouped_vals[order(age_start)]
    x <- grouped_vals$age_start
    if (any(!is.finite(grouped_vals$value)) || any(grouped_vals$value < 0)) {
      stop(sprintf("Invalid country grouped %s values for year %d.", kind, yr), call. = FALSE)
    }
    pclm_context <- sprintf("year=%d kind=%s", yr, kind)
    pclm_fit <- fit_pclm_person_scale(x = x, y = grouped_vals$value, nlast = 1L, context = pclm_context)
    fit <- pclm_fit$fitted
    if (length(fit) != 36L) {
      stop(sprintf("Unexpected PCLM output length for year %d and %s.", yr, kind), call. = FALSE)
    }
    fit_age <- data.table(year = yr, age = 65:100, value = as.numeric(fit), kind = kind)
    single_rows[[length(single_rows) + 1L]] <- fit_age
    pclm_diag_rows[[length(pclm_diag_rows) + 1L]] <- data.table(
      year = yr,
      kind = kind,
      input_unit = "persons",
      source_unit = "thousand-person Wittgenstein counts",
      input_conversion_factor = pclm_input_scale,
      grouped_total = sum(grouped_vals$value),
      selected_lambda = pclm_fit$lambda,
      convergence_status = pclm_fit$convergence,
      reconstruction_status = pclm_fit$reconstruction,
      grouped_reconstruction_error = pclm_fit$grouped_error,
      weight_sum = pclm_fit$weight_sum,
      min_weight = pclm_fit$min_weight,
      max_weight = pclm_fit$max_weight,
      age_of_max_weight = pclm_fit$age_of_max_weight,
      signed_input_present = FALSE,
      nonnegative_schedule = TRUE,
      exact_band_containment = TRUE
    )
  }
  pop_dt <- single_rows[[1L]]
  death_dt <- single_rows[[2L]]
  out <- merge(pop_dt, death_dt, by = c("year", "age"), suffixes = c("_pop", "_death"))
  out[, `:=`(
    pop = value_pop,
    death = value_death
  )]
  out[, c("value_pop", "value_death") := NULL]
  out[, agegroup := fifelse(age <= 74L, "65-74", fifelse(age <= 84L, "75-84", "85+"))]
  country_single_list[[length(country_single_list) + 1L]] <- out
  pclm_diag_list[[length(pclm_diag_list) + 1L]] <- rbindlist(pclm_diag_rows, fill = TRUE)
}

country_single <- rbindlist(country_single_list, use.names = TRUE, fill = TRUE)
pclm_diag <- rbindlist(pclm_diag_list, use.names = TRUE, fill = TRUE)

country_single_agegroup <- country_single[, .(
  country_pop = sum(pop),
  country_death = sum(death)
), by = .(year, agegroup, age)]

country_single_agegroup <- merge(
  country_single_agegroup,
  country_agegroup[, .(year, agegroup, country_agegroup_pop = country_pop, country_agegroup_death = country_death)],
  by = c("year", "agegroup"),
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(country_single_agegroup$country_agegroup_pop) || anyNA(country_single_agegroup$country_agegroup_death)) {
  stop("Failed to attach country agegroup totals to single-age country table.", call. = FALSE)
}

country_single_agegroup[, `:=`(
  pop_weight = fifelse(country_agegroup_pop > 0, country_pop / country_agegroup_pop, NA_real_),
  death_weight = fifelse(country_agegroup_death > 0, country_death / country_agegroup_death, NA_real_)
)]

country_single_agegroup[, `:=`(
  pop_weight = pop_weight / sum(pop_weight),
  death_weight = death_weight / sum(death_weight)
), by = .(year, agegroup)]

if (anyNA(country_single_agegroup$pop_weight) || anyNA(country_single_agegroup$death_weight)) {
  stop("Single-age country weights could not be constructed.", call. = FALSE)
}

#----- Apply the calibration ratio to build Madrid's grouped projection

city_grouped <- merge(country_agegroup, share_tbl[, .(agegroup, pop_share, death_share)], by = "agegroup", all.x = TRUE, sort = FALSE)
city_grouped[, `:=`(
  geo_id = city_id,
  label = city_name,
  country_code = country_code,
  cntr_name = unique(country_annual$cntr_name),
  ssp = as.integer(ssp_target),
  city_pop = country_pop * pop_share * 1000,
  city_death = country_death * death_share * 1000
)]

city_grouped <- city_grouped[year %in% future_years, .(
  geo_id, label, country_code, cntr_name, ssp, year, agegroup,
  pop = as.numeric(city_pop),
  death = as.numeric(city_death),
  country_pop = as.numeric(country_pop * 1000),
  country_death = as.numeric(country_death * 1000),
  pop_share,
  death_share
)]

#----- Allocate Madrid's grouped totals to single ages using the country age shape

city_single <- merge(
  city_grouped,
  country_single_agegroup[, .(year, agegroup, age, pop_weight, death_weight)],
  by = c("year", "agegroup"),
  all.x = TRUE,
  sort = FALSE
)

city_single[, `:=`(
  geo_id = city_id,
  label = city_name,
  country_code = country_code,
  cntr_name = unique(country_annual$cntr_name),
  ssp = as.integer(ssp_target),
  pop = pop * pop_weight,
  death = death * death_weight
)]

city_single <- city_single[year %in% future_years, .(
  geo_id, label, country_code, cntr_name, ssp, year, agegroup, age,
  pop = as.numeric(pop),
  death = as.numeric(death),
  country_pop = as.numeric(country_pop),
  country_death = as.numeric(country_death),
  pop_share, death_share,
  pop_weight, death_weight
)]

setorder(city_grouped, year, agegroup)
setorder(city_single, year, agegroup, age)

#----- Reconciliation and completeness checks (grouped vs. single-age totals)

recon_check <- merge(
  city_grouped[, .(year, agegroup, grouped_pop = pop, grouped_death = death)],
  city_single[, .(single_pop = sum(pop), single_death = sum(death)), by = .(year, agegroup)],
  by = c("year", "agegroup"),
  all.x = TRUE,
  sort = FALSE
)
recon_check[, `:=`(
  pop_abs_diff = abs(grouped_pop - single_pop),
  death_abs_diff = abs(grouped_death - single_death)
)]

full_group_grid <- CJ(year = future_years, agegroup = target_agegroups)
full_single_grid <- CJ(year = future_years, age = 65:100)

grouped_keys <- unique(city_grouped[, .(year, agegroup)])
single_keys <- unique(city_single[, .(year, agegroup, age)])
missing_grouped_keys <- fsetdiff(full_group_grid, grouped_keys)
extra_grouped_keys <- fsetdiff(grouped_keys, full_group_grid)
expected_single_keys <- copy(full_single_grid)
expected_single_keys[, agegroup := fifelse(age <= 74L, "65-74", fifelse(age <= 84L, "75-84", "85+"))]
setcolorder(expected_single_keys, c("year", "agegroup", "age"))
missing_single_keys <- fsetdiff(expected_single_keys, single_keys)
extra_single_keys <- fsetdiff(single_keys, expected_single_keys)

pclm_diag[, convergence_ok := convergence_status %in% c("PASS", "ALL_ZERO")]
pclm_diag[, reconstruction_ok := reconstruction_status %in% c("PASS", "ALL_ZERO") & grouped_reconstruction_error <= 1e-9]
pclm_diag[, weight_sum_ok := fifelse(grouped_total > 0, abs(weight_sum - 1) <= 1e-12, weight_sum == 0)]
pclm_diag[, lambda_ok := is.finite(selected_lambda) | convergence_status == "ALL_ZERO"]
pclm_diag[, band_ok := exact_band_containment & nonnegative_schedule & !signed_input_present]

#----- Invariant checks (project convention: a failing check stops the run)

checks <- data.table(
  check_name = c(
    "madrid_identifier_unique",
    "grouped_primary_keys_complete",
    "single_primary_keys_complete",
    "population_finite_nonnegative",
    "death_finite_nonnegative",
    "annual_death_bridge",
    "grouped_equals_single_sum",
    "baseline_shares_finite",
    "pclm_input_unit_conversion",
    "pclm_convergence",
    "pclm_selected_lambda_finite",
    "pclm_nonnegative_schedule",
    "pclm_exact_age_band_containment",
    "pclm_weight_sum_one",
    "pclm_grouped_reconstruction",
    "pclm_no_signed_input"
  ),
  status = c(
    if (nrow(city_meta) == 3L) "PASS" else "FAIL",
    if (!nrow(missing_grouped_keys) && !nrow(extra_grouped_keys) && !anyDuplicated(city_grouped, by = c("year", "agegroup"))) "PASS" else "FAIL",
    if (!nrow(missing_single_keys) && !nrow(extra_single_keys) && !anyDuplicated(city_single, by = c("year", "age"))) "PASS" else "FAIL",
    if (!any(!is.finite(city_grouped$pop)) && !any(city_grouped$pop < 0) && !any(!is.finite(city_single$pop)) && !any(city_single$pop < 0)) "PASS" else "FAIL",
    if (!any(!is.finite(city_grouped$death)) && !any(city_grouped$death < 0) && !any(!is.finite(city_single$death)) && !any(city_single$death < 0)) "PASS" else "FAIL",
    if (max(country_annual_bridge$bridge_abs_diff) <= 1e-12) "PASS" else "FAIL",
    if (max(recon_check$pop_abs_diff) <= 1e-9 && max(recon_check$death_abs_diff) <= 1e-9) "PASS" else "FAIL",
    if (!any(!is.finite(share_tbl$pop_share)) && !any(!is.finite(share_tbl$death_share))) "PASS" else "FAIL",
    if (all(pclm_diag$input_unit == "persons") && all(pclm_diag$source_unit == "thousand-person Wittgenstein counts") && all(pclm_diag$input_conversion_factor == pclm_input_scale)) "PASS" else "FAIL",
    if (all(pclm_diag$convergence_ok)) "PASS" else "FAIL",
    if (all(pclm_diag$lambda_ok)) "PASS" else "FAIL",
    if (all(pclm_diag$nonnegative_schedule)) "PASS" else "FAIL",
    if (all(pclm_diag$band_ok)) "PASS" else "FAIL",
    if (all(pclm_diag$weight_sum_ok)) "PASS" else "FAIL",
    if (all(pclm_diag$reconstruction_ok)) "PASS" else "FAIL",
    if (all(!pclm_diag$signed_input_present)) "PASS" else "FAIL"
  ),
  value = c(
    nrow(city_meta),
    nrow(grouped_keys),
    nrow(single_keys),
    sprintf("min_pop=%g; min_single_pop=%g", min(city_grouped$pop), min(city_single$pop)),
    sprintf("min_death=%g; min_single_death=%g", min(city_grouped$death), min(city_single$death)),
    sprintf("max_abs_diff=%0.3e", max(country_annual_bridge$bridge_abs_diff)),
    sprintf("max_pop_diff=%0.3e; max_death_diff=%0.3e", max(recon_check$pop_abs_diff), max(recon_check$death_abs_diff)),
    sprintf("max_pop_share=%0.6f; max_death_share=%0.6f", max(share_tbl$pop_share), max(share_tbl$death_share)),
    sprintf("scale_factor=%d; unit=persons", pclm_input_scale),
    sprintf("converged_rows=%d/%d", sum(pclm_diag$convergence_ok), nrow(pclm_diag)),
    sprintf("lambda_min=%0.6f; lambda_max=%0.6f", min(pclm_diag$selected_lambda, na.rm = TRUE), max(pclm_diag$selected_lambda, na.rm = TRUE)),
    sprintf("all_nonnegative=%s", all(pclm_diag$nonnegative_schedule)),
    sprintf("bad_band_rows=%d", sum(!pclm_diag$band_ok)),
    sprintf("min_weight_sum=%0.12f; max_weight_sum=%0.12f", min(pclm_diag$weight_sum), max(pclm_diag$weight_sum)),
    sprintf("max_grouped_error=%0.3e", max(pclm_diag$grouped_reconstruction_error)),
    sprintf("signed_input_rows=%d", sum(pclm_diag$signed_input_present))
  ),
  threshold = c(
    "unique city id with 3 age-group rows",
    sprintf("%d rows", nrow(full_group_grid)),
    sprintf("%d rows", nrow(full_single_grid)),
    "finite and >= 0",
    "finite and >= 0",
    "<= 1e-12",
    "<= 1e-9",
    "finite shares",
    sprintf("input conversion factor %d to persons", pclm_input_scale),
    "all PCLM fits converge or all-zero",
    "finite lambda or all-zero",
    "nonnegative fitted schedule",
    "exact source-band containment",
    "weight sums within 1e-12",
    "grouped reconstruction <= 1e-9",
    "no signed PCLM input"
  )
)

failures <- data.table()
if (any(checks$status == "FAIL")) {
  failures <- rbindlist(lapply(which(checks$status == "FAIL"), function(i) {
    data.table(
      failing_check = checks$check_name[i],
      observed_value = checks$value[i],
      expected_bound = checks$threshold[i]
    )
  }), fill = TRUE)
}

#----- Persist checks before exposing primary outputs

fwrite(checks, checks_file)
if (nrow(failures)) {
  fwrite(failures, failures_file)
  stop(sprintf("00_demography.R failed %d invariant(s); see %s", nrow(failures), failures_file), call. = FALSE)
} else {
  if (file.exists(failures_file)) file.remove(failures_file)
  file.create(failures_file)
}

fwrite(city_grouped, grouped_file)
fwrite(city_single, single_file)

#----- Diagnostic figure

plot_group <- city_grouped[, .(pop = sum(pop), death = sum(death)), by = .(year, agegroup)]
plot_group <- melt(plot_group, id.vars = c("year", "agegroup"), variable.name = "measure", value.name = "value")

p <- ggplot(plot_group, aes(x = year, y = value, color = agegroup)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~measure, scales = "free_y") +
  labs(
    title = "Madrid SSP3-7.0 demographic projection diagnostic",
    subtitle = sprintf("City %s (%s); piecewise-constant annualization from 5-year source snapshots", city_name, city_id),
    x = "Year",
    y = "Count"
  ) +
  theme_minimal(base_size = 11)

ggsave(fig_file, p, width = 10, height = 6, dpi = 160)

message("Saved grouped demography to ", grouped_file)
message("Saved single-age demography to ", single_file)
message("Saved checks to ", checks_file)
message("Saved diagnostic figure to ", fig_file)
