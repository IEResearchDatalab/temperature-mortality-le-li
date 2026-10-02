#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# Pipeline Part 03: Dataset for analysis (Lloyd et al. 2024, Fig S1)
#   Produces the same population-and-deaths-by-cause table at every geography.
#   City mode combines Parts 00 and 02. Collect mode builds Object 1 from all
#   completed cities. Geography mode sums city population and deaths by cause
#   before any mortality rate or life table is calculated (Simon, 25 Sep 2026).
#
#   ANALYSIS_MODE=city       one city/GCM (default)
#   ANALYSIS_MODE=collect    collect completed ensemble city tables as Object 1
#   ANALYSIS_MODE=geography  one country, region or Europe from Object 1
#
################################################################################

source("pipeline/00_pkg_params.R")

analysis_mode <- tolower(Sys.getenv("ANALYSIS_MODE", "city"))
if (!analysis_mode %in% c("city", "collect", "geography")) {
  stop("ANALYSIS_MODE must be city, collect or geography.", call. = FALSE)
}

#------------------------
# COLLECT CITY DATASETS (OBJECT 1)
#------------------------

if (analysis_mode == "collect") {
  root <- Sys.getenv("BATCH_ROOT", "results/europe")
  collected_dir <- file.path(root, "collected")
  dir.create(collected_dir, recursive = TRUE, showWarnings = FALSE)

  meta <- unique(fread("data/city_results.csv")[, .(
    city = URAU_CODE,
    city_name = LABEL,
    country = CNTR_CODE,
    region
  )])
  ensemble_dirs <- Sys.glob(file.path(root, "ssp*", "*", "ENSEMBLE"))
  ensemble_dirs <- ensemble_dirs[file.exists(file.path(ensemble_dirs, ".done"))]
  if (!length(ensemble_dirs)) stop("No completed ensemble runs under ", root, call. = FALSE)

  runs <- data.table(
    dir = ensemble_dirs,
    city = basename(dirname(ensemble_dirs)),
    ssp = as.integer(sub("ssp", "", basename(dirname(dirname(ensemble_dirs)))))
  )
  object1 <- rbindlist(lapply(seq_len(nrow(runs)), function(i) {
    master <- fread(file.path(runs$dir[i], "03_master_table.csv"))
    master[, .(
      city = geo_id,
      ssp = runs$ssp[i],
      scenario = branch,
      year,
      age,
      cause = range,
      deaths = deaths_component,
      pop
    )]
  }))
  object1 <- merge(meta, object1, by = "city")
  setorder(object1, city, ssp, scenario, year, age, cause)
  write_parquet(object1, file.path(collected_dir, "object1_dataset.parquet"))
  message("Saved Object 1 for ", uniqueN(object1$city), " cities to ", collected_dir)
  quit(save = "no")
}

#------------------------
# AGGREGATE CITY COUNTS TO A REPORTING GEOGRAPHY
#------------------------

if (analysis_mode == "geography") {
  if (!geo_level %in% c("country", "region", "europe")) {
    stop("Geography mode requires GEO_LEVEL=country, region or europe.", call. = FALSE)
  }

  object1_file <- Sys.getenv(
    "OBJECT1_FILE",
    file.path("results/europe/collected", "object1_dataset.parquet")
  )
  if (!file.exists(object1_file)) stop("Object 1 not found: ", object1_file, call. = FALSE)

  meta <- unique(fread("data/city_results.csv")[, .(
    city = URAU_CODE,
    city_name = LABEL,
    country = CNTR_CODE,
    region
  )])
  expected_cities <- switch(
    geo_level,
    europe = meta$city,
    region = meta[region == geo_id, city],
    country = meta[country == geo_id, city]
  )
  if (!length(expected_cities)) stop("Unknown or empty geography: ", geo_id, call. = FALSE)

  query <- open_dataset(object1_file) %>%
    filter(ssp == !!as.integer(ssp_name), city %in% !!expected_cities) %>%
    select(city, scenario, year, age, cause, deaths, pop)
  city_data <- as.data.table(collect(query))
  included_cities <- unique(city_data$city)
  if (!setequal(included_cities, expected_cities)) {
    stop("City coverage does not match the requested geography.", call. = FALSE)
  }

  deaths <- city_data[, .(
    deaths_component = sum(deaths)
  ), by = .(scenario, year, age, cause)]
  population <- unique(city_data[cause == "rest", .(
    city, scenario, year, age, pop
  )])[, .(pop = sum(pop)), by = .(scenario, year, age)]
  master <- merge(deaths, population, by = c("scenario", "year", "age"))
  setnames(master, c("scenario", "cause"), c("branch", "range"))

  baseline <- master[branch == "without_cc", .(
    death = sum(deaths_component)
  ), by = .(year, age)]
  rest_mortality <- master[range == "rest", .(
    rest = deaths_component
  ), by = .(branch, year, age)]
  temperature_mortality <- master[range != "rest", .(
    temp_deaths = sum(deaths_component)
  ), by = .(branch, year, age)]
  master <- merge(master, baseline, by = c("year", "age"))
  master <- merge(master, rest_mortality, by = c("branch", "year", "age"))
  master <- merge(master, temperature_mortality, by = c("branch", "year", "age"))
  master[, `:=`(
    geo_id = geo_id,
    label = city_name,
    ssp = as.integer(ssp_name),
    gcm = "ENSEMBLE",
    an = fifelse(range == "rest", 0, deaths_component)
  )]
  setcolorder(master, c(
    "geo_id", "label", "ssp", "gcm", "branch", "year", "age", "range",
    "pop", "death", "an", "temp_deaths", "rest", "deaths_component"
  ))
  setorder(master, branch, year, age, range)

  expected_rows <- length(branch_levels) * length(future_years) *
    length(age_levels) * length(cause_levels)
  checks <- data.table(
    check_name = c(
      "city_coverage", "complete_grid", "finite_values", "positive_population",
      "nonnegative_deaths", "population_identical_between_branches",
      "rest_identical_between_branches"
    ),
    status = c(
      if (setequal(included_cities, expected_cities)) "PASS" else "FAIL",
      if (nrow(master) == expected_rows && !anyDuplicated(master, by = c("branch", "year", "age", "range"))) "PASS" else "FAIL",
      if (all(vapply(
        master[, .(pop, death, an, temp_deaths, rest, deaths_component)],
        function(value) all(is.finite(value)),
        logical(1)
      ))) "PASS" else "FAIL",
      if (all(master$pop > 0)) "PASS" else "FAIL",
      if (all(master$deaths_component >= -1e-9)) "PASS" else "FAIL",
      if (master[, max(pop) - min(pop), by = .(year, age, range)][, max(V1)] <= 1e-6) "PASS" else "FAIL",
      if (master[, max(rest) - min(rest), by = .(year, age, range)][, max(V1)] <= 1e-6) "PASS" else "FAIL"
    ),
    value = c(
      sprintf("%d/%d cities", length(included_cities), length(expected_cities)),
      sprintf("%d/%d rows", nrow(master), expected_rows),
      "required fields finite",
      sprintf("minimum population = %.3f", min(master$pop)),
      sprintf("minimum deaths = %.3f", min(master$deaths_component)),
      sprintf("maximum difference = %.3e", master[, max(pop) - min(pop), by = .(year, age, range)][, max(V1)]),
      sprintf("maximum difference = %.3e", master[, max(rest) - min(rest), by = .(year, age, range)][, max(V1)])
    )
  )
  fwrite(checks, file.path(check_dir, "03_master_checks.csv"))
  if (any(checks$status == "FAIL")) {
    stop("Pooled master table failed its invariant checks.", call. = FALSE)
  }
  fwrite(master, file.path(out_dir, "03_master_table.csv"))
  message("Saved pooled master table for ", city_name, " (", length(included_cities), " cities).")
  quit(save = "no")
}

message(sprintf("\n[03] Assembling %s master analysis table...", city_name))

master_file <- file.path(out_dir, "03_master_table.csv")
checks_file <- file.path(check_dir, "03_master_checks.csv")
failures_file <- file.path(check_dir, "03_master_failures.csv")

#----- Load Part 00/02 outputs and restrict to the city/SSP/GCM domain

grouped_dem <- fread(file.path(dem_dir, "00_demography_grouped.csv"))
single_dem <- fread(file.path(dem_dir, "00_demography_single_age.csv"))
single_an <- fread(file.path(out_dir, "02_single_age_an.csv"))

grouped_dem <- grouped_dem[geo_id == city_id & ssp == as.integer(ssp_name)]
single_dem <- single_dem[geo_id == city_id & ssp == as.integer(ssp_name)]
single_an <- single_an[geo_id == city_id & ssp == as.integer(ssp_name) & gcm == gcm_name]

if (!nrow(grouped_dem)) stop("Grouped demography missing for the city.", call. = FALSE)
if (!nrow(single_dem)) stop("Single-age demography missing for the city.", call. = FALSE)
if (!nrow(single_an)) stop("Single-age AN missing for the city.", call. = FALSE)

grouped_dem <- grouped_dem[agegroup %in% c("65-74", "75-84", "85+")]
single_dem <- single_dem[age %in% age_levels]
single_an <- single_an[range %in% range_levels & branch %in% branch_levels]

if (any(!is.finite(grouped_dem$pop)) || any(!is.finite(grouped_dem$death)) || any(grouped_dem$pop < 0) || any(grouped_dem$death < 0)) {
  stop("Grouped demography contains invalid values.", call. = FALSE)
}
if (any(!is.finite(single_dem$pop)) || any(!is.finite(single_dem$death)) || any(single_dem$pop < 0) || any(single_dem$death < 0)) {
  stop("Single-age demography contains invalid values.", call. = FALSE)
}
if (any(!is.finite(single_an$an)) || any(!is.finite(single_an$weight))) {
  stop("Single-age AN contains invalid values.", call. = FALSE)
}

#----- Build the full (branch x year x age x range) domain

dem_grid <- CJ(year = future_years, age = age_levels)
an_grid <- CJ(year = future_years, age = age_levels, range = range_levels, branch = branch_levels)

dem_keys <- single_dem[, .(year, age)]
an_keys <- single_an[, .(year, branch, age, range)]
dem_keys_unique <- unique(dem_keys)
setcolorder(dem_keys_unique, names(dem_grid))
an_keys_unique <- unique(an_keys)
setcolorder(an_keys_unique, names(an_grid))
if (nrow(fsetdiff(dem_grid, dem_keys_unique)) || nrow(fsetdiff(dem_keys_unique, dem_grid)) || anyDuplicated(dem_keys)) {
  stop("Single-age demographic primary keys are incomplete.", call. = FALSE)
}
if (nrow(fsetdiff(an_grid, an_keys_unique)) || nrow(fsetdiff(an_keys_unique, an_grid)) || anyDuplicated(an_keys)) {
  stop("Single-age AN primary keys are incomplete.", call. = FALSE)
}

#----- Merge attributable numbers (Part 02) with single-age demography (Part 00)
master <- merge(
  an_grid,
  single_an,
  by = c("year", "branch", "age", "range"),
  all.x = TRUE,
  sort = FALSE
)

master <- merge(
  master,
  single_dem[, .(year, age, agegroup = agegroup, pop, death, country_pop, country_death, pop_share, death_share, pop_weight, death_weight)],
  by = c("year", "age"),
  all.x = TRUE,
  sort = FALSE
)

master <- merge(
  master,
  grouped_dem[, .(year, agegroup, grouped_pop = pop, grouped_death = death)],
  by = c("year", "agegroup"),
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(master$source_agegroup) || anyNA(master$agegroup) || any(master$source_agegroup != master$agegroup)) {
  stop("Part 00 and Part 02 age-band mappings disagree.", call. = FALSE)
}

without_cc_temp <- master[branch == "without_cc", .(without_cc_temp_deaths = sum(an)), by = .(year, age)]
master <- merge(master, without_cc_temp, by = c("year", "age"), all.x = TRUE, sort = FALSE)

master[, age_temp_deaths := sum(an), by = .(year, branch, age)]
master[, rest := death - without_cc_temp_deaths]
master[branch == "with_cc", adjusted_death := death + (age_temp_deaths - without_cc_temp_deaths)]
master[branch == "without_cc", adjusted_death := death]

#----- Attach identifying columns
master[, `:=`(
  geo_id = city_id,
  label = city_name,
  ssp = as.integer(ssp_name),
  gcm = gcm_name,
  temp_deaths = age_temp_deaths
)]

master[, deaths_component := an]
rest_rows <- unique(master[, .(
  geo_id, label, ssp, gcm, branch, year, age, agegroup, source_agegroup,
  pop, death, grouped_pop, grouped_death, country_pop, country_death,
  pop_share, death_share, pop_weight, death_weight, temp_deaths, rest, adjusted_death, without_cc_temp_deaths
)])
rest_rows[, `:=`(
  range = "rest",
  group_an = 0,
  weight = 0,
  an = 0,
  deaths_component = rest
)]
master <- rbindlist(list(master, rest_rows), use.names = TRUE, fill = TRUE)

setcolorder(master, c(
  "geo_id", "label", "ssp", "gcm", "branch", "year", "age", "agegroup", "source_agegroup",
  "range", "pop", "death", "grouped_pop", "grouped_death", "country_pop", "country_death",
  "pop_share", "death_share", "pop_weight", "death_weight", "group_an", "weight", "an",
  "temp_deaths", "rest", "deaths_component"
))

master[, agegroup := as.character(agegroup)]
master[, source_agegroup := as.character(source_agegroup)]

#----- Invariant checks (project convention: a failing check stops the run)

full_grid <- CJ(branch = branch_levels, year = future_years, age = age_levels, range = cause_levels)

key_cols <- c("geo_id", "label", "ssp", "gcm", "branch", "year", "age", "range")
duplicate_rows <- master[duplicated(master, by = key_cols) | duplicated(master, by = key_cols, fromLast = TRUE)]
observed_grid <- unique(master[, .(branch, year, age, range)])
missing_grid <- fsetdiff(full_grid, observed_grid)
extra_grid <- fsetdiff(observed_grid, full_grid)
dem_branch_cmp <- master[, .(
  pop = unique(pop),
  death = unique(death),
  grouped_pop = unique(grouped_pop),
  grouped_death = unique(grouped_death)
), by = .(branch, year, age)]

branch_delta <- dem_branch_cmp[, .(
  pop_delta = max(pop) - min(pop),
  death_delta = max(death) - min(death),
  grouped_pop_delta = max(grouped_pop) - min(grouped_pop),
  grouped_death_delta = max(grouped_death) - min(grouped_death)
), by = .(year, age)]

required_cols <- c("pop", "death", "grouped_pop", "grouped_death", "an", "temp_deaths", "rest", "deaths_component")
nonfinite_rows <- master[!apply(master[, ..required_cols], 1L, function(x) all(is.finite(x)))]
conservation <- master[, .(component_sum = sum(deaths_component), death = unique(death)), by = .(branch, year, age)]
conservation[, abs_diff := abs(component_sum - death)]
rest_delta_tbl <- master[, .(rest_delta = max(rest) - min(rest)), by = .(year, age)]

checks <- data.table(
  check_name = c(
    "primary_keys_unique",
    "complete_grid",
    "required_columns_finite",
    "without_cc_components_conserve_death",
    "with_cc_components_conserve_adjusted_death",
    "rest_mortality_nonnegative",
    "demographics_identical_across_branches",
    "rest_identical_across_branches"
  ),
  status = c(
    if (!nrow(duplicate_rows)) "PASS" else "FAIL",
    if (!nrow(missing_grid) && !nrow(extra_grid)) "PASS" else "FAIL",
    if (!nrow(nonfinite_rows)) "PASS" else "FAIL",
    if (max(conservation[branch == "without_cc", abs(abs_diff)], na.rm = TRUE) <= 1e-9) "PASS" else "FAIL",
    if (max(master[, .(component_sum = sum(deaths_component), adjusted_death = unique(adjusted_death)), by = .(branch, year, age)][branch == "with_cc", abs(component_sum - adjusted_death)], na.rm = TRUE) <= 1e-9) "PASS" else "FAIL",
    if (min(master$rest) >= -1e-9) "PASS" else "FAIL",
    if (max(branch_delta$pop_delta) == 0 && max(branch_delta$death_delta) == 0 && max(branch_delta$grouped_pop_delta) == 0 && max(branch_delta$grouped_death_delta) == 0) "PASS" else "FAIL",
    if (max(abs(rest_delta_tbl$rest_delta)) <= 1e-12) "PASS" else "FAIL"
  ),
  value = c(
    nrow(duplicate_rows),
    sprintf("missing=%d; extra=%d", nrow(missing_grid), nrow(extra_grid)),
    nrow(nonfinite_rows),
    sprintf("max_abs_diff=%0.3e", max(conservation[branch == "without_cc", abs_diff], na.rm = TRUE)),
    sprintf("max_abs_diff=%0.3e", max(master[, .(component_sum = sum(deaths_component), adjusted_death = unique(adjusted_death)), by = .(branch, year, age)][branch == "with_cc", abs(component_sum - adjusted_death)], na.rm = TRUE)),
    sprintf("min_rest=%g", min(master$rest)),
    sprintf("max_pop_delta=%g; max_death_delta=%g", max(branch_delta$pop_delta), max(branch_delta$death_delta)),
    sprintf("max_abs_diff=%0.3e", max(abs(rest_delta_tbl$rest_delta)))
  ),
  threshold = c(
    "0 duplicate rows",
    sprintf("exactly %d keys", nrow(full_grid)),
    "0 rows with NA/NaN/Inf",
    "<= 1e-9 for without_cc",
    "<= 1e-9 for with_cc",
    ">= -1e-9",
    "zero difference across branches",
    "<= 1e-12 across branches"
  )
)

failures <- data.table()
if (any(checks$status == "FAIL")) {
  failures <- rbindlist(list(
    if (checks$status[checks$check_name == "primary_keys_unique"] == "FAIL") {
      duplicate_rows[, .(
        geo_id, label, ssp, gcm, branch, year, age, range,
        failing_check = "primary_keys_unique",
        observed_value = "duplicate row",
        expected_bound = "unique keys"
      )]
    } else NULL,
    if (checks$status[checks$check_name == "complete_grid"] == "FAIL") {
      missing_grid[, .(
        branch, year, age, range,
        failing_check = "complete_grid",
        observed_value = "missing row",
        expected_bound = "present"
      )]
    } else NULL,
    if (checks$status[checks$check_name == "required_columns_finite"] == "FAIL") {
      nonfinite_rows[, .(
        geo_id, label, ssp, gcm, branch, year, age, range,
        failing_check = "required_columns_finite",
        observed_value = "NA/NaN/Inf present",
        expected_bound = "all required fields finite"
      )]
    } else NULL,
    if (checks$status[checks$check_name == "without_cc_components_conserve_death"] == "FAIL") {
      conservation[branch == "without_cc" & abs(component_sum - death) > 1e-9, .(
        branch, year, age,
        failing_check = "without_cc_components_conserve_death",
        observed_value = component_sum,
        expected_bound = death
      )]
    } else NULL,
    if (checks$status[checks$check_name == "with_cc_components_conserve_adjusted_death"] == "FAIL") {
      master[branch == "with_cc", .(component_sum = sum(deaths_component), adjusted_death = unique(adjusted_death)), by = .(branch, year, age)][abs(component_sum - adjusted_death) > 1e-9, .(
        branch, year, age,
        failing_check = "with_cc_components_conserve_adjusted_death",
        observed_value = component_sum,
        expected_bound = adjusted_death
      )]
    } else NULL,
    if (checks$status[checks$check_name == "rest_mortality_nonnegative"] == "FAIL") {
      master[rest < -1e-9, .(
        geo_id, label, ssp, gcm, branch, year, age, range,
        failing_check = "rest_mortality_nonnegative",
        observed_value = rest,
        expected_bound = ">= -1e-9"
      )]
    } else NULL,
    if (checks$status[checks$check_name == "demographics_identical_across_branches"] == "FAIL") {
      branch_delta[pop_delta != 0 | death_delta != 0 | grouped_pop_delta != 0 | grouped_death_delta != 0, .(
        year, age,
        failing_check = "demographics_identical_across_branches",
        observed_value = sprintf("pop_delta=%g; death_delta=%g; grouped_pop_delta=%g; grouped_death_delta=%g", pop_delta, death_delta, grouped_pop_delta, grouped_death_delta),
        expected_bound = "zero difference"
      )]
    } else NULL,
    if (checks$status[checks$check_name == "rest_identical_across_branches"] == "FAIL") {
      rest_delta_tbl[abs(rest_delta) > 1e-12, .(
        year, age,
        failing_check = "rest_identical_across_branches",
        observed_value = rest_delta,
        expected_bound = "zero difference"
      )]
    } else NULL
  ), fill = TRUE)
}

#----- Persist checks before exposing primary outputs

fwrite(checks, checks_file)
if (nrow(failures)) {
  fwrite(failures, failures_file)
  stop(sprintf("03_master_table.R failed %d invariant(s); see %s", nrow(failures), failures_file), call. = FALSE)
} else {
  if (file.exists(failures_file)) file.remove(failures_file)
  invisible(file.create(failures_file))
}

fwrite(master, master_file)

message("Saved master table to ", master_file)
message("Saved checks to ", checks_file)
