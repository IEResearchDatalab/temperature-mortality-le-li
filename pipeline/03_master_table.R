#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# Pipeline Part 03: Dataset for analysis (Lloyd et al. 2024, Fig S1)
#   Population and deaths by cause (4 temperature ranges + rest) by single age,
#   year and branch. Rest = deaths - AN(without CC), identical in both branches;
#   with CC, deaths = projected deaths + AN(with CC) - AN(without CC)
#   (Simon, 18 Sep 2026; methods draft 2.5).
#
################################################################################

source("pipeline/00_pkg_params.R")

message(sprintf("\n[03] Assembling %s master analysis table...", city_name))

master_file <- file.path(out_dir, "03_master_table.csv")
checks_file <- file.path(check_dir, "03_master_checks.csv")
failures_file <- file.path(check_dir, "03_master_failures.csv")
fig_file <- file.path(fig_dir, "03_master_table_diagnostic.png")

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

plot_totals <- master[, .(
  deaths = unique(death),
  temp_deaths = unique(temp_deaths),
  rest = sum(rest)
), by = .(branch, year, age)]
plot_totals <- plot_totals[, .(
  deaths = sum(deaths),
  temp_deaths = sum(temp_deaths),
  rest = sum(rest)
), by = .(branch, year)]

plot_dt <- melt(plot_totals, id.vars = c("branch", "year"), variable.name = "measure", value.name = "value")

p <- ggplot(plot_dt, aes(x = year, y = value, color = measure)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~branch, scales = "free_y") +
  labs(
    title = sprintf("%s master table diagnostic", city_name),
    subtitle = sprintf("Demographics, AN, and rest mortality; GCM=%s", gcm_name),
    x = "Year",
    y = "Count"
  ) +
  theme_minimal(base_size = 11)

ggsave(fig_file, p, width = 11, height = 6, dpi = 160)

message("Saved master table to ", master_file)
message("Saved checks to ", checks_file)
message("Saved diagnostic figure to ", fig_file)
