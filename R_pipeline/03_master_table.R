#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

message("\n[03] Assembling Madrid master analysis table...")

out_dir <- "results/phase1_madrid"
check_dir <- "results/checks"
fig_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

master_file <- file.path(out_dir, "03_master_table.csv")
checks_file <- file.path(check_dir, "03_master_checks.csv")
failures_file <- file.path(check_dir, "03_master_failures.csv")
fig_file <- file.path(fig_dir, "03_master_table_diagnostic.png")

city_id <- "ES001C"
city_name <- "Madrid"
ssp_target <- 3L
gcm_target <- "GFDL_ESM4"
branch_levels <- c("with_cc", "without_cc")
range_levels <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
age_levels <- 65:100

grouped_dem <- fread(file.path(out_dir, "00_demography_grouped.csv"))
single_dem <- fread(file.path(out_dir, "00_demography_single_age.csv"))
single_an <- fread(file.path(out_dir, "02_single_age_an.csv"))

grouped_dem <- grouped_dem[geo_id == city_id & ssp == ssp_target]
single_dem <- single_dem[geo_id == city_id & ssp == ssp_target]
single_an <- single_an[geo_id == city_id & ssp == ssp_target & gcm == gcm_target]

if (!nrow(grouped_dem)) stop("Grouped demography missing for Madrid.", call. = FALSE)
if (!nrow(single_dem)) stop("Single-age demography missing for Madrid.", call. = FALSE)
if (!nrow(single_an)) stop("Single-age AN missing for Madrid.", call. = FALSE)

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

dem_grid <- CJ(year = sort(unique(single_dem$year)), age = age_levels)
an_grid <- CJ(year = sort(unique(single_an$year)), age = age_levels, range = range_levels, branch = branch_levels)

if (nrow(unique(single_dem, by = c("year", "age"))) != nrow(dem_grid)) {
  stop("Single-age demographic primary keys are incomplete.", call. = FALSE)
}
if (nrow(unique(single_an, by = c("year", "branch", "age", "range"))) != nrow(an_grid)) {
  stop("Single-age AN primary keys are incomplete.", call. = FALSE)
}

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

master[, age_temp_deaths := sum(an), by = .(year, branch, age)]
master[, rest := death - age_temp_deaths]

master[, `:=`(
  geo_id = city_id,
  label = city_name,
  ssp = ssp_target,
  gcm = gcm_target,
  source_agegroup = agegroup,
  temp_deaths = age_temp_deaths
)]

setcolorder(master, c(
  "geo_id", "label", "ssp", "gcm", "branch", "year", "age", "agegroup", "source_agegroup",
  "range", "pop", "death", "grouped_pop", "grouped_death", "country_pop", "country_death",
  "pop_share", "death_share", "pop_weight", "death_weight", "group_an", "weight", "an",
  "temp_deaths", "rest"
))

master[, agegroup := as.character(agegroup)]
master[, source_agegroup := as.character(source_agegroup)]

full_grid <- CJ(branch = branch_levels, year = sort(unique(master$year)), age = age_levels, range = range_levels)

unique_keys <- unique(master[, .(geo_id, label, ssp, gcm, branch, year, age, range)])
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

required_cols <- c("pop", "death", "grouped_pop", "grouped_death", "an", "temp_deaths", "rest")
na_rows <- master[!complete.cases(master[, ..required_cols])]

checks <- data.table(
  check_name = c(
    "primary_keys_unique",
    "complete_grid",
    "required_columns_finite",
    "rest_mortality_nonnegative",
    "demographics_identical_across_branches"
  ),
  status = c(
    if (nrow(unique_keys) == nrow(master)) "PASS" else "FAIL",
    if (nrow(unique(master, by = c("branch", "year", "age", "range"))) == nrow(full_grid)) "PASS" else "FAIL",
    if (nrow(na_rows) == 0L) "PASS" else "FAIL",
    if (min(master$rest) >= -1e-9) "PASS" else "FAIL",
    if (max(branch_delta$pop_delta) == 0 && max(branch_delta$death_delta) == 0 && max(branch_delta$grouped_pop_delta) == 0 && max(branch_delta$grouped_death_delta) == 0) "PASS" else "FAIL"
  ),
  value = c(
    nrow(unique_keys),
    nrow(full_grid),
    nrow(na_rows),
    sprintf("min_rest=%g", min(master$rest)),
    sprintf("max_pop_delta=%g; max_death_delta=%g", max(branch_delta$pop_delta), max(branch_delta$death_delta))
  ),
  threshold = c(
    sprintf("%d rows", nrow(master)),
    sprintf("%d rows", nrow(full_grid)),
    "0 rows with NA/NaN/Inf",
    ">= -1e-9",
    "zero difference across branches"
  )
)

failures <- data.table()
if (any(checks$status == "FAIL")) {
  failures <- rbindlist(list(
    if (checks$status[checks$check_name == "primary_keys_unique"] == "FAIL") {
      unique_keys[duplicated(unique_keys) | duplicated(unique_keys, fromLast = TRUE), .(
        geo_id, label, ssp, gcm, branch, year, age, range,
        failing_check = "primary_keys_unique",
        observed_value = "duplicate row",
        expected_bound = "unique keys"
      )]
    } else NULL,
    if (checks$status[checks$check_name == "complete_grid"] == "FAIL") {
      merged_grid <- merge(full_grid, unique(master[, .(branch, year, age, range)]), by = c("branch", "year", "age", "range"), all.x = TRUE, sort = FALSE)
      merged_grid[is.na(V1), .(
        branch, year, age, range,
        failing_check = "complete_grid",
        observed_value = "missing row",
        expected_bound = "present"
      )]
    } else NULL,
    if (checks$status[checks$check_name == "required_columns_finite"] == "FAIL") {
      na_rows[, .(
        geo_id, label, ssp, gcm, branch, year, age, range,
        failing_check = "required_columns_finite",
        observed_value = "NA/NaN/Inf present",
        expected_bound = "all required fields finite"
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
    } else NULL
  ), fill = TRUE)
}

fwrite(master, master_file)
fwrite(checks, checks_file)
if (nrow(failures)) {
  fwrite(failures, failures_file)
} else {
  if (file.exists(failures_file)) file.remove(failures_file)
  invisible(file.create(failures_file))
}

plot_totals <- master[, .(
  deaths = unique(death),
  temp_deaths = unique(temp_deaths),
  rest = unique(rest)
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
    title = "Madrid master table diagnostic",
    subtitle = sprintf("Demographics, AN, and rest mortality; GCM=%s", gcm_target),
    x = "Year",
    y = "Count"
  ) +
  theme_minimal(base_size = 11)

ggsave(fig_file, p, width = 11, height = 6, dpi = 160)

if (nrow(failures)) {
  stop(sprintf("03_master_table.R failed %d invariant(s); see %s", nrow(failures), failures_file), call. = FALSE)
}

message("Saved master table to ", master_file)
message("Saved checks to ", checks_file)
message("Saved diagnostic figure to ", fig_file)
