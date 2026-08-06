#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

message("\n[20] Building projected-demography adapter for Script 16...")

projdata_file <- trimws(Sys.getenv(
  "PROJDATA_FILE",
  unset = "results/projdata/projdata_prototype.csv"
))

template_file <- trimws(Sys.getenv(
  "TEMPLATE_FILE",
  unset = "results/le_li_input/lloyd_fig_s1_future_constrained_city_year_single_age.csv"
))

baseline_file <- trimws(Sys.getenv(
  "BASELINE_FILE",
  unset = "results/le_li_input/lloyd_fig_s1_baseline_city_year_single_age.csv"
))

output_file <- trimws(Sys.getenv(
  "OUTPUT_FILE",
  unset = "results/le_li_input/projdata_adapter/future_demography_for_script16.csv"
))

check_dir <- trimws(Sys.getenv(
  "CHECK_DIR",
  unset = "results/le_li_input/projdata_adapter"
))

country_filter <- trimws(Sys.getenv("COUNTRY_FILTER", unset = ""))
ssp_filter <- trimws(Sys.getenv("SSP_FILTER", unset = ""))
gcm_filter <- trimws(Sys.getenv("GCM_FILTER", unset = ""))
sim_filter <- trimws(Sys.getenv("SIM_FILTER", unset = ""))
year_filter <- trimws(Sys.getenv("YEAR_FILTER", unset = ""))

max_year <- suppressWarnings(as.integer(trimws(Sys.getenv("MAX_YEAR", unset = "2099"))))
if (is.na(max_year)) max_year <- 2099L

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

parse_int_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

parse_chr_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

agegroup_to_ages <- function(agegroup) {
  if (agegroup == "65-74") return(65:74)
  if (agegroup == "75-84") return(75:84)
  if (agegroup == "85+") return(85:100)
  integer(0)
}

wanted_countries <- parse_chr_filter(country_filter)
wanted_ssps <- parse_int_filter(ssp_filter)
wanted_sims <- parse_int_filter(sim_filter)
wanted_gcms <- parse_chr_filter(gcm_filter)
wanted_years <- parse_int_filter(year_filter)

tol <- 1e-7

if (!file.exists(projdata_file)) {
  stop(sprintf("PROJDATA_FILE not found: %s", projdata_file), call. = FALSE)
}
if (!file.exists(template_file)) {
  stop(sprintf("TEMPLATE_FILE not found: %s", template_file), call. = FALSE)
}
if (!file.exists(baseline_file)) {
  stop(sprintf("BASELINE_FILE not found: %s", baseline_file), call. = FALSE)
}

proj <- fread(projdata_file)
tmpl <- fread(template_file)
base <- fread(baseline_file)

if ("year" %in% names(tmpl)) {
  tmpl_year_max <- suppressWarnings(max(as.integer(tmpl$year), na.rm = TRUE))
  if (is.finite(tmpl_year_max)) {
    max_year <- max(max_year, as.integer(tmpl_year_max))
  }
}

proj_required <- c("CNTR_CODE", "agegroup", "ssp", "year5", "pop", "death")
missing_proj <- setdiff(proj_required, names(proj))
if (length(missing_proj)) {
  stop(sprintf("PROJDATA_FILE missing columns: %s", paste(missing_proj, collapse = ", ")), call. = FALSE)
}

tmpl_required <- c("CNTR_CODE", "year", "ssp", "gcm", "sim", "age")
missing_tmpl <- setdiff(tmpl_required, names(tmpl))
if (length(missing_tmpl)) {
  stop(sprintf("TEMPLATE_FILE missing columns: %s", paste(missing_tmpl, collapse = ", ")), call. = FALSE)
}

base_required <- c("CNTR_CODE", "year", "age", "pop", "death_baseline")
missing_base <- setdiff(base_required, names(base))
if (length(missing_base)) {
  stop(sprintf("BASELINE_FILE missing columns: %s", paste(missing_base, collapse = ", ")), call. = FALSE)
}

proj[, `:=`(
  CNTR_CODE = as.character(CNTR_CODE),
  agegroup = as.character(agegroup),
  ssp = as.character(ssp),
  year5 = as.integer(year5),
  pop = as.numeric(pop),
  death = as.numeric(death)
)]

# Keep only future SSPs and 65+ groups for the Script 16 bridge.
proj <- proj[ssp %in% c("1", "2", "3") & agegroup %in% c("65-74", "75-84", "85+")]

if (!is.null(wanted_countries)) proj <- proj[CNTR_CODE %in% wanted_countries]
if (!is.null(wanted_ssps)) proj <- proj[ssp %in% as.character(wanted_ssps)]

if (!nrow(proj)) {
  stop("Projected demographic table is empty after SSP/country filtering.", call. = FALSE)
}

# City-level -> country-level aggregation.
country_group_5y <- proj[, .(
  pop = sum(pop),
  death = sum(death)
), by = .(CNTR_CODE, ssp, year5, agegroup)]

# Expand each 5-year point into annual years as a piecewise-constant annual schedule.
annual_list <- lapply(seq_len(nrow(country_group_5y)), function(i) {
  row <- country_group_5y[i]
  years <- row$year5:min(row$year5 + 4L, max_year)
  data.table(
    CNTR_CODE = row$CNTR_CODE,
    ssp = row$ssp,
    year5 = row$year5,
    year = years,
    agegroup = row$agegroup,
    pop = row$pop,
    death = row$death
  )
})
annual_group <- rbindlist(annual_list, use.names = TRUE)

if (!is.null(wanted_years)) annual_group <- annual_group[year %in% wanted_years]
if (!nrow(annual_group)) {
  stop("Annualized demographic table is empty after year filtering.", call. = FALSE)
}

# Build country-specific baseline age weights for 65+ age groups.
base <- base[age %between% c(65L, 100L)]
base[, `:=`(
  CNTR_CODE = as.character(CNTR_CODE),
  year = as.integer(year),
  age = as.integer(age),
  pop = as.numeric(pop),
  death_baseline = as.numeric(death_baseline)
)]

if (!is.null(wanted_countries)) base <- base[CNTR_CODE %in% wanted_countries]

base_country_age <- base[, .(
  pop = sum(pop),
  death = sum(death_baseline)
), by = .(CNTR_CODE, age)]

base_country_age[, agegroup := fifelse(
  age <= 74L, "65-74",
  fifelse(age <= 84L, "75-84", "85+")
)]

weights_country <- base_country_age[, {
  pop_den <- sum(pop)
  death_den <- sum(death)
  .(
    age = age,
    w_pop = if (pop_den > 0) pop / pop_den else NA_real_,
    w_death = if (death_den > 0) death / death_den else NA_real_
  )
}, by = .(CNTR_CODE, agegroup)]

# Global fallback weights for any country/group that lacks baseline support.
base_global_age <- base_country_age[, .(
  pop = sum(pop),
  death = sum(death)
), by = .(agegroup, age)]

weights_global <- base_global_age[, {
  pop_den <- sum(pop)
  death_den <- sum(death)
  .(
    age = age,
    w_pop_global = if (pop_den > 0) pop / pop_den else NA_real_,
    w_death_global = if (death_den > 0) death / death_den else NA_real_
  )
}, by = .(agegroup)]

# Prepare annual country-group rows for single-age expansion.
annual_group_exp <- copy(annual_group)
annual_group_exp[, row_id := .I]

age_map <- rbindlist(list(
  data.table(agegroup = "65-74", age = 65:74),
  data.table(agegroup = "75-84", age = 75:84),
  data.table(agegroup = "85+", age = 85:100)
))

annual_single <- merge(
  annual_group_exp,
  age_map,
  by = "agegroup",
  allow.cartesian = TRUE,
  sort = FALSE
)

annual_single <- merge(
  annual_single,
  weights_country,
  by = c("CNTR_CODE", "agegroup", "age"),
  all.x = TRUE,
  sort = FALSE
)

annual_single <- merge(
  annual_single,
  weights_global,
  by = c("agegroup", "age"),
  all.x = TRUE,
  sort = FALSE
)

annual_single[, `:=`(
  w_pop_final = fifelse(is.na(w_pop), w_pop_global, w_pop),
  w_death_final = fifelse(is.na(w_death), w_death_global, w_death)
)]

if (anyNA(annual_single$w_pop_final) || anyNA(annual_single$w_death_final)) {
  stop("Missing age-disaggregation weights after country + global fallback merge.", call. = FALSE)
}

annual_single[, `:=`(
  pop_single = pop * w_pop_final,
  death_single = death * w_death_final
)]

# Row-level preservation checks after age disaggregation.
age_preservation <- annual_single[, .(
  pop_group = unique(pop),
  pop_single_sum = sum(pop_single),
  death_group = unique(death),
  death_single_sum = sum(death_single),
  pop_abs_diff = abs(unique(pop) - sum(pop_single)),
  death_abs_diff = abs(unique(death) - sum(death_single))
), by = .(CNTR_CODE, ssp, year, year5, agegroup, row_id)]

max_pop_diff <- max(age_preservation$pop_abs_diff)
max_death_diff <- max(age_preservation$death_abs_diff)

if (max_pop_diff > tol || max_death_diff > tol) {
  stop(
    sprintf(
      "Age disaggregation failed preservation checks (max pop diff=%g, max death diff=%g).",
      max_pop_diff,
      max_death_diff
    ),
    call. = FALSE
  )
}

# Aggregate to the Script 16 geo/year/age/ssp space before scenario propagation.
ssp_age_schedule <- annual_single[, .(
  pop = sum(pop_single),
  death = sum(death_single)
), by = .(CNTR_CODE, ssp, year, age)]

# Scenario template from future AN table. We propagate SSP-specific demography
# across GCM/sim for matching country-year-SSP keys.
tmpl[, `:=`(
  CNTR_CODE = as.character(CNTR_CODE),
  ssp = as.character(ssp),
  year = as.integer(year),
  gcm = as.character(gcm),
  sim = as.integer(sim),
  age = as.integer(age)
)]

if (!is.null(wanted_countries)) tmpl <- tmpl[CNTR_CODE %in% wanted_countries]
if (!is.null(wanted_ssps)) tmpl <- tmpl[ssp %in% as.character(wanted_ssps)]
if (!is.null(wanted_gcms)) tmpl <- tmpl[gcm %in% wanted_gcms]
if (!is.null(wanted_sims)) tmpl <- tmpl[sim %in% wanted_sims]
if (!is.null(wanted_years)) tmpl <- tmpl[year %in% wanted_years]

tmpl_scenarios <- unique(tmpl[, .(CNTR_CODE, ssp, year, gcm, sim)])

if (!nrow(tmpl_scenarios)) {
  stop("Template scenario table is empty after filtering.", call. = FALSE)
}

adapter_out <- merge(
  tmpl_scenarios,
  ssp_age_schedule,
  by = c("CNTR_CODE", "ssp", "year"),
  all.x = TRUE,
  sort = FALSE,
  allow.cartesian = TRUE
)

coverage <- adapter_out[, .(
  n_rows = .N,
  missing_pop = sum(is.na(pop)),
  missing_death = sum(is.na(death))
), by = .(CNTR_CODE, ssp, year, gcm, sim)]

if (any(coverage$missing_pop > 0L | coverage$missing_death > 0L)) {
  bad <- coverage[missing_pop > 0L | missing_death > 0L]
  fwrite(bad, file.path(check_dir, "adapter_missing_coverage.csv"))
  stop(
    sprintf(
      "Adapter output has missing projected demography for %d scenario rows; see %s",
      nrow(bad),
      file.path(check_dir, "adapter_missing_coverage.csv")
    ),
    call. = FALSE
  )
}

adapter_out <- adapter_out[age %between% c(65L, 100L)]
adapter_out[, geo_id := CNTR_CODE]

setcolorder(adapter_out, c("geo_id", "CNTR_CODE", "year", "age", "ssp", "gcm", "sim", "pop", "death"))

# Validation: year expansion consistency (collapse annual back to year5 by mean).
year_roundtrip <- annual_group[, .(
  pop_roundtrip = mean(pop),
  death_roundtrip = mean(death)
), by = .(CNTR_CODE, ssp, year5, agegroup)]

year_roundtrip <- merge(
  year_roundtrip,
  country_group_5y,
  by = c("CNTR_CODE", "ssp", "year5", "agegroup"),
  suffixes = c("_annual", "_source")
)

year_roundtrip[, `:=`(
  pop_abs_diff = abs(pop_roundtrip - pop),
  death_abs_diff = abs(death_roundtrip - death)
)]

# Validation: scenario propagation consistency.
propagation_check <- adapter_out[, .(
  n_ages = uniqueN(age),
  min_age = min(age),
  max_age = max(age),
  total_pop = sum(pop),
  total_death = sum(death),
  any_negative_pop = any(pop < -tol),
  any_negative_death = any(death < -tol)
), by = .(CNTR_CODE, ssp, year, gcm, sim)]

summary_checks <- data.table(
  metric = c(
    "source_rows_65plus",
    "country_group_rows_5y",
    "annual_group_rows",
    "single_age_rows_before_scenario",
    "template_scenarios",
    "adapter_rows",
    "adapter_country_count",
    "adapter_ssp_count",
    "adapter_gcm_count",
    "adapter_sim_count",
    "adapter_year_min",
    "adapter_year_max",
    "age_preservation_max_pop_abs_diff",
    "age_preservation_max_death_abs_diff",
    "year_roundtrip_max_pop_abs_diff",
    "year_roundtrip_max_death_abs_diff",
    "adapter_negative_pop_rows",
    "adapter_negative_death_rows"
  ),
  value = c(
    nrow(proj),
    nrow(country_group_5y),
    nrow(annual_group),
    nrow(ssp_age_schedule),
    nrow(tmpl_scenarios),
    nrow(adapter_out),
    uniqueN(adapter_out$CNTR_CODE),
    uniqueN(adapter_out$ssp),
    uniqueN(adapter_out$gcm),
    uniqueN(adapter_out$sim),
    min(adapter_out$year),
    max(adapter_out$year),
    max_pop_diff,
    max_death_diff,
    max(year_roundtrip$pop_abs_diff),
    max(year_roundtrip$death_abs_diff),
    sum(adapter_out$pop < -tol),
    sum(adapter_out$death < -tol)
  )
)

fwrite(adapter_out, output_file)
fwrite(summary_checks, file.path(check_dir, "adapter_summary_checks.csv"))
fwrite(age_preservation, file.path(check_dir, "adapter_age_disaggregation_checks.csv"))
fwrite(year_roundtrip, file.path(check_dir, "adapter_year_roundtrip_checks.csv"))
fwrite(propagation_check, file.path(check_dir, "adapter_scenario_propagation_checks.csv"))

message("Saved adapter output to ", output_file)
message("Saved checks to ", check_dir)
message(sprintf(
  "Adapter rows: %d | Countries: %d | SSPs: %d | GCMs: %d | Sims: %d | Years: %d-%d",
  nrow(adapter_out), uniqueN(adapter_out$CNTR_CODE), uniqueN(adapter_out$ssp),
  uniqueN(adapter_out$gcm), uniqueN(adapter_out$sim),
  min(adapter_out$year), max(adapter_out$year)
))
