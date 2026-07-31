#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

message("\n[16] Building Lloyd-compatible future life-table input with fixed baseline rest...")

# ------------------------------------------------------------------------------
# Inputs and outputs
# ------------------------------------------------------------------------------

future_file <- trimws(Sys.getenv(
  "INPUT_FILE",
  unset = "results/le_li_input/lloyd_fig_s1_future_city_year_single_age.csv"
))

baseline_file <- trimws(Sys.getenv(
  "BASELINE_FILE",
  unset = "results/le_li_input/lloyd_fig_s1_baseline_city_year_single_age.csv"
))

output_wide_file <- trimws(Sys.getenv(
  "OUTPUT_WIDE_FILE",
  unset = "results/le_li_input/future_temperature_deaths_wide.csv"
))

output_long_file <- trimws(Sys.getenv(
  "OUTPUT_LONG_FILE",
  unset = "results/le_li_input/future_deaths_input_100.csv"
))

check_file <- trimws(Sys.getenv(
  "CHECK_FILE",
  unset = "results/le_li_input/future_deaths_input_100_checks.csv"
))

rest_file <- trimws(Sys.getenv(
  "REST_FILE",
  unset = "results/le_li_input/fixed_rest_reference.csv"
))

geo_level <- trimws(Sys.getenv("GEO_LEVEL", unset = "country"))
country_filter <- trimws(Sys.getenv("COUNTRY_FILTER", unset = ""))
ssp_filter <- trimws(Sys.getenv("SSP_FILTER", unset = ""))
gcm_filter <- trimws(Sys.getenv("GCM_FILTER", unset = ""))
sim_filter <- trimws(Sys.getenv("SIM_FILTER", unset = ""))
year_filter <- trimws(Sys.getenv("YEAR_FILTER", unset = ""))
sex_label <- trimws(Sys.getenv("SEX_LABEL", unset = "pooled"))
reth_label <- trimws(Sys.getenv("RETH_LABEL", unset = "all"))

dir.create(dirname(output_wide_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_long_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(check_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(rest_file), recursive = TRUE, showWarnings = FALSE)

parse_int_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

parse_chr_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

wanted_years <- parse_int_filter(year_filter)
wanted_ssps <- parse_int_filter(ssp_filter)
wanted_sims <- parse_int_filter(sim_filter)
wanted_gcms <- parse_chr_filter(gcm_filter)
wanted_countries <- parse_chr_filter(country_filter)

floating_point_tol <- 1e-7

# ------------------------------------------------------------------------------
# Read future and baseline inputs
# ------------------------------------------------------------------------------

future_dt <- fread(future_file)
baseline_dt <- fread(baseline_file)

future_required_cols <- c(
  "URAU_CODE", "CNTR_CODE", "cntr_name", "year", "ssp", "gcm", "sim", "age",
  "pop", "AN_ExtrCold", "AN_ModCold", "AN_ModHeat", "AN_ExtrHeat"
)

baseline_required_cols <- c(
  "URAU_CODE", "CNTR_CODE", "cntr_name", "year", "age",
  "death_baseline", "AN_ExtrCold", "AN_ModCold", "AN_ModHeat", "AN_ExtrHeat"
)

missing_future <- setdiff(future_required_cols, names(future_dt))
if (length(missing_future)) {
  stop(
    sprintf("Future input is missing required columns: %s", paste(missing_future, collapse = ", ")),
    call. = FALSE
  )
}

missing_baseline <- setdiff(baseline_required_cols, names(baseline_dt))
if (length(missing_baseline)) {
  stop(
    sprintf("Baseline input is missing required columns: %s", paste(missing_baseline, collapse = ", ")),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Apply filters
# ------------------------------------------------------------------------------

if (!is.null(wanted_countries)) {
  future_dt <- future_dt[CNTR_CODE %in% wanted_countries]
  baseline_dt <- baseline_dt[CNTR_CODE %in% wanted_countries]
}

if (!is.null(wanted_years)) future_dt <- future_dt[year %in% wanted_years]
if (!is.null(wanted_ssps)) future_dt <- future_dt[ssp %in% wanted_ssps]
if (!is.null(wanted_gcms)) future_dt <- future_dt[gcm %in% wanted_gcms]
if (!is.null(wanted_sims)) future_dt <- future_dt[sim %in% wanted_sims]

if (!nrow(future_dt)) {
  stop("All selected future rows were filtered out; nothing to write.", call. = FALSE)
}

if (!nrow(baseline_dt)) {
  stop("All selected baseline rows were filtered out; cannot build fixed rest reference.", call. = FALSE)
}

if (geo_level == "country") {
  future_dt[, `:=`(geo_id = CNTR_CODE, geo_name = cntr_name)]
  baseline_dt[, `:=`(geo_id = CNTR_CODE, geo_name = cntr_name)]
} else if (geo_level == "europe") {
  future_dt[, `:=`(geo_id = "EUROPE", geo_name = "Europe")]
  baseline_dt[, `:=`(geo_id = "EUROPE", geo_name = "Europe")]
} else {
  stop("GEO_LEVEL must be either 'country' or 'europe'.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Build fixed baseline rest reference
# ------------------------------------------------------------------------------

baseline_geo_year_age <- baseline_dt[, .(
  total_deaths = sum(death_baseline),
  extr_cold = sum(AN_ExtrCold),
  mod_cold = sum(AN_ModCold),
  mod_heat = sum(AN_ModHeat),
  extr_heat = sum(AN_ExtrHeat)
), by = .(geo_id, geo_name, year, age)]

baseline_geo_year_age[, temp_deaths := extr_cold + mod_cold + mod_heat + extr_heat]

rest_ref <- baseline_geo_year_age[, .(
  rest = mean(total_deaths - temp_deaths)
), by = .(geo_id, geo_name, age)]

if (anyNA(rest_ref$rest)) {
  stop("Baseline rest reference contains NA values.", call. = FALSE)
}

n_negative <- sum(rest_ref$rest < -floating_point_tol)
if (n_negative > 0L) {
  negative_rows <- rest_ref[rest < -floating_point_tol]
  message(sprintf(
    "  Clamping %d baseline geo-age rows with negative rest, mostly at age %s.",
    n_negative,
    paste(sort(unique(negative_rows$age)), collapse = ", ")
  ))
  rest_ref[rest < -floating_point_tol, rest := 0]
}

# Save clamped rows for reproducibility
clamped_file <- sub(
  "\\.csv$",
  "_clamped_rows.csv",
  rest_file
)

negative_rows <- copy(rest_ref[rest < -floating_point_tol])

message(sprintf(
  "  Clamping %d baseline geo-age rows with negative rest, mostly at age %s.",
  nrow(negative_rows),
  paste(sort(unique(negative_rows$age)), collapse = ", ")
))

fwrite(negative_rows, clamped_file)
message("  Saved clamped rows to ", clamped_file)

# Temporary workaround until baseline is rebuilt with constrained redistribution
rest_ref[rest < -floating_point_tol, rest := 0]

fwrite(rest_ref, rest_file)
message("Saved fixed baseline rest reference to ", rest_file)

# ------------------------------------------------------------------------------
# Aggregate future temperature causes and merge fixed rest
# ------------------------------------------------------------------------------

wide <- future_dt[, .(
  pop = sum(pop),
  extr_cold = sum(AN_ExtrCold),
  mod_cold = sum(AN_ModCold),
  mod_heat = sum(AN_ModHeat),
  extr_heat = sum(AN_ExtrHeat)
), by = .(geo_id, geo_name, year, ssp, gcm, sim, age)]

wide <- merge(
  wide,
  rest_ref[, .(geo_id, age, rest)],
  by = c("geo_id", "age"),
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(wide$rest)) {
  missing_keys <- wide[is.na(rest), unique(.(
    geo_id, age
  ))]
  stop(
    sprintf("Missing fixed rest reference for %d geo-age combinations.", nrow(missing_keys)),
    call. = FALSE
  )
}

wide[, total_deaths := extr_cold + mod_cold + mod_heat + extr_heat + rest]

wide[, `:=`(
  sex = sex_label,
  reth = reth_label
)]

setcolorder(wide, c(
  "geo_id", "geo_name", "year", "ssp", "gcm", "sim", "sex", "reth", "age", "pop",
  "total_deaths", "extr_cold", "mod_cold", "mod_heat", "extr_heat", "rest"
))

# ------------------------------------------------------------------------------
# Long format for Script 18
# ------------------------------------------------------------------------------

long <- melt(
  wide,
  id.vars = c("geo_id", "geo_name", "year", "ssp", "gcm", "sim", "sex", "reth", "age", "pop", "total_deaths"),
  measure.vars = c("extr_cold", "mod_cold", "mod_heat", "extr_heat", "rest"),
  variable.name = "cause",
  value.name = "deaths_cause"
)

long[, mx_cause := deaths_cause / pop]
long[, mx_total := total_deaths / pop]

# ------------------------------------------------------------------------------
# Checks
# ------------------------------------------------------------------------------

expected_ages <- 65:100

checks <- wide[, .(
  n_ages = uniqueN(age),
  min_age = min(age),
  max_age = max(age),
  missing_age_count = length(setdiff(expected_ages, age)),
  max_abs_accounting_diff = max(abs(total_deaths - (extr_cold + mod_cold + mod_heat + extr_heat + rest))),
  any_na = anyNA(.SD),
  min_rest = min(rest),
  max_rest = max(rest),
  negative_rest_rows = sum(rest < -floating_point_tol)
), by = .(geo_id, geo_name, year, ssp, gcm, sim), .SDcols = c(
  "pop", "total_deaths", "extr_cold", "mod_cold", "mod_heat", "extr_heat", "rest"
)]

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

fwrite(wide, output_wide_file)
fwrite(long, output_long_file)
fwrite(checks, check_file)

message("Saved Lloyd-compatible wide future table to ", output_wide_file)
message("Saved Lloyd-compatible long future table to ", output_long_file)
message("Saved future input checks to ", check_file)
message(sprintf(
  "Wide rows: %d | Long rows: %d | Geographies: %d | Years: %d | SSPs: %d | GCMs: %d | Sims: %d | Any NA checks: %s | Max abs accounting diff: %g",
  nrow(wide),
  nrow(long),
  uniqueN(wide$geo_id),
  uniqueN(wide$year),
  uniqueN(wide$ssp),
  uniqueN(wide$gcm),
  uniqueN(wide$sim),
  any(checks$any_na),
  max(checks$max_abs_accounting_diff)
))