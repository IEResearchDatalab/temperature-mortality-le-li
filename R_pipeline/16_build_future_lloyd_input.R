#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

message("\n[16] Building Lloyd-compatible future life-table input...")

input_file <- trimws(Sys.getenv("INPUT_FILE", unset = "results/le_li_input/lloyd_fig_s1_future_city_year_single_age.csv"))
output_wide_file <- trimws(Sys.getenv("OUTPUT_WIDE_FILE", unset = "results/le_li_input/future_temperature_deaths_wide.csv"))
output_long_file <- trimws(Sys.getenv("OUTPUT_LONG_FILE", unset = "results/le_li_input/future_deaths_input_100.csv"))
check_file <- trimws(Sys.getenv("CHECK_FILE", unset = "results/le_li_input/future_deaths_input_100_checks.csv"))
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

dt <- fread(input_file)

required_cols <- c(
  "URAU_CODE", "CNTR_CODE", "cntr_name", "year", "ssp", "gcm", "sim", "age",
  "pop", "death_baseline", "AN_ExtrCold", "AN_ModCold", "AN_ModHeat", "AN_ExtrHeat"
)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols)) {
  stop(sprintf("Input file is missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

if (!is.null(wanted_countries)) dt <- dt[CNTR_CODE %in% wanted_countries]
if (!is.null(wanted_years)) dt <- dt[year %in% wanted_years]
if (!is.null(wanted_ssps)) dt <- dt[ssp %in% wanted_ssps]
if (!is.null(wanted_gcms)) dt <- dt[gcm %in% wanted_gcms]
if (!is.null(wanted_sims)) dt <- dt[sim %in% wanted_sims]

if (!nrow(dt)) {
  stop("All selected rows were filtered out; nothing to write.", call. = FALSE)
}

if (geo_level == "country") {
  dt[, `:=`(geo_id = CNTR_CODE, geo_name = cntr_name)]
} else if (geo_level == "europe") {
  dt[, `:=`(geo_id = "EUROPE", geo_name = "Europe")]
} else {
  stop("GEO_LEVEL must be either 'country' or 'europe'.", call. = FALSE)
}

wide <- dt[, .(
  pop = sum(pop),
  total_deaths = sum(death_baseline),
  extr_cold = sum(AN_ExtrCold),
  mod_cold = sum(AN_ModCold),
  mod_heat = sum(AN_ModHeat),
  extr_heat = sum(AN_ExtrHeat)
), by = .(geo_id, geo_name, year, ssp, gcm, sim, age)]

wide[, rest := total_deaths - (extr_cold + mod_cold + mod_heat + extr_heat)]

# Remove tiny floating-point artefacts (e.g. -1e-8 deaths)
floating_point_tol <- 1e-7
wide[abs(rest) < floating_point_tol, rest := 0]

wide[, `:=`(sex = sex_label, reth = reth_label)]
setcolorder(wide, c(
  "geo_id", "geo_name", "year", "ssp", "gcm", "sim", "sex", "reth", "age", "pop",
  "total_deaths", "extr_cold", "mod_cold", "mod_heat", "extr_heat", "rest"
))

long <- melt(
  wide,
  id.vars = c("geo_id", "geo_name", "year", "ssp", "gcm", "sim", "sex", "reth", "age", "pop", "total_deaths"),
  measure.vars = c("extr_cold", "mod_cold", "mod_heat", "extr_heat", "rest"),
  variable.name = "cause",
  value.name = "deaths_cause"
)
long[, mx_cause := deaths_cause / pop]
long[, mx_total := total_deaths / pop]

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
), by = .(geo_id, geo_name, year, ssp, gcm, sim), .SDcols = c("pop", "total_deaths", "extr_cold", "mod_cold", "mod_heat", "extr_heat", "rest")]

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