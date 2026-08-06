#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

message("\n[16] Building Lloyd-compatible future life-table input with fixed baseline rest...")
message("NOTE: baseline rest is temporarily clamped at ages with negative values.")
message("      Final analyses will use a constrained baseline redistribution.")

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

projected_demography_file <- trimws(Sys.getenv(
  "FUTURE_DEMOGRAPHIC_FILE",
  unset = ""
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

read_projected_demography <- function(path) {
  if (!nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    stop(sprintf("Projected demographic input not found: %s", path), call. = FALSE)
  }

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    dem <- readRDS(path)
  } else {
    dem <- fread(path)
  }

  if (!is.data.table(dem)) {
    dem <- as.data.table(dem)
  }

  if ("CNTR_CODE" %in% names(dem)) {
    dem[, geo_id := as.character(CNTR_CODE)]
  } else if ("geo_id" %in% names(dem)) {
    dem[, geo_id := as.character(geo_id)]
  } else if ("country" %in% names(dem)) {
    dem[, geo_id := as.character(country)]
  } else {
    stop("Projected demographic input must contain CNTR_CODE, geo_id, or country.", call. = FALSE)
  }

  if ("year" %in% names(dem)) {
    dem[, year := as.integer(year)]
  } else if ("year5" %in% names(dem)) {
    dem[, year := as.integer(year5)]
  } else {
    stop("Projected demographic input must contain a year column.", call. = FALSE)
  }

  if ("age" %in% names(dem)) {
    dem[, age := as.integer(age)]
  } else if ("agegroup" %in% names(dem)) {
    dem[, age := NA_integer_]
  } else {
    stop("Projected demographic input must contain an age column.", call. = FALSE)
  }

  if (all(is.na(dem$age))) {
    stop(
      paste(
        "Projected demographic input does not contain single-age values.",
        "Script 16 expects age-specific deaths for ages 65:100.",
        "Please transform age-group rows into single-age rows before using FUTURE_DEMOGRAPHIC_FILE.",
        sep = " "
      ),
      call. = FALSE
    )
  }

  if ("death" %in% names(dem)) {
    dem[, death := as.numeric(death)]
  } else {
    stop("Projected demographic input must contain a death column.", call. = FALSE)
  }

  if ("pop" %in% names(dem)) {
    dem[, pop := as.numeric(pop)]
  } else {
    stop("Projected demographic input must contain a pop column.", call. = FALSE)
  }

  if (!is.null(wanted_countries)) {
    dem <- dem[geo_id %in% wanted_countries]
  }
  if (!is.null(wanted_years)) dem <- dem[year %in% wanted_years]
  if ("ssp" %in% names(dem) && !is.null(wanted_ssps)) dem <- dem[ssp %in% wanted_ssps]
  if ("gcm" %in% names(dem) && !is.null(wanted_gcms)) dem <- dem[gcm %in% wanted_gcms]
  if ("sim" %in% names(dem) && !is.null(wanted_sims)) dem <- dem[sim %in% wanted_sims]

  if (!nrow(dem)) {
    stop("Projected demographic input filtered to zero rows.", call. = FALSE)
  }

  dem
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
# Build the future rest reference
# ------------------------------------------------------------------------------

future_geo_year_age <- future_dt[, .(
  extr_cold = sum(AN_ExtrCold),
  mod_cold = sum(AN_ModCold),
  mod_heat = sum(AN_ModHeat),
  extr_heat = sum(AN_ExtrHeat)
), by = .(geo_id, geo_name, year, ssp, gcm, sim, age)]

future_geo_year_age[, temp_deaths := extr_cold + mod_cold + mod_heat + extr_heat]

projected_demography <- read_projected_demography(projected_demography_file)

if (!is.null(projected_demography)) {
  message("Using projected demographic deaths to build the future rest reference.")

  projected_cols <- c("geo_id", "year", "age")
  join_cols <- projected_cols
  if ("ssp" %in% names(projected_demography)) join_cols <- c(join_cols, "ssp")
  if ("gcm" %in% names(projected_demography)) join_cols <- c(join_cols, "gcm")
  if ("sim" %in% names(projected_demography)) join_cols <- c(join_cols, "sim")

  if (!all(join_cols %in% names(future_geo_year_age))) {
    future_geo_year_age[, `:=`(
      ssp = as.integer(ssp),
      gcm = as.character(gcm),
      sim = as.integer(sim)
    )]
  }

  projected_demography[, `:=`(
    geo_id = as.character(geo_id),
    age = as.integer(age),
    year = as.integer(year)
  )]

  join_cols <- intersect(join_cols, names(projected_demography))
  join_cols <- intersect(join_cols, names(future_geo_year_age))

  dem_keep_cols <- unique(c(join_cols, "death", "pop"))

  rest_ref <- merge(
    future_geo_year_age,
    projected_demography[, ..dem_keep_cols],
    by = join_cols,
    all.x = TRUE,
    sort = FALSE
  )

  if (anyNA(rest_ref$death)) {
    missing_keys <- rest_ref[is.na(death), unique(.(
      geo_id, year, age
    ))]
    stop(
      sprintf("Projected demographic deaths are missing for %d geo-year-age combinations.", nrow(missing_keys)),
      call. = FALSE
    )
  }

  rest_ref[, rest := death - temp_deaths]
  rest_ref <- rest_ref[, .(rest = mean(rest)), by = .(geo_id, geo_name, age)]
} else {
  message("Projected demographic file not supplied; falling back to the baseline rest reference.")

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
}

if (anyNA(rest_ref$rest)) {
  stop("Rest reference contains NA values.", call. = FALSE)
}

# ------------------------------------------------------------------
# TEMPORARY VALIDATION WORKAROUND
#
# Negative rest values occur because the baseline single-age
# attributable deaths were generated with the original unconstrained
# PCLM redistribution, whereas future attributable deaths use the
# constrained redistribution (Script 17).
#
# Until Script 11 is rebuilt using the constrained algorithm,
# negative baseline rest values are truncated to zero solely to
# validate the end-to-end LE/LI pipeline.
# ------------------------------------------------------------------

negative_rows <- copy(rest_ref[rest < -floating_point_tol])

if (nrow(negative_rows) > 0L) {
  clamped_file <- sub("\\.csv$", "_clamped_rows.csv", rest_file)

  message(sprintf(
    "  Clamping %d geo-age rows (minimum rest = %.2f deaths).",
    nrow(negative_rows),
    min(negative_rows$rest)
  ))

  fwrite(negative_rows, clamped_file)
  message("  Saved clamped rows to ", clamped_file)

  rest_ref[rest < -floating_point_tol, rest := 0]
}

rest_ref[abs(rest) < floating_point_tol, rest := 0]

fwrite(rest_ref, rest_file)
message("Saved future rest reference to ", rest_file)

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