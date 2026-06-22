################################################################################
#
# Excess mortality attributed to heat and cold: 
#   a health impact assessment study in 854 cities in Europe
#
# The Lancet Planetary Health, 2023
# https://doi.org/10.1016/S2542-5196(23)00023-2
#
# Data loading functions
#
################################################################################

library(data.table)
library(arrow)
library(dplyr)
library(eurostat)

#---------------------------
# Load projected temperatures
#---------------------------

load_projected_temperatures <- function(city_code,
  parquet_path = "data/tmeanproj.gz.parquet",
  gcmexcl = character(0)) {

  proj_data <- open_dataset(parquet_path) %>%
    dplyr::filter(URAU_CODE == city_code) %>%
    dplyr::collect() %>%
    as.data.table()

  proj_data[, year := year(date)]
  proj_data[, doy := as.integer(format(date, "%j"))]
  proj_data[doy > 365, doy := 365L]

  gcm_cols <- names(proj_data)[grepl("^tas_", names(proj_data))]
  gcm_cols <- gcm_cols[!gsub("tas_", "", gcm_cols) %in% gcmexcl]

  list(proj_data = proj_data, gcm_cols = gcm_cols)
}

#---------------------------
# Extract historical temperatures
#---------------------------

extract_hist_temps <- function(proj_data, gcm_cols) {
  hist_data <- proj_data[ssp == "hist"]
  hist_temps <- unlist(hist_data[, ..gcm_cols], use.names = FALSE)
  hist_temps[!is.na(hist_temps)]
}

#---------------------------
# Pool baseline temperatures
#---------------------------

pool_baseline_temperatures <- function(proj_data, gcm_cols, ssp_codes,
  baseline_temp_period) {

  baseline_hist <- proj_data[ssp == "hist" & year %in% baseline_temp_period]
  baseline_proj <- proj_data[ssp %in% ssp_codes &
    year %in% baseline_temp_period & year > 2014]

  temps <- c(
    unlist(baseline_hist[, ..gcm_cols], use.names = FALSE),
    unlist(baseline_proj[, ..gcm_cols], use.names = FALSE))

  doys <- c(
    rep(baseline_hist$doy, length(gcm_cols)),
    rep(baseline_proj$doy, length(gcm_cols)))

  valid <- !is.na(temps)
  list(temps = temps[valid], doys = doys[valid])
}

#---------------------------
# Load seasonal weights
#---------------------------

load_seasonal_weights <- function(city_name_lower) {
  sw_file <- sprintf("results_csv/seasonal_weights_daily_%s.csv",
    city_name_lower)

  if (!file.exists(sw_file))
    return(list(sw_matrix = NULL, available = FALSE))

  sw_dt <- fread(sw_file)
  sw_matrix <- matrix(1 / 365, nrow = 81, ncol = 365,
    dimnames = list(20:100, 1:365))

  for (i in seq_len(nrow(sw_dt))) {
    a <- sw_dt$age[i]
    d <- sw_dt$doy[i]
    sw_matrix[as.character(a), d] <- sw_dt$weight[i]
  }

  list(sw_matrix = sw_matrix, available = TRUE)
}

#---------------------------
# Load projected mortality from Eurostat (EUROPOP)
#---------------------------

load_eurostat_mortality <- function(country_code, sex = "M",
  age_min = 20L, age_max = 100L, year_min = 2020L, year_max = 2100L,
  baseline_year = 2022L, pop_cache = NULL, lt_cache = NULL) {

  if (is.null(pop_cache)) {
    pop <- get_eurostat("proj_19np", time_format = "num", cache = TRUE,
      stringsAsFactors = FALSE)
    setDT(pop)
  } else {
    pop <- copy(pop_cache)
    if (!is.data.table(pop)) setDT(pop)
  }
  pop <- pop[projection == "BSL"]

  pop <- pop[geo == country_code & sex == sex]
  pop <- pop[grepl("^Y[0-9]+$", age)]
  pop[, age := as.integer(sub("Y", "", age))]
  pop <- pop[age >= age_min & age <= age_max + 1L &
    TIME_PERIOD >= year_min & TIME_PERIOD <= year_max]
  setnames(pop, c("TIME_PERIOD", "values"), c("year", "pop"))
  pop <- pop[, .(year, age, pop)]
  setorder(pop, year, age)

  # Cohort-based survival ratio: s(a,t) = P(a+1, t+1) / P(a, t)
  pop[, cohort := year - age]
  setorder(pop, cohort, year, age)
  pop[, `:=`(pop_next = shift(pop, type = "lead"),
    year_next = shift(year, type = "lead"),
    age_next = shift(age, type = "lead"))]
  surv <- pop[age_next == age + 1L & year_next == year + 1L]
  surv <- surv[pop > 0 & pop_next > 0]
  surv[, s := pop_next / pop]
  surv <- surv[s > 0 & s <= 1]
  surv[, mx_trend := -log(s)]
  surv <- unique(surv[, .(age, year, mx_trend)], by = c("age", "year"))

  # Deduplicate (should not happen, but safeguard)
  surv <- unique(surv, by = c("age", "year"))

  # Trend ratio: future mx_trend / baseline mx_trend
  base_trend <- surv[year == baseline_year, .(age, mx_trend_base = mx_trend)]
  surv <- merge(surv, base_trend, by = "age")
  surv[, mx_ratio := mx_trend / mx_trend_base]

  # Historical mx from life tables
  if (is.null(lt_cache)) {
    lt <- get_eurostat("demo_mlifetable", time_format = "num", cache = TRUE,
      stringsAsFactors = FALSE)
    setDT(lt)
  } else {
    lt <- copy(lt_cache)
    if (!is.data.table(lt)) setDT(lt)
  }
  lt <- lt[geo == country_code & sex == sex & indic_de == "DEATHRATE"]
  lt <- lt[grepl("^Y[0-9]+$", age)]
  lt[, age := as.integer(sub("Y", "", age))]
  lt <- lt[age >= age_min & age <= age_max]
  setnames(lt, c("TIME_PERIOD", "values"), c("year", "mx_hist"))
  lt <- lt[, .(year, age, mx_hist)]
  lt <- unique(lt, by = c("year", "age"))
  setorder(lt, year, age)

  # Latest historical mx as baseline (most recent year with data)
  lt_last <- lt[year == max(year), .(age, mx_base = mx_hist)]
  # Extrapolate mx_base for ages beyond available data (forward fill)
  all_ages <- data.table(age = age_min:age_max)
  lt_last <- merge(all_ages, lt_last, by = "age", all.x = TRUE)
  lt_last[, mx_base := nafill(mx_base, type = "locf")]
  lt_last[mx_base <= 0 | is.na(mx_base), mx_base := 1e-10]

  # Build full age-year grid with projected mx = mx_base * mx_ratio
  age_year_grid <- CJ(age = age_min:age_max, year = year_min:year_max)
  result <- merge(age_year_grid, lt_last, by = "age", allow.cartesian = TRUE)
  result <- merge(result, unique(surv[, .(age, year, mx_ratio)]),
    by = c("age", "year"), all.x = TRUE, allow.cartesian = TRUE)
  result[, mx := mx_base * mx_ratio]
  result <- unique(result[, .(age, year, mx, mx_ratio, mx_base)],
    by = c("age", "year"))

  # Forward-fill mx_ratio across ages for missing oldest ages
  setorder(result, year, age)
  result[, mx_ratio := nafill(mx_ratio, type = "locf"), by = year]

  # Ensure all age-year combinations are covered
  full <- CJ(age = age_min:age_max, year = year_min:year_max)
  result <- merge(full, result, by = c("age", "year"), all.x = TRUE)
  result[, mx_ratio := nafill(mx_ratio, type = "locf"), by = age]
  result[, mx_ratio := nafill(mx_ratio, type = "nocb"), by = age]
  result[mx_ratio <= 0 | is.infinite(mx_ratio), mx_ratio := 1]
  result[, mx_base := nafill(mx_base, type = "locf"), by = age]
  result[is.na(mx), mx := mx_base * mx_ratio]

  # Interpolate mx_ratio for any missing years (edge cases)
  result[, mx_ratio := nafill(mx_ratio, type = "locf"), by = age]
  result[, mx_ratio := nafill(mx_ratio, type = "nocb"), by = age]
  result[mx_ratio <= 0 | is.infinite(mx_ratio), mx_ratio := 1]

  result[, mx := mx_base * mx_ratio]
  result[, .(age, year, mx, mx_ratio, mx_base)]
}