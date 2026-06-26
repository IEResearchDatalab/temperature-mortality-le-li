################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 01_PrepData.R
#   Load raw ERF inputs, city metadata, and raw Eurostat demographic tables
#   No mortality projection model is built here.
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")

#---------------------------
# Load ERF coefficients
#---------------------------

coefs_all <- fread(file.path(path_data, "coefs.csv"))
if (!is.null(city_subset)) {
  coefs_all <- coefs_all[URAU_CODE %in% city_subset]
}

city_codes <- sort(unique(coefs_all$URAU_CODE))
n_cities <- length(city_codes)

#---------------------------
# Load pre-simulated coefficients (for uncertainty)
#---------------------------

coef_simu <- fread(file.path(path_data, "coef_simu.csv"))
if (!is.null(city_subset)) {
  coef_simu <- coef_simu[URAU_CODE %in% city_subset]
}

#---------------------------
# Load city metadata and baseline demography
#---------------------------

city_res <- fread(file.path(path_data, "city_results.csv"))
if (!is.null(city_subset)) {
  city_res <- city_res[URAU_CODE %in% city_subset]
}

# Keep only cities present in both ERF and city metadata inputs
common_cities <- intersect(city_codes, unique(city_res$URAU_CODE))
coefs_all <- coefs_all[URAU_CODE %in% common_cities]
coef_simu <- coef_simu[URAU_CODE %in% common_cities]
city_res  <- city_res[URAU_CODE %in% common_cities]
city_codes <- sort(common_cities)
n_cities <- length(city_codes)

#---------------------------
# Load raw Eurostat tables
#---------------------------
# Important: this script only loads raw Eurostat data.
# Any construction of future mortality schedules should be done later,
# in a dedicated script with explicit assumptions.

countries <- sort(unique(city_res$CNTR_CODE))

pop_eu_raw <- eurostat::get_eurostat(
  "proj_19np",
  time_format = "num",
  cache = TRUE,
  stringsAsFactors = FALSE
)
setDT(pop_eu_raw)

lt_eu_raw <- eurostat::get_eurostat(
  "demo_mlifetable",
  time_format = "num",
  cache = TRUE,
  stringsAsFactors = FALSE
)
setDT(lt_eu_raw)

#---------------------------
# Filter Eurostat projected population
#---------------------------

pop_eu <- copy(pop_eu_raw)
pop_eu <- pop_eu[
  projection == "BSL" &
    geo %in% countries &
    sex == "T" &
    grepl("^Y[0-9]+$", age)
]
pop_eu[, age := as.integer(sub("Y", "", age))]
pop_eu <- pop_eu[
  age >= min(single_ages) & age <= max(single_ages) + 1L &
    TIME_PERIOD >= proj_year_min & TIME_PERIOD <= proj_year_max
]
setnames(pop_eu, c("TIME_PERIOD", "values"), c("year", "pop"))
pop_eu <- unique(pop_eu[, .(geo, year, age, pop)], by = c("geo", "year", "age"))
setorder(pop_eu, geo, year, age)
setnames(pop_eu, "geo", "country")

#---------------------------
# Filter Eurostat historical life-table death rates
#---------------------------

lt_eu <- copy(lt_eu_raw)
lt_eu <- lt_eu[
  geo %in% countries &
    sex == "T" &
    indic_de == "DEATHRATE" &
    grepl("^Y[0-9]+$", age)
]
lt_eu[, age := as.integer(sub("Y", "", age))]
lt_eu <- lt_eu[
  age >= min(single_ages) & age <= max(single_ages) &
    TIME_PERIOD >= hist_year_min & TIME_PERIOD <= hist_year_max
]
setnames(lt_eu, c("TIME_PERIOD", "values"), c("year", "mx"))
lt_eu <- unique(lt_eu[, .(geo, year, age, mx)], by = c("geo", "year", "age"))
setorder(lt_eu, geo, year, age)
setnames(lt_eu, "geo", "country")

#---------------------------
# Basic validation checks
#---------------------------

missing_city_meta <- setdiff(unique(coefs_all$URAU_CODE), unique(city_res$URAU_CODE))
missing_coef_simu <- setdiff(unique(coefs_all$URAU_CODE), unique(coef_simu$URAU_CODE))
missing_countries_pop <- setdiff(countries, unique(pop_eu$country))
missing_countries_lt  <- setdiff(countries, unique(lt_eu$country))

if (length(missing_city_meta) > 0) {
  stop(sprintf("Missing city metadata for %d ERF cities", length(missing_city_meta)))
}
if (length(missing_coef_simu) > 0) {
  warning(sprintf("Missing simulated coefficients for %d ERF cities", length(missing_coef_simu)))
}
if (length(missing_countries_pop) > 0) {
  stop(sprintf("Missing Eurostat projected population for countries: %s",
               paste(missing_countries_pop, collapse = ", ")))
}
if (length(missing_countries_lt) > 0) {
  stop(sprintf("Missing Eurostat life-table data for countries: %s",
               paste(missing_countries_lt, collapse = ", ")))
}

#---------------------------
# Save prepared raw inputs
#---------------------------

saveRDS(coefs_all, file.path(path_out, "coefs_all.rds"))
saveRDS(coef_simu, file.path(path_out, "coef_simu.rds"))
saveRDS(city_res, file.path(path_out, "city_res.rds"))
saveRDS(pop_eu, file.path(path_out, "pop_eu_raw_filtered.rds"))
saveRDS(lt_eu, file.path(path_out, "lt_eu_raw_filtered.rds"))

cat(sprintf("Cities retained: %d\n", n_cities))
cat(sprintf("Countries retained: %d\n", length(countries)))
cat(sprintf("ERF rows: %d | Simulated-coef rows: %d | City metadata rows: %d\n",
            nrow(coefs_all), nrow(coef_simu), nrow(city_res)))
cat(sprintf("Projected population rows: %d | Historical life-table rows: %d\n",
            nrow(pop_eu), nrow(lt_eu)))
