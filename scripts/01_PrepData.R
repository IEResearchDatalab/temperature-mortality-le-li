################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 01_PrepData.R
#   Load ERF coefficients, population, and mortality projections
#
################################################################################

if (length(ls()) == 0) source("00_Packages_Parameters.R")

#---------------------------
# Load ERF coefficients
#---------------------------

coefs_all <- fread(file.path(path_data, "coefs.csv"))

if (!is.null(city_subset))
  coefs_all <- coefs_all[URAU_CODE %in% city_subset]

city_codes <- unique(coefs_all$URAU_CODE)
n_cities <- length(city_codes)

cat(sprintf("Cities: %d | Age groups: %d\n", n_cities, length(agelabs)))

#---------------------------
# Load pre-simulated coefficients (for CIs)
#---------------------------

coef_simu <- fread(file.path(path_data, "coef_simu.csv"))
if (!is.null(city_subset))
  coef_simu <- coef_simu[URAU_CODE %in% city_subset]

#---------------------------
# Load city metadata and demography
#---------------------------

city_res <- fread(file.path(path_data, "city_results.csv"))
if (!is.null(city_subset))
  city_res <- city_res[URAU_CODE %in% city_subset]

#---------------------------
# Load national mortality projections (downloaded on the fly from Eurostat)
#---------------------------

source("../R/load_data.R")

# Download Eurostat data once, process per country
pop_eu <- get_eurostat("proj_19np", time_format = "num", cache = TRUE,
  stringsAsFactors = FALSE)
lt_eu <- get_eurostat("demo_mlifetable", time_format = "num", cache = TRUE,
  stringsAsFactors = FALSE)

countries <- unique(city_res$CNTR_CODE)
all_mort <- rbindlist(lapply(countries, function(cc) {
  tryCatch({
    load_eurostat_mortality(cc, sex = "T", year_min = proj_year_min,
      year_max = proj_year_max, pop_cache = pop_eu, lt_cache = lt_eu)
  }, error = function(e) NULL)
}))

cat(sprintf("Eurostat mortality: %d countries, %d rows\n",
  length(unique(all_mort$age)), nrow(all_mort)))