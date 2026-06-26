################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 09_Aggregation.R
#   Aggregate city-level outputs to country and Europe level.
#
# Inputs:
#   - analysis_dataset_all_cities.csv
#   - lifespan_inequality_all_cities.csv
#   - decomposition_period_all_cities.csv
#   - city_res.rds (for country mapping)
#
# Outputs:
#   - country_summary.csv
#   - europe_summary.csv
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")

#---------------------------
# Load inputs
#---------------------------

city_res <- readRDS(file.path(path_out, "city_res.rds"))
city_codes <- unique(city_res$URAU_CODE)

# City-to-country mapping
city_country <- unique(city_res[, .(city_code = URAU_CODE, country = CNTR_CODE)])

#---------------------------
# Aggregate analysis dataset
#---------------------------

analysis_all <- fread(file.path(path_out, "analysis_dataset_all_cities.csv"))

# Merge country code
analysis_all <- merge(analysis_all, city_country, by = "city_code", all.x = TRUE)

# Aggregate to country level: sum deaths and ANs by country-year-age-GCM-SSP
country_analysis <- analysis_all[,
  .(population = sum(population),
    deaths = sum(deaths),
    residual_deaths = sum(residual_deaths),
    extreme_cold = sum(extreme_cold),
    moderate_cold = sum(moderate_cold),
    moderate_heat = sum(moderate_heat),
    extreme_heat = sum(extreme_heat),
    n_cities = .N),
  by = .(country, year, age, gcm, ssp)]

# Aggregate to Europe level
europe_analysis <- analysis_all[,
  .(population = sum(population),
    deaths = sum(deaths),
    residual_deaths = sum(residual_deaths),
    extreme_cold = sum(extreme_cold),
    moderate_cold = sum(moderate_cold),
    moderate_heat = sum(moderate_heat),
    extreme_heat = sum(extreme_heat),
    n_cities = .N,
    n_countries = uniqueN(country)),
  by = .(year, age, gcm, ssp)]

europe_analysis[, country := "Europe"]

#---------------------------
# Aggregate life tables
#---------------------------

lt_all <- fread(file.path(path_out, "lifespan_inequality_all_cities.csv"))

# Merge country
lt_all <- merge(lt_all, city_country, by = "city_code", all.x = TRUE)

# Country: population-weighted average
city_pops <- analysis_all[, .(total_pop = sum(population)), by = .(city_code, year)]
lt_all <- merge(lt_all, city_pops, by = c("city_code", "year"), all.x = TRUE)
lt_all[is.na(total_pop), total_pop := 0]

country_lt <- lt_all[,
  .(e65 = weighted.mean(e65, total_pop, na.rm = TRUE),
    sd = weighted.mean(sd, total_pop, na.rm = TRUE),
    n_cities = .N),
  by = .(country, year, gcm, ssp)]

# Europe: population-weighted average
europe_lt <- lt_all[,
  .(e65 = weighted.mean(e65, total_pop, na.rm = TRUE),
    sd = weighted.mean(sd, total_pop, na.rm = TRUE),
    n_cities = .N,
    n_countries = uniqueN(country)),
  by = .(year, gcm, ssp)]

europe_lt[, country := "Europe"]

#---------------------------
# Aggregate decomposition
#---------------------------

decomp_period <- fread(file.path(path_out, "decomposition_period_all_cities.csv"))
decomp_period <- merge(decomp_period, city_country, by = "city_code", all.x = TRUE)

country_decomp <- decomp_period[,
  .(delta_e65 = sum(delta_e65),
    delta_sd = sum(delta_sd)),
  by = .(country, period, age, cause, gcm, ssp)]

europe_decomp <- decomp_period[,
  .(delta_e65 = sum(delta_e65),
    delta_sd = sum(delta_sd)),
  by = .(period, age, cause, gcm, ssp)]

europe_decomp[, country := "Europe"]

#---------------------------
# Save
#---------------------------

fwrite(country_analysis, file.path(path_out, "country_analysis.csv"))
fwrite(europe_analysis, file.path(path_out, "europe_analysis.csv"))
fwrite(country_lt, file.path(path_out, "country_lifespan_inequality.csv"))
fwrite(europe_lt, file.path(path_out, "europe_lifespan_inequality.csv"))
fwrite(country_decomp, file.path(path_out, "country_decomposition.csv"))
fwrite(europe_decomp, file.path(path_out, "europe_decomposition.csv"))

cat(sprintf("Country analysis: %s rows\n", format(nrow(country_analysis), big.mark = ",")))
cat(sprintf("Europe analysis: %s rows\n", format(nrow(europe_analysis), big.mark = ",")))
cat(sprintf("Country life tables: %s rows\n", format(nrow(country_lt), big.mark = ",")))
cat(sprintf("Europe life tables: %s rows\n", format(nrow(europe_lt), big.mark = ",")))