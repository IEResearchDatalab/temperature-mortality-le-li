################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 04_CityMortalitySchedules.R
#   Translate the national mortality backbone into city-level all-cause
#   mortality schedules.
#
# Method:
#   Cities inherit national mortality improvement ratios applied to their
#   own baseline mortality (disaggregated from grouped data via PCLM).
#
#   city_mx(year, age) = city_baseline_mx(age) * ratio(year, age)
#   where ratio(year, age) = national_mx_proj(year, age) / national_baseline_mx(age)
#
#   This preserves city-specific baseline differences while applying
#   national mortality trends.
#
# Output: city × year × age all-cause mx and population (single ages)
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")

library(ungroup)

#---------------------------
# Load inputs
#---------------------------

city_res <- readRDS(file.path(path_out, "city_res.rds"))
backbone <- readRDS(file.path(path_out, "mortality_backbone.rds"))

# Filter to selected cities if in demo mode
if (!is.null(city_subset)) {
  city_res <- city_res[URAU_CODE %in% city_subset]
}

city_codes <- sort(unique(city_res$URAU_CODE))
countries <- sort(unique(city_res$CNTR_CODE))

cat(sprintf("Building city mortality schedules for %d cities in %d countries...\n",
    length(city_codes), length(countries)))

#---------------------------
# PCLM disaggregation helper
#---------------------------

disaggregate_to_single <- function(group_values, age_breaks, nlast) {
  M <- pclm(x = age_breaks, y = group_values, nlast = nlast,
    out.step = 1, control = list(lambda = NA, kr = 2, deg = 3,
      max.iter = 200, tol = 1e-8, opt.method = "BIC"))
  as.numeric(fitted(M))
}

#---------------------------
# National baseline mx from last historical Eurostat year
# (backbone stores projection years only; use the raw life-table data as anchor)
#---------------------------

lt_hist <- readRDS(file.path(path_out, "lt_eu_raw_filtered.rds"))

# Most recent historical year per country (last year with complete mx data)
nat_baseline_yr <- lt_hist[, .(baseline_year = max(year)), by = country]

national_baseline <- merge(lt_hist, nat_baseline_yr, by = "country")[
  year == baseline_year,
  .(country, age, mx_nat_base = mx)]

cat(sprintf("National baseline years by country:\n"))
print(nat_baseline_yr)
rm(lt_hist)

#---------------------------
# Build city schedules
#---------------------------

all_city_schedules <- list()

for (cc in city_codes) {

  cntr <- unique(city_res[URAU_CODE == cc, CNTR_CODE])
  if (length(cntr) != 1) {
    cat(sprintf("  Skipping %s: no single country code\n", cc))
    next
  }

  # City baseline grouped data
  city_base <- city_res[URAU_CODE == cc & agegroup %in% agelabs]
  base_pop <- setNames(city_base$agepop, city_base$agegroup)
  base_death <- setNames(city_base$death, city_base$agegroup)

  # Disaggregate to single ages
  pop_single <- disaggregate_to_single(base_pop, agebreaks, pclm_nlast)
  pop_single <- pmax(pop_single, 0)
  pop_single <- pop_single / sum(pop_single) * sum(base_pop)

  death_single <- disaggregate_to_single(base_death, agebreaks, pclm_nlast)
  death_single <- pmax(death_single, 0)
  death_single <- death_single / sum(death_single) * sum(base_death)

  # City baseline mx at single ages
  mx_city_base <- death_single / pop_single
  mx_city_base[is.nan(mx_city_base) | is.infinite(mx_city_base)] <- 1e-10

  # National improvement ratio for all projection years
  natr_backbone <- backbone[country == cntr & year >= proj_year_min & year <= proj_year_max]
  natr_backbone <- merge(natr_backbone, national_baseline[country == cntr],
    by = "age", all.x = TRUE)
  natr_backbone[, mx_ratio := mx_proj / mx_nat_base]
  natr_backbone[is.na(mx_ratio) | mx_ratio <= 0, mx_ratio := 1]
  natr_backbone[is.infinite(mx_ratio), mx_ratio := 1]

  # Build city schedule for each year
  years <- sort(unique(natr_backbone$year))

  for (yr in years) {
    yr_ratio <- natr_backbone[year == yr, .(age, mx_ratio)]
    yr_ratio <- yr_ratio[age %in% single_ages]

    schedule <- data.table(
      city_code = cc, country = cntr,
      year = yr, age = single_ages,
      population = round(pop_single),
      mx = mx_city_base * yr_ratio$mx_ratio)

    schedule[mx <= 0 | is.na(mx), mx := 1e-10]
    schedule[, deaths := mx * population]

    all_city_schedules[[length(all_city_schedules) + 1]] <- schedule
  }

  if (which(cc == city_codes) %% 50 == 0) {
    cat(sprintf("  %d / %d cities processed\n", which(cc == city_codes), length(city_codes)))
  }
}

city_mx_all <- rbindlist(all_city_schedules)
setorder(city_mx_all, city_code, year, age)

#---------------------------
# Save
#---------------------------

saveRDS(city_mx_all, file.path(path_out, "city_mortality_schedules.rds"))

cat(sprintf("City mortality schedules: %d rows\n", nrow(city_mx_all)))
cat(sprintf("Cities: %d | Ages: %d-%d | Years: %d-%d\n",
    uniqueN(city_mx_all$city_code), min(city_mx_all$age), max(city_mx_all$age),
    min(city_mx_all$year), max(city_mx_all$year)))