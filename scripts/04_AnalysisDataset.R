################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 04_AnalysisDataset.R
#   Combine ANs with population and all-cause mortality
#   Build analysis dataset for all cities
#
################################################################################

if (length(ls()) == 0) source("03_Disaggregate.R")

library(ungroup)

#---------------------------
# Load single-age ANs
#---------------------------

single_an <- fread(file.path(path_out, "ans_single_age_all_cities.csv"))

#---------------------------
# Build analysis dataset per city
#---------------------------

all_datasets <- list()

for (cc in unique(single_an$city_code)) {

  cat(sprintf("Building dataset: %s\n", cc))

  city_res <- city_res[URAU_CODE == cc]
  city_base <- city_res[agegroup %in% agelabs]

  base_pop <- setNames(city_base$agepop, city_base$agegroup)
  base_death <- setNames(city_base$death, city_base$agegroup)

  # PCLM disaggregation of population
  pop_pclm <- pclm(x = agebreaks, y = base_pop,
    nlast = pclm_nlast, out.step = 1,
    control = list(lambda = NA, kr = 2, deg = 3))
  pop_single <- pmax(as.numeric(fitted(pop_pclm)), 0)
  pop_single <- pop_single / sum(pop_single) * sum(base_pop)

  # PCLM disaggregation of deaths
  death_pclm <- pclm(x = agebreaks, y = base_death,
    nlast = pclm_nlast, out.step = 1,
    control = list(lambda = NA, kr = 2, deg = 3))
  death_single <- pmax(as.numeric(fitted(death_pclm)), 0)
  death_single <- death_single / sum(death_single) * sum(base_death)

  mx_single_base <- death_single / pop_single
  mx_single_base[is.nan(mx_single_base) | is.infinite(mx_single_base)] <- 0

  pop_dt <- data.table(age = single_ages, population = round(pop_single),
    mx_baseline = mx_single_base)

  # For each SSP, GCM, year: build dataset
  city_an <- single_an[city_code == cc]

  # Build projected mortality by merging Eurostat trend with city baseline
  cntr <- unique(city_res[URAU_CODE == cc, CNTR_CODE])[1]
  cntr_mort <- all_mort[age %in% single_ages]  # from 01_PrepData
  mx_ratio_dt <- cntr_mort[, .(age, year = year, mx_ratio)]

  for (gc in unique(city_an$gcm)) {
    for (ss in unique(city_an$ssp)) {
      gc_an <- city_an[gcm == gc & ssp == ss]

      years <- unique(gc_an$year)

      for (yr in years) {
        yr_an <- gc_an[year == yr]

        # Projected mx = city baseline * national trend ratio
        yr_mx <- merge(data.table(age = single_ages, mx_base = mx_single_base),
          mx_ratio_dt[year == yr, .(age, mx_ratio)], by = "age", all.x = TRUE)
        yr_mx[is.na(mx_ratio), mx_ratio := 1]
        yr_mx[, mx := mx_base * mx_ratio]
        an_wide <- dcast(yr_an, year + age ~ temp_range,
          value.var = "AN", fill = 0)

        dataset <- data.table(
          city_code = cc, gcm = gc, ssp = ss,
          year = yr, age = single_ages,
          population = pop_dt$population,
          mx = yr_mx$mx,
          deaths = yr_mx$mx * pop_dt$population,
          an_extreme_cold = an_wide$extreme_cold[match(single_ages, an_wide$age)],
          an_moderate_cold = an_wide$moderate_cold[match(single_ages, an_wide$age)],
          an_moderate_heat = an_wide$moderate_heat[match(single_ages, an_wide$age)],
          an_extreme_heat = an_wide$extreme_heat[match(single_ages, an_wide$age)])

        dataset[is.na(an_extreme_cold), an_extreme_cold := 0]
        dataset[is.na(an_moderate_cold), an_moderate_cold := 0]
        dataset[is.na(an_moderate_heat), an_moderate_heat := 0]
        dataset[is.na(an_extreme_heat), an_extreme_heat := 0]

        dataset[, an_total := an_extreme_cold + an_moderate_cold +
          an_moderate_heat + an_extreme_heat]
        dataset[, deaths_baseline := deaths - an_total]
        dataset[deaths_baseline < 0, deaths_baseline := 0]

        all_datasets[[length(all_datasets) + 1]] <- dataset
      }
    }
  }
}

analysis_all <- rbindlist(all_datasets)
fwrite(analysis_all, file.path(path_out, "analysis_dataset_all_cities.csv"))
cat(sprintf("Analysis dataset: %s rows\n",
  format(nrow(analysis_all), big.mark = ",")))