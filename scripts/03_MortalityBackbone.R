################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#   a multi-city, multi-GCM, multi-SSP health impact assessment
#
# Pipeline: 03_MortalityBackbone.R
#   Build future country × year × age non-temperature mortality schedules
#   using an explicit, documented method based on Eurostat life-table inputs.
#
# Method:
#   1. Load historical Eurostat death rates (demo_mlifetable) by country-age-year
#   2. For each country-age combo, fit log-linear trend: ln(mx) ~ year
#      over the available historical period to estimate the annual
#      mortality improvement rate
#   3. Choose the most recent historical year with data as baseline
#   4. Project forward: mx_future(y) = mx_baseline * exp(b * (y - baseline_year))
#      where b is the age-specific log-linear slope
#   5. Output a clean country × year × age mx schedule for all projection years
#
# Assumptions:
#   - Mortality improvement follows a constant log-linear trend per age
#   - Trends are estimated from the full historical period (1990-2019)
#   - Open age interval is 100+ (mx at age 100 is carried forward as constant)
#   - If historical data is insufficient, borrow trend from adjacent ages
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")

#---------------------------
# Load historical life tables
#---------------------------

lt_hist <- readRDS(file.path(path_out, "lt_eu_raw_filtered.rds"))
pop_proj <- readRDS(file.path(path_out, "pop_eu_raw_filtered.rds"))

countries <- sort(unique(lt_hist$country))

cat(sprintf("Building mortality backbone for %d countries...\n", length(countries)))
cat(sprintf("Historical years: %d-%d\n", min(lt_hist$year), max(lt_hist$year)))

#---------------------------
# Helper: fit log-linear trend and project
#---------------------------

project_mx_age <- function(mx_hist_df, proj_years, min_years = 5) {
  # mx_hist_df: data.table with columns year, mx
  # proj_years: integer vector of years to project to
  # returns: data.table with columns year, mx, mx_proj, method_label

  mx_hist_df <- mx_hist_df[!is.na(mx) & mx > 0]
  n_hist <- nrow(mx_hist_df)

  if (n_hist >= min_years) {
    # Log-linear trend
    lm_fit <- lm(log(mx) ~ year, data = mx_hist_df)
    b <- coef(lm_fit)[["year"]]
  } else {
    b <- 0
  }

  baseline_year <- max(mx_hist_df$year)
  baseline_mx <- mx_hist_df[year == baseline_year, mx]
  if (length(baseline_mx) == 0) baseline_mx <- mx_hist_df[which.max(year), mx]

  proj_dt <- data.table(year = proj_years)
  proj_dt[, mx_proj := baseline_mx * exp(b * (year - baseline_year))]
  proj_dt[, mx_proj := pmax(mx_proj, 1e-10)]
  proj_dt[, method := "log_linear_trend"]
  proj_dt[, baseline_year := baseline_year]
  proj_dt[, b := b]

  proj_dt
}

#---------------------------
# Build backbone per country
#---------------------------

proj_years <- proj_year_min:proj_year_max
all_backbone <- list()

for (cntr in countries) {
  cntr_lt <- lt_hist[country == cntr]
  ages_in_data <- sort(unique(cntr_lt$age))
  min_hist_age <- min(ages_in_data)
  max_hist_age <- max(ages_in_data)

  for (ag in min_hist_age:max_hist_age) {
    age_data <- cntr_lt[age == ag]
    proj_dt <- project_mx_age(
      age_data[, .(year, mx)],
      proj_years
    )
    proj_dt[, country := cntr]
    proj_dt[, age := ag]
    all_backbone[[length(all_backbone) + 1]] <- proj_dt
  }
}

backbone_dt <- rbindlist(all_backbone)
setcolorder(backbone_dt, c("country", "year", "age", "mx_proj", "method", "baseline_year", "b"))

#---------------------------
# Handle missing ages by borrowing from adjacent ages
#---------------------------

for (cntr in countries) {
  cntr_ages <- sort(unique(backbone_dt[country == cntr]$age))
  all_ages <- single_ages
  missing_ages <- setdiff(all_ages, cntr_ages)

  for (ag in missing_ages) {
    nearest <- all_ages[which.min(abs(all_ages - ag))]
    if (nearest %in% cntr_ages) {
      template <- backbone_dt[country == cntr & age == nearest]
      template[, age := ag]
      backbone_dt <- rbind(backbone_dt, template)
    }
  }
}

setorder(backbone_dt, country, year, age)

#---------------------------
# Ensure all countries have complete age coverage
#---------------------------

full_grid <- CJ(country = countries, year = proj_years, age = single_ages)
backbone_dt <- merge(full_grid, backbone_dt, by = c("country", "year", "age"), all.x = TRUE)

# Forward-fill mx_proj within each country-age for any missing years
backbone_dt[, mx_proj := nafill(mx_proj, type = "locf"), by = .(country, age)]
backbone_dt[, mx_proj := nafill(mx_proj, type = "nocb"), by = .(country, age)]
backbone_dt[is.na(mx_proj) | mx_proj <= 0, mx_proj := 1e-10]

#---------------------------
# Save backbone
#---------------------------

saveRDS(backbone_dt, file.path(path_out, "mortality_backbone.rds"))

cat(sprintf("Mortality backbone: %d rows\n", nrow(backbone_dt)))
cat(sprintf("Countries: %d | Ages: %d-%d | Years: %d-%d\n",
    uniqueN(backbone_dt$country), min(backbone_dt$age), max(backbone_dt$age),
    min(backbone_dt$year), max(backbone_dt$year)))