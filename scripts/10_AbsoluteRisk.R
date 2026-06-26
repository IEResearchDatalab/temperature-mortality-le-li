################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 10_AbsoluteRisk.R
#   Absolute risk / attributable burden analysis (Branch A).
#   Implements the Lloyd et al. (2024) fixed vs rising longevity logic.
#
# Logic (per Lloyd reference):
#   AN_fixed(y,a,k)  = AF(y,a,k) × D_fixed(y,a)
#     where D_fixed = m_x(y0,a) × P(y,a)   [baseline mortality frozen]
#   AN_rising(y,a,k) = AF(y,a,k) × D_rising(y,a)
#     where D_rising = m_x(y,a)  × P(y,a)  [projected mortality]
#
#   Since 02_ComputeAN.R already computed AN_fixed using city_results.csv
#   deaths (frozen at ~2020), the rising-longevity AN is obtained by rescaling:
#
#     AN_rising = AN_fixed × (m_x(y,a) / m_x(y0,a))
#
#   Absolute-risk change (using standardised baseline population):
#     ΔAR_fixed(k)  = (AN_fixed(y1,k) − AN_fixed(y0,k)) / AN_fixed(y0,k) × 100
#     ΔAR_rising(k) = (AN_rising(y1,k) − AN_fixed(y0,k)) / AN_fixed(y0,k) × 100
#
#   Contribution of rising longevity:
#     Δ_longevity = ΔAR_fixed − ΔAR_rising
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")

#---------------------------
# Load inputs
#---------------------------

# Grouped ANs (climate-only, fixed baseline) from 02_ComputeAN.R
annual_ans <- fread(file.path(path_out, "ans_annual_all_cities.csv"))

# City mortality schedules (contain projected mx by age)
city_mx <- readRDS(file.path(path_out, "city_mortality_schedules.rds"))

# City metadata for baseline deaths/population
city_res <- readRDS(file.path(path_out, "city_res.rds"))

cat("Computing fixed vs rising longevity AN comparison...\n")

#---------------------------
# Step 1: AF = AN_fixed / baseline_deaths per age group
#---------------------------

# Baseline deaths from city_results.csv (frozen ~2020 level used in 02)
city_baseline <- city_res[agegroup %in% agelabs,
  .(city_code = URAU_CODE, age_group = agegroup, baseline_deaths = death)]

ans_with_af <- merge(annual_ans, city_baseline,
  by = c("city_code", "age_group"), all.x = TRUE)

ans_with_af[baseline_deaths > 0, AF := an_est / baseline_deaths]
ans_with_af[baseline_deaths == 0 | is.na(baseline_deaths), AF := 0]

#---------------------------
# Step 2: projected deaths per age group from city schedules
#---------------------------

# Map single ages to Masselot age groups
city_mx[, age_group := fcase(
  age >= 20 & age <= 44, "20-44",
  age >= 45 & age <= 64, "45-64",
  age >= 65 & age <= 74, "65-74",
  age >= 75 & age <= 84, "75-84",
  age >= 85,             "85+"
)]

proj_deaths <- city_mx[age >= 20,
  .(deaths_proj = sum(deaths)),
  by = .(city_code, year, age_group)]

#---------------------------
# Step 3: compute AN_rising
#---------------------------

ans_full <- merge(ans_with_af, proj_deaths,
  by = c("city_code", "year", "age_group"), all.x = TRUE)

# Rising: same AF applied to projected deaths
ans_full[, an_rising := AF * deaths_proj]
ans_full[is.na(an_rising), an_rising := 0]

#---------------------------
# Step 4: baseline year absolute risk (use first projection year)
#---------------------------

year0 <- proj_year_min  # 2020

ans_y0 <- ans_full[year == year0,
  .(an_fixed_y0 = sum(an_est), an_rising_y0 = sum(an_rising)),
  by = .(city_code, gcm, ssp, temp_range)]

ans_summary <- ans_full[,
  .(an_fixed = sum(an_est), an_rising = sum(an_rising)),
  by = .(city_code, gcm, ssp, year, temp_range)]

ans_summary <- merge(ans_summary, ans_y0,
  by = c("city_code", "gcm", "ssp", "temp_range"), all.x = TRUE)

# % change relative to fixed baseline at year0
ans_summary[an_fixed_y0 > 0,
  delta_ar_fixed  := (an_fixed  - an_fixed_y0) / an_fixed_y0 * 100]
ans_summary[an_fixed_y0 > 0,
  delta_ar_rising := (an_rising - an_fixed_y0) / an_fixed_y0 * 100]

# Contribution of rising longevity in percentage points
ans_summary[, contribution_longevity_pp := delta_ar_fixed - delta_ar_rising]

#---------------------------
# Save
#---------------------------

fwrite(ans_summary,
  file.path(path_out, "absolute_risk_summary.csv"))

fwrite(ans_full[, .(city_code, gcm, ssp, year, age_group, temp_range,
  an_fixed = an_est, an_rising, AF)],
  file.path(path_out, "absolute_risk_detail.csv"))

cat(sprintf("Absolute risk summary: %s rows\n",
  format(nrow(ans_summary), big.mark = ",")))

# Quick sanity print for Madrid
madrid_ar <- ans_summary[city_code == "ES001C" & ssp == "2" &
  gcm == unique(gcm)[1] & year %in% c(2030, 2050, 2090)]
if (nrow(madrid_ar) > 0) {
  cat("\nMadrid sample (SSP3-7.0):\n")
  print(madrid_ar[, .(year, temp_range,
    an_fixed = round(an_fixed, 1), an_rising = round(an_rising, 1),
    delta_ar_fixed = round(delta_ar_fixed, 1),
    delta_ar_rising = round(delta_ar_rising, 1),
    contribution_longevity_pp = round(contribution_longevity_pp, 1))])
}