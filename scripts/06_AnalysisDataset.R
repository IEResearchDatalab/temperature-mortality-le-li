################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 06_AnalysisDataset.R (rewritten)
#   Build the analysis dataset by combining projected all-cause mortality
#   with temperature-attributable deaths.
#
# Key design: the climate-only ANs from 02_ComputeAN.R are computed under
# fixed baseline deaths (~2020). To make them commensurate with projected
# mortality, we:
#   1. Compute AF = AN_fixed / baseline_deaths (per single age)
#   2. Apply AF to projected deaths → AN_projected
#   3. Residual = projected_deaths − AN_projected
#
# This enforces: all-cause = residual + extreme_cold + moderate_cold
#                             + moderate_heat + extreme_heat
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")
library(ungroup)

#---------------------------
# Load inputs
#---------------------------

city_schedules <- readRDS(file.path(path_out, "city_mortality_schedules.rds"))
city_res <- readRDS(file.path(path_out, "city_res.rds"))
single_an <- fread(file.path(path_out, "ans_single_age_all_cities.csv"))

common_cities <- intersect(
  unique(city_schedules$city_code),
  unique(single_an$city_code)
)
city_schedules <- city_schedules[city_code %in% common_cities]
city_res <- city_res[URAU_CODE %in% common_cities]
single_an <- single_an[city_code %in% common_cities]

cat(sprintf("Cities in analysis: %d\n", length(common_cities)))

#---------------------------
# Step 1 — AF from fixed-baseline ANs and city baseline deaths
#---------------------------

cat("Computing AF from fixed-baseline ANs and city baseline deaths...\n")

# Baseline deaths at single ages via PCLM
disaggregate_to_single <- function(grp_vals, brks, nlast) {
  M <- pclm(x = brks, y = grp_vals, nlast = nlast,
    out.step = 1, control = list(lambda = NA, kr = 2, deg = 3,
      max.iter = 200, tol = 1e-8, opt.method = "BIC"))
  as.numeric(fitted(M))
}

baseline_deaths_single <- rbindlist(lapply(common_cities, function(cc) {
  cd <- city_res[URAU_CODE == cc & agegroup %in% agelabs]
  grp_d <- setNames(cd$death, cd$agegroup)
  ds <- disaggregate_to_single(grp_d, agebreaks, pclm_nlast)
  ds <- pmax(ds, 0)
  ds <- ds / sum(ds) * sum(grp_d)
  data.table(city_code = cc, age = single_ages, baseline_deaths = ds)
}))

# Merge AF
single_an <- merge(single_an, baseline_deaths_single, by = c("city_code", "age"), all.x = TRUE)
single_an[baseline_deaths > 0, AF := AN / baseline_deaths]
single_an[baseline_deaths <= 0 | is.na(baseline_deaths), AF := 0]
single_an[is.na(AF), AF := 0]

#---------------------------
# Step 2 — Apply AF to projected deaths
#---------------------------

cat("Applying AF to projected deaths...\n")

# Merge projected deaths from city_schedules
proj_deaths <- city_schedules[, .(city_code, year, age, population, mx,
  deaths_proj = deaths)]
comb <- merge(single_an, proj_deaths, by = c("city_code", "year", "age"),
  all.x = TRUE, allow.cartesian = TRUE)

comb[, AN_proj := AF * deaths_proj]

# Cast projected ANs wide by temp_range
wide_proj <- dcast(comb, city_code + gcm + ssp + year + age + population + mx +
  deaths_proj ~ temp_range, value.var = "AN_proj", fill = 0)
for (tr in c("extreme_cold", "moderate_cold", "moderate_heat", "extreme_heat")) {
  if (!tr %in% names(wide_proj)) wide_proj[, (tr) := 0]
}

wide_proj[, an_total := extreme_cold + moderate_cold + moderate_heat + extreme_heat]
wide_proj[, residual_deaths := deaths_proj - an_total]

#---------------------------
# Step 3 — Check identity
#---------------------------

n_neg <- sum(wide_proj$residual_deaths < -0.001)
if (n_neg > 0) {
  worst <- wide_proj[which.min(residual_deaths)]
  cat(sprintf("FAIL: %d rows with residual < −0.001\n", n_neg))
  cat(sprintf("Worst: city=%s yr=%d age=%d gcm=%s ssp=%s\n",
    worst$city_code, worst$year, worst$age, worst$gcm, worst$ssp))
  cat(sprintf("  deaths=%.4f an_total=%.4f residual=%.4f\n",
    worst$deaths_proj, worst$an_total, worst$residual_deaths))
  stop("Bookkeeping identity failed.")
}

wide_proj[residual_deaths < 0, residual_deaths := 0]

#---------------------------
# Step 4 — Save
#---------------------------

setnames(wide_proj, "deaths_proj", "deaths")
setorder(wide_proj, city_code, gcm, ssp, year, age)

fwrite(wide_proj, file.path(path_out, "analysis_dataset_all_cities.csv"))
cat(sprintf("Analysis dataset: %s rows\nIdentity OK.\n",
  format(nrow(wide_proj), big.mark = ",")))