#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality / life expectancy pipeline -- Madrid pilot
#
# R Pipeline Step 01: Grouped attributable numbers under SSP3-7.0
#   Follows the exposure-response and ISIMIP3 bias-correction logic of
#   Masselot & Gasparrini (2025, R Code Part 3), narrowed to one city (Madrid),
#   one GCM and central (point-estimate) coefficients only -- no simulation
#   draws, so no empirical confidence interval is produced at this stage.
#   Two branches are contrasted: `with_cc` (full climate-change signal, ISIMIP3
#   bias-corrected) and `without_cc` (observed 2000-2019 daily sequences
#   repeated without an additional warming trend).
#
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(arrow)
  library(dlnm)
  library(splines)
  library(ggplot2)
})

#----- Global parameters and paths (inlined from Masselot & Gasparrini R Code Part 1;
# only the subset this point-estimate, single-GCM script actually uses)

# DLNM basis function: Masselot uses bs (B-spline, degree 2); our original used ns
varfun <- "bs"
vardegree <- 2

# Internal knots for the natural cubic spline cross-basis
knots_percentiles <- c(10, 75, 90)

path_tmean <- "data/tmeanproj.gz.parquet"

message("\n[01] Building Madrid SSP3-7.0 grouped attributable numbers...")

out_dir <- "results/phase1_madrid"
check_dir <- "results/checks"
fig_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

grouped_file <- file.path(out_dir, "01_attribution_grouped.csv")
checks_file <- file.path(check_dir, "01_attribution_checks.csv")
failures_file <- file.path(check_dir, "01_attribution_failures.csv")
fig_file <- file.path(fig_dir, "01_attribution_diagnostic.png")

city_id <- "ES001C"
city_name <- "Madrid"
gcm_name <- "GFDL_ESM4"
# Masselot (2025) 01_pkg_params.R: these two GCMs are excluded, leaving 19 of
# the 21 in tmeanproj.gz.parquet
gcm_excluded <- c("CMCC_CM2_SR5", "TaiESM1")
if (gcm_name %in% gcm_excluded) stop(sprintf("GCM %s is excluded in Masselot (2025).", gcm_name), call. = FALSE)
ssp_name <- "3"
variant_levels <- c("with_cc", "without_cc")
range_levels <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
branch_labels <- c(
  with_cc = "Projected climate change",
  without_cc = "No additional warming (Masselot recalibration to 2010-2014)"
)
range_labels <- c(
  ExtrCold = "Extreme cold",
  ModCold = "Moderate cold",
  ModHeat = "Moderate heat",
  ExtrHeat = "Extreme heat"
)
range_colors <- c(
  ExtrCold = "#2166ac",
  ModCold = "#67a9cf",
  ModHeat = "#ef8a62",
  ExtrHeat = "#b2182b"
)
hist_years_bias <- 2000:2014
# Calibration periods for the projections (Masselot 2025 `projrange`)
calib_breaks <- c(2015, seq(2030, 2100, by = 10))
# Without-climate-change counterfactual. Default follows Masselot (2025) and
# Simon's methods draft 2.4.3: each 5-year block of the calibrated GCM series is
# re-mapped (ISIMIP3) onto the calibrated 2010-2014 distribution, preserving the
# GCM's day-to-day weather but removing the warming trend.
# "era5_cycle" (observed 2000-2019 repeated forward) is kept as a sensitivity option.
counterfactual <- "masselot_demo"   # or "era5_cycle"
counterfactual_ref_year5 <- 2010L
hist_years_counterfactual <- 2000:2019
future_years <- 2020:2099
# Masselot (2025) removes 29 February from all daily series and uses 365-day years.
annualization_rule_text <- "sum(daily AN) / 365 (29 Feb removed)"

#----- Load prepared thresholds/observations and verify Madrid identity

load("data/prep_data.RData")
setDT(thresholds)
setDT(obs_data)

city_meta <- fread("data/city_results.csv")
city_meta <- unique(city_meta[URAU_CODE == city_id & LABEL == city_name & agegroup %in% c("65-74", "75-84", "85+")])
if (nrow(city_meta) != 3L) {
  stop(sprintf("Madrid identifier verification failed for %s.", city_id), call. = FALSE)
}
city_meta <- city_meta[match(c("65-74", "75-84", "85+"), agegroup)]

#----- Load the demographic domain built in Step 00

demography <- fread(file.path(out_dir, "00_demography_grouped.csv"))
demography <- demography[geo_id == city_id & ssp == as.integer(ssp_name)]
if (!nrow(demography)) {
  stop("Demographic domain for Madrid is empty; run 00_demography.R first.", call. = FALSE)
}

city_thresholds <- thresholds[URAU_CODE == city_id & agegroup %in% c("65-74", "75-84", "85+")]
if (nrow(city_thresholds) != 3L) stop(sprintf("Missing Madrid threshold rows in prep_data.RData for %s.", city_id), call. = FALSE)
if (anyNA(city_thresholds)) stop("Madrid thresholds contain NA values.", call. = FALSE)

#----- Load Madrid's historical (ERA5) observed temperature series

obs_city <- obs_data[URAU_CODE == city_id]
if (!nrow(obs_city)) stop("Historical observed temperatures for Madrid are missing.", call. = FALSE)
obs_city[, `:=`(
  year = as.integer(format(date, "%Y")),
  month = as.integer(format(date, "%m")),
  month_day = format(date, "%m-%d")
)]
#----- Load and bias-correct projected (GCM) temperature

proj_ds <- open_dataset(path_tmean)
tmean_all <- proj_ds %>%
  filter(URAU_CODE == city_id, ssp %in% c("hist", ssp_name)) %>%
  collect() %>%
  as.data.table()

if (!(paste0("tas_", gcm_name) %in% names(tmean_all))) {
  stop(sprintf("Projected temperature column tas_%s not found.", gcm_name), call. = FALSE)
}

tmean_all <- tmean_all[, .(date, ssp, tmean = get(paste0("tas_", gcm_name)))]
tmean_all[, `:=`(
  year = as.integer(format(date, "%Y")),
  month = as.integer(format(date, "%m")),
  month_day = format(date, "%m-%d")
)]

tmean_all <- tmean_all[month_day != "02-29"]
# Masselot (2025) 03_attribution.R: IITM_ESM has no SSP3 values for 2099; they
# are filled with the 2098 series.
if (gcm_name == "IITM_ESM" && ssp_name == "3") {
  fill_2098 <- tmean_all[ssp == ssp_name & year == 2098L, .(month_day, fill = tmean)]
  tmean_all <- merge(tmean_all, fill_2098, by = "month_day", all.x = TRUE, sort = FALSE)
  tmean_all[ssp == ssp_name & year == 2099L, tmean := fill]
  tmean_all[, fill := NULL]
  setorder(tmean_all, date)
}
tmean_future <- tmean_all[ssp == ssp_name & year %in% future_years]
hist_sim <- tmean_all[ssp == "hist" & year %in% hist_years_bias]
if (!nrow(tmean_future)) stop("Projected future temperature data are empty after city/GCM/year filtering.", call. = FALSE)
if (!nrow(hist_sim)) stop("Historical model temperatures for bias correction are missing from the projection table.", call. = FALSE)

#----- ISIMIP3 trend-preserving quantile mapping (Masselot & Gasparrini, 2025)

isimip3 <- function(obshist, simhist, simfut, yearobshist, yearsimhist, yearsimfut, detrend = TRUE) {
  if (detrend) {
    obstrend <- lm(obshist ~ yearobshist, na.action = na.exclude) |> predict() |> scale(scale = FALSE)
    simhisttrend <- lm(simhist ~ yearsimhist, na.action = na.exclude) |> predict() |> scale(scale = FALSE)
    simfuttrend <- lm(simfut ~ yearsimfut, na.action = na.exclude) |> predict() |> scale(scale = FALSE)

    obshist <- obshist - obstrend
    simhist <- simhist - simhisttrend
    simfut <- simfut - simfuttrend
  }

  ecdfobs <- ecdf(obshist)(obshist)
  deltaadd <- quantile(simfut, ecdfobs, na.rm = TRUE) - quantile(simhist, ecdfobs, na.rm = TRUE)
  obsfut <- deltaadd + obshist
  simfutcdf <- pnorm(simfut, mean(simfut, na.rm = TRUE), sd(simfut, na.rm = TRUE))
  calsimfut <- qnorm(p = simfutcdf, mean = mean(obsfut, na.rm = TRUE), sd = sd(obsfut, na.rm = TRUE))

  if (detrend) calsimfut <- calsimfut + simfuttrend
  calsimfut
}

#----- Exposure-response basis (bs, shared across variants and age groups)

obs_hist <- obs_city[year %in% hist_years_bias & month_day != "02-29"]
obs_repeat <- obs_city[year %in% hist_years_counterfactual, .(
  source_year = year,
  month_day,
  tmean_hist = tmean_obs
)]
if (anyDuplicated(obs_repeat, by = c("source_year", "month_day"))) {
  stop("Historical observed temperatures contain duplicate calendar dates.", call. = FALSE)
}
if (any(!is.finite(obs_repeat$tmean_hist))) {
  stop("Historical observed temperatures contain non-finite values.", call. = FALSE)
}

# ERF basis: knots and boundaries must be those used to estimate the ERFs, i.e.
# percentiles of the city's FULL observed ERA5-Land series (1990-2019), as in
# Masselot & Gasparrini (2025) 03_attribution.R (`tper`). Using any other window
# changes the shape of the published exposure-response function.
predper <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))
tper <- quantile(obs_city$tmean_obs, predper / 100, na.rm = TRUE)
knots <- tper[paste0(knots_percentiles, ".0%")]
bound <- range(tper)

#----- Validation fixture: reproduce Masselot et al. (2023) published historical
# excess deaths (city_results.csv) for this city's 65+ age groups from ERA5,
# the published coefficients and MMT. Checks the ERF reconstruction end to end.
fixture_rows <- list()
coef_fix <- fread("data/coefs.csv")[URAU_CODE == city_id]
for (agegrp in c("65-74", "75-84", "85+")) {
  cm <- city_meta[agegroup == agegrp]
  bfix <- as.matrix(coef_fix[agegroup == agegrp, .(b1, b2, b3, b4, b5)])
  bx <- suppressWarnings(onebasis(obs_city$tmean_obs, fun = varfun, degree = vardegree, knots = knots, Bound = bound))
  bc <- scale(bx, center = onebasis(cm$mmt, fun = varfun, degree = vardegree, knots = knots, Bound = bound), scale = FALSE)
  rr_fix <- exp(as.numeric(bc %*% t(bfix)))
  an_fix <- (1 - 1 / rr_fix) * cm$death / 365.25
  n_years <- uniqueN(obs_city$year)
  fixture_rows[[agegrp]] <- data.table(
    agegroup = agegrp,
    heat = sum(an_fix[obs_city$tmean_obs > cm$mmt]) / n_years,
    cold = sum(an_fix[obs_city$tmean_obs <= cm$mmt]) / n_years,
    pub_heat = cm$excess_heat_est,
    pub_cold = cm$excess_cold_est
  )
}
fixture <- rbindlist(fixture_rows)
fixture[, max_rel_error := pmax(abs(heat / pub_heat - 1), abs(cold / pub_cold - 1))]
fixture_tolerance <- 1e-3

coef_dt <- fread("data/coefs.csv")[URAU_CODE == city_id & agegroup %in% c("65-74", "75-84", "85+")]

city_results <- list()

#----- Build the two temperature variants, then loop age groups within each

#----- Calibrated series (Masselot 2025): "full" = with climate change,
# "demo" = recalibrated without additional warming

cal_src <- rbind(
  tmean_all[ssp == "hist" & year %in% hist_years_bias],
  tmean_all[ssp == ssp_name & year >= min(calib_breaks)]
)
cal_src[, calperiod := cut(year, c(min(hist_years_bias), calib_breaks), right = FALSE)]
if (anyNA(cal_src$calperiod)) stop("Temperature years fall outside the calibration periods.", call. = FALSE)
cal_src[, full := {
  m <- .BY$month
  obs_m <- obs_hist[month == m]
  sim_h_m <- hist_sim[month == m]
  isimip3(
    obshist = obs_m$tmean_obs, simhist = sim_h_m$tmean, simfut = tmean,
    yearobshist = obs_m$year, yearsimhist = sim_h_m$year, yearsimfut = year
  )
}, by = .(month, calperiod)]
cal_src[, year5 := (year %/% 5L) * 5L]
cal_src[, demo := {
  m <- .BY$month
  ref <- cal_src[month == m & year5 == counterfactual_ref_year5]
  isimip3(
    obshist = ref$full, simhist = full, simfut = full,
    yearobshist = ref$year, yearsimhist = year, yearsimfut = year
  )
}, by = .(month, year5)]

# Diagnostic: warming signal by decade in each branch (vs the 2010-2014 reference)
ref_mean <- cal_src[year5 == counterfactual_ref_year5, mean(full)]
warming_tbl <- cal_src[year >= min(future_years), .(
  full_minus_ref = mean(full) - ref_mean,
  demo_minus_ref = mean(demo) - ref_mean
), by = .(decade = (year %/% 10L) * 10L)]
max_demo_drift <- max(abs(warming_tbl$demo_minus_ref))

for (variant in variant_levels) {
  t_work <- copy(tmean_future)
  t_work[, days_in_year := 365L]
  coverage <- t_work[, .(rows = .N, unique_dates = uniqueN(date), expected_days = unique(days_in_year)), by = year]
  if (any(coverage$rows != coverage$expected_days) || any(coverage$unique_dates != coverage$expected_days)) {
    bad <- coverage[rows != expected_days | unique_dates != expected_days][1L]
    stop(sprintf("Incomplete daily temperature coverage for %s: year=%d, rows=%d, unique_dates=%d, expected=%d.", variant, bad$year, bad$rows, bad$unique_dates, bad$expected_days), call. = FALSE)
  }
  if (variant == "with_cc") {
    # Masselot (2025): calibrated by month x calibration period against ERA5 2000-2014
    t_work <- merge(t_work, cal_src[, .(date, tmean_variant = full)], by = "date", all.x = TRUE, sort = FALSE)
    t_work[, temp_branch := "with_cc"]
  } else if (counterfactual == "masselot_demo") {
    t_work <- merge(t_work, cal_src[, .(date, tmean_variant = demo)], by = "date", all.x = TRUE, sort = FALSE)
    t_work[, temp_branch := "without_cc"]
  } else if (counterfactual == "era5_cycle") {
    t_work[, source_year := min(hist_years_counterfactual) + ((year - min(future_years)) %% length(hist_years_counterfactual))]
    t_work <- merge(t_work, obs_repeat, by = c("source_year", "month_day"), all.x = TRUE, sort = FALSE)
    t_work[, tmean_variant := tmean_hist]
    t_work[, temp_branch := "without_cc"]
    t_work[, c("source_year", "tmean_hist") := NULL]
  } else {
    stop(sprintf("Unknown counterfactual '%s'.", counterfactual), call. = FALSE)
  }

  if (anyNA(t_work$tmean_variant)) {
    stop(sprintf("Temperature variant %s contains NA values after construction.", variant), call. = FALSE)
  }

  #----- Per age group: classify temperature range, compute AF/AN, annualise
  for (agegrp in c("65-74", "75-84", "85+")) {
    age_row <- city_thresholds[agegroup == agegrp]
    if (nrow(age_row) != 1L) stop(sprintf("Missing Madrid threshold row for age group %s.", agegrp), call. = FALSE)

    p2_5 <- age_row[["p2_5"]]
    p97_5 <- age_row[["p97_5"]]
    mmt <- age_row[["mmt"]]
    death_annual <- demography[agegroup == agegrp, .(year, death)]
    setorder(death_annual, year)
    age_work <- merge(copy(t_work), death_annual, by = "year", all.x = TRUE, sort = FALSE)
    if (nrow(age_work) != nrow(t_work) || !identical(names(age_work)[names(age_work) == "death"], "death")) {
      stop(sprintf("Annual-death join changed the temperature domain for %s/%s.", city_id, agegrp), call. = FALSE)
    }
    if (anyNA(age_work$death)) stop(sprintf("Missing annual deaths for %s/%s.", city_id, agegrp), call. = FALSE)

    # Split at the MMT first (cold vs heat), then at the fixed percentiles
    # (Lloyd et al. 2024). If the MMT lies above p97.5, all heat is extreme;
    # if below p2.5, all cold is extreme (Simon's methods draft 2.4.2).
    tv <- age_work$tmean_variant
    range_idx <- fifelse(
      tv < mmt,
      fifelse(tv < p2_5, "ExtrCold", "ModCold"),
      fifelse(tv >= p97_5, "ExtrHeat", "ModHeat")
    )

    b_fut <- onebasis(age_work$tmean_variant, fun = varfun, degree = vardegree, knots = knots, Bound = bound)
    b_mmt <- onebasis(mmt, fun = varfun, degree = vardegree, knots = knots, Bound = bound)
    b_centered <- scale(b_fut, center = b_mmt, scale = FALSE)

    age_coefs <- as.matrix(coef_dt[agegroup == agegrp, .(b1, b2, b3, b4, b5)])
    if (nrow(age_coefs) != 1L) stop(sprintf("Missing central coefficients for %s/%s.", city_id, agegrp), call. = FALSE)

    log_rr <- as.numeric(b_centered %*% t(age_coefs))
    # Masselot (2025): RR not allowed below 1 (`rr <- pmax(exp(bcen %*% coef), 1)`),
    # so extrapolated tails cannot produce negative ANs. Simon's methods draft 2.4.1.
    rr <- pmax(exp(log_rr), 1)
    af <- 1 - 1 / rr
    an_daily <- af * age_work$death

    annual <- data.table(
      year = age_work$year,
      agegroup = agegrp,
      range = range_idx,
      an = an_daily,
      branch = variant,
      days_in_year = age_work$days_in_year
    )
    annual <- annual[, .(an = sum(an) / unique(days_in_year)), by = .(year, agegroup, range, branch)]
    annual[, `:=`(geo_id = city_id, label = city_name, ssp = as.integer(ssp_name), gcm = gcm_name)]
    city_results[[length(city_results) + 1L]] <- annual
  }
}

#----- Assemble the full (branch x year x agegroup x range) domain, zero-filled

grouped <- rbindlist(city_results, use.names = TRUE)
setorder(grouped, branch, year, agegroup, range)

if (any(!grouped$range %in% range_levels)) stop("Unexpected temperature range labels produced.", call. = FALSE)

full_domain <- CJ(branch = variant_levels, year = future_years, agegroup = c("65-74", "75-84", "85+"), range = range_levels)
grouped <- merge(full_domain, grouped, by = c("branch", "year", "agegroup", "range"), all.x = TRUE, sort = FALSE)
grouped[is.na(an), an := 0]
grouped[, `:=`(geo_id = city_id, label = city_name, ssp = as.integer(ssp_name), gcm = gcm_name)]

year_days_tbl <- data.table(year = sort(unique(grouped$year)))
year_days_tbl[, days_in_year := 365L]

grouped <- merge(grouped, year_days_tbl, by = "year", all.x = TRUE, sort = FALSE)
grouped[, annualization_rule := annualization_rule_text]

#----- Invariant checks (project convention: a failing check stops the run)

checks <- data.table(
  check_name = c(
    "city_identifier_unique",
    "projected_gcm_present",
    "annualization_rule_documented",
    "year_specific_day_counts",
    "domain_complete_with_cc",
    "domain_complete_without_cc",
    "finite_signed_an",
    "masselot2023_fixture",
    "counterfactual_no_warming"
  ),
  status = c(
    if (nrow(city_meta) == 3L) "PASS" else "FAIL",
    if (uniqueN(grouped$gcm) == 1L && unique(grouped$gcm) == gcm_name) "PASS" else "FAIL",
    if (all(grouped$annualization_rule == annualization_rule_text)) "PASS" else "FAIL",
    if (all(year_days_tbl$days_in_year == 365L)) "PASS" else "FAIL",
    if (nrow(unique(grouped[branch == "with_cc"], by = c("year", "agegroup", "range"))) == nrow(full_domain[branch == "with_cc"]) ) "PASS" else "FAIL",
    if (nrow(unique(grouped[branch == "without_cc"], by = c("year", "agegroup", "range"))) == nrow(full_domain[branch == "without_cc"]) ) "PASS" else "FAIL",
    if (!any(!is.finite(grouped$an)) && all(grouped$an >= 0)) "PASS" else "FAIL",
    if (all(fixture$max_rel_error <= fixture_tolerance)) "PASS" else "FAIL",
    if (counterfactual != "masselot_demo" || max_demo_drift <= 0.25) "PASS" else "FAIL"
  ),
  value = c(
    nrow(city_meta),
    unique(grouped$gcm),
    annualization_rule_text,
    paste(sort(unique(year_days_tbl$days_in_year)), collapse = ","),
    nrow(grouped[branch == "with_cc"]),
    nrow(grouped[branch == "without_cc"]),
    sprintf("min_an=%g; max_an=%g", min(grouped$an), max(grouped$an)),
    sprintf("max_rel_error=%.2e", max(fixture$max_rel_error)),
    sprintf("%s; max |decadal mean - 2010-14 mean| = %.3f C (with_cc 2090s: %+.2f C)", counterfactual, max_demo_drift, warming_tbl[decade == 2090, full_minus_ref])
  ),
  threshold = c(
    "unique Madrid city identifier",
    gcm_name,
    annualization_rule_text,
    "365-day years (29 Feb removed)",
    sprintf("%d rows", nrow(full_domain[branch == "with_cc"])),
    sprintf("%d rows", nrow(full_domain[branch == "without_cc"])),
    "finite and non-negative AN (RR clamped at 1)",
    sprintf("<= %g vs published heat/cold excess (Masselot 2023)", fixture_tolerance),
    "<= 0.25 C decadal drift in the without-CC series"
  )
)

failures <- data.table()
if (any(checks$status == "FAIL")) {
  failures <- rbindlist(lapply(which(checks$status == "FAIL"), function(i) {
    data.table(
      failing_check = checks$check_name[i],
      observed_value = checks$value[i],
      expected_bound = checks$threshold[i]
    )
  }), fill = TRUE)
}

#----- Persist checks before exposing primary outputs

fwrite(checks, checks_file)
if (nrow(failures)) {
  fwrite(failures, failures_file)
  stop(sprintf("01_attribution.R failed %d invariant(s); see %s", nrow(failures), failures_file), call. = FALSE)
} else {
  if (file.exists(failures_file)) file.remove(failures_file)
  invisible(file.create(failures_file))
}

fwrite(grouped, grouped_file)

plot_dt <- grouped[, .(an = sum(an)), by = .(year, branch, range)]
plot_dt[, range := factor(range, levels = range_levels)]
setorder(plot_dt, branch, year, range)
plot_dt[, lower := cumsum(an) - an, by = .(year, branch)]
plot_dt[, upper := cumsum(an), by = .(year, branch)]
axis_values <- c(plot_dt$lower, plot_dt$upper)
y_min <- min(axis_values, na.rm = TRUE)
y_max <- max(axis_values, na.rm = TRUE)
p <- ggplot(plot_dt, aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = range), colour = NA, alpha = 0.85) +
  geom_line(aes(y = upper, colour = range), linewidth = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  facet_wrap(~branch, labeller = as_labeller(branch_labels)) +
  scale_fill_manual(values = range_colors, labels = range_labels, name = NULL) +
  scale_color_manual(values = range_colors, guide = "none") +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(
    title = "Madrid SSP3-7.0 annual temperature-attributable deaths by temperature range",
    subtitle = sprintf("One GCM (%s); shaded areas accumulate from extreme cold to extreme heat", gcm_name),
    x = "Year",
    y = "Annual temperature-attributable deaths"
  ) +
  theme_minimal(base_size = 11)
ggsave(fig_file, p, width = 11, height = 6, dpi = 160)

message("Saved grouped AN to ", grouped_file)
message("Saved checks to ", checks_file)
message("Saved diagnostic figure to ", fig_file)
