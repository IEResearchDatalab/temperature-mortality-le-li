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

source("R_pipeline/01_initialize.R")

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
ssp_name <- "3"
variant_levels <- c("with_cc", "without_cc")
range_levels <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
branch_labels <- c(
  with_cc = "Projected climate change",
  without_cc = "No additional warming"
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
hist_years_counterfactual <- 2000:2019
future_years <- 2020:2099
annualization_rule_text <- "sum(daily AN) / actual days in year"

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

obs_hist <- obs_city[year %in% hist_years_bias]
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

obs_temp_vals <- obs_hist$tmean_obs
knots <- quantile(obs_temp_vals, knots_percentiles / 100, na.rm = TRUE)
bound <- range(obs_temp_vals, na.rm = TRUE)

coef_dt <- fread("data/coefs.csv")[URAU_CODE == city_id & agegroup %in% c("65-74", "75-84", "85+")]

city_results <- list()

#----- Build the two temperature variants, then loop age groups within each

for (variant in variant_levels) {
  t_work <- copy(tmean_future)
  t_work[, days_in_year := as.integer(as.Date(sprintf("%d-12-31", year)) - as.Date(sprintf("%d-01-01", year)) + 1L)]
  coverage <- t_work[, .(rows = .N, unique_dates = uniqueN(date), expected_days = unique(days_in_year)), by = year]
  if (any(coverage$rows != coverage$expected_days) || any(coverage$unique_dates != coverage$expected_days)) {
    bad <- coverage[rows != expected_days | unique_dates != expected_days][1L]
    stop(sprintf("Incomplete daily temperature coverage for %s: year=%d, rows=%d, unique_dates=%d, expected=%d.", variant, bad$year, bad$rows, bad$unique_dates, bad$expected_days), call. = FALSE)
  }
  if (variant == "with_cc") {
    t_work[, tmean_variant := {
      m <- .BY$month
      obs_m <- obs_hist[month == m]
      sim_h_m <- hist_sim[month == m]
      if (nrow(obs_m) < 10L || nrow(sim_h_m) < 10L) return(as.numeric(NA))
      isimip3(
        obshist = obs_m$tmean_obs,
        simhist = sim_h_m$tmean,
        simfut = tmean,
        yearobshist = obs_m$year,
        yearsimhist = sim_h_m$year,
        yearsimfut = year
      )
    }, by = month]
    t_work[, temp_branch := "with_cc"]
  } else {
    t_work[, source_year := min(hist_years_counterfactual) + ((year - min(future_years)) %% length(hist_years_counterfactual))]
    t_work <- merge(t_work, obs_repeat, by = c("source_year", "month_day"), all.x = TRUE, sort = FALSE)
    t_work[, tmean_variant := tmean_hist]
    t_work[, temp_branch := "without_cc"]
    t_work[, c("source_year", "tmean_hist") := NULL]
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

    range_idx <- fcase(
      age_work$tmean_variant < p2_5, "ExtrCold",
      age_work$tmean_variant < mmt, "ModCold",
      age_work$tmean_variant < p97_5, "ModHeat",
      default = "ExtrHeat"
    )

    b_fut <- onebasis(age_work$tmean_variant, fun = varfun, degree = vardegree, knots = knots, Bound = bound)
    b_mmt <- onebasis(mmt, fun = varfun, degree = vardegree, knots = knots, Bound = bound)
    b_centered <- scale(b_fut, center = b_mmt, scale = FALSE)

    age_coefs <- as.matrix(coef_dt[agegroup == agegrp, .(b1, b2, b3, b4, b5)])
    if (nrow(age_coefs) != 1L) stop(sprintf("Missing central coefficients for %s/%s.", city_id, agegrp), call. = FALSE)

    log_rr <- as.numeric(b_centered %*% t(age_coefs))
    af <- 1 - exp(-log_rr)
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
year_days_tbl[, days_in_year := as.integer(as.Date(sprintf("%d-12-31", year)) - as.Date(sprintf("%d-01-01", year)) + 1L)]

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
    "finite_signed_an"
  ),
  status = c(
    if (nrow(city_meta) == 3L) "PASS" else "FAIL",
    if (uniqueN(grouped$gcm) == 1L && unique(grouped$gcm) == gcm_name) "PASS" else "FAIL",
    if (all(grouped$annualization_rule == annualization_rule_text)) "PASS" else "FAIL",
    if (length(unique(year_days_tbl$days_in_year)) > 1L && all(year_days_tbl$days_in_year %in% c(365L, 366L))) "PASS" else "FAIL",
    if (nrow(unique(grouped[branch == "with_cc"], by = c("year", "agegroup", "range"))) == nrow(full_domain[branch == "with_cc"]) ) "PASS" else "FAIL",
    if (nrow(unique(grouped[branch == "without_cc"], by = c("year", "agegroup", "range"))) == nrow(full_domain[branch == "without_cc"]) ) "PASS" else "FAIL",
    if (!any(!is.finite(grouped$an))) "PASS" else "FAIL"
  ),
  value = c(
    nrow(city_meta),
    unique(grouped$gcm),
    annualization_rule_text,
    paste(sort(unique(year_days_tbl$days_in_year)), collapse = ","),
    nrow(grouped[branch == "with_cc"]),
    nrow(grouped[branch == "without_cc"]),
    sprintf("min_an=%g; max_an=%g", min(grouped$an), max(grouped$an))
  ),
  threshold = c(
    "unique Madrid city identifier",
    gcm_name,
    annualization_rule_text,
    "actual calendar days by year",
    sprintf("%d rows", nrow(full_domain[branch == "with_cc"])),
    sprintf("%d rows", nrow(full_domain[branch == "without_cc"])),
    "finite signed AN"
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
y_min <- min(plot_dt$an, na.rm = TRUE)
y_max <- max(plot_dt$an, na.rm = TRUE)
p <- ggplot(plot_dt, aes(x = year, y = an, color = range)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~branch, labeller = as_labeller(branch_labels)) +
  scale_color_manual(values = range_colors, labels = range_labels, name = NULL) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(
    title = "Madrid SSP3-7.0 annual temperature-attributable deaths",
    subtitle = sprintf("One GCM (%s); central coefficients; annualized with actual calendar days", gcm_name),
    x = "Year",
    y = "Annual temperature-attributable deaths"
  ) +
  theme_minimal(base_size = 11)
ggsave(fig_file, p, width = 11, height = 6, dpi = 160)

message("Saved grouped AN to ", grouped_file)
message("Saved checks to ", checks_file)
message("Saved diagnostic figure to ", fig_file)
