#!/usr/bin/env Rscript

################################################################################
# Diagnostic requested by Simon Lloyd (25 Sep 2026): inspect whether the
# time-block decomposition pattern is explained by temperature trajectories.
# For one city, reconstruct the exact Masselot with- and without-climate-change
# daily series for all 19 GCMs, then report five-year temperature mean,
# variance, and frequency of days above the fixed ERA5 1990-2019 p97.5.
################################################################################

source("pipeline/00_pkg_params.R")

path_tmean <- "data/tmeanproj.gz.parquet"
diag_dir <- Sys.getenv("DIAG_DIR", file.path("results/europe/diagnostics", city_id))
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(path_tmean)) stop("Missing ", path_tmean, call. = FALSE)

load("data/prep_data.RData")
setDT(obs_data); setDT(thresholds)
obs_city <- obs_data[URAU_CODE == city_id]
if (!nrow(obs_city)) stop("Missing observed temperatures for ", city_id, call. = FALSE)
obs_city[, `:=`(
  year = as.integer(format(date, "%Y")),
  month = as.integer(format(date, "%m")),
  month_day = format(date, "%m-%d")
)]
obs_hist <- obs_city[year %between% histrange & month_day != "02-29"]
p97_5 <- unique(thresholds[URAU_CODE == city_id & agegroup %in% agelabs, p97_5])
if (length(p97_5) != 1L) stop("City p97.5 threshold is not unique.", call. = FALSE)

raw <- open_dataset(path_tmean) %>%
  filter(URAU_CODE == !!city_id, ssp %in% c("hist", !!ssp_name)) %>%
  collect() %>%
  as.data.table()
needed_cols <- paste0("tas_", gcmlist)
if (!all(needed_cols %in% names(raw))) stop("Projection data lack one or more required GCM columns.", call. = FALSE)

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

one_gcm <- function(gcm) {
  x <- raw[, .(date, ssp, tmean = get(paste0("tas_", gcm)))]
  x[, `:=`(
    year = as.integer(format(date, "%Y")), month = as.integer(format(date, "%m")),
    month_day = format(date, "%m-%d")
  )]
  x <- x[month_day != "02-29"]
  if (gcm == "IITM_ESM" && ssp_name == "3") {
    fill <- x[ssp == ssp_name & year == 2098L, .(month_day, fill = tmean)]
    x <- merge(x, fill, by = "month_day", all.x = TRUE, sort = FALSE)
    x[ssp == ssp_name & year == 2099L, tmean := fill]
    x[, fill := NULL]
  }
  hist_sim <- x[ssp == "hist" & year %between% histrange]
  cal <- rbind(x[ssp == "hist" & year %between% histrange], x[ssp == ssp_name & year >= min(projrange)])
  cal[, calperiod := cut(year, c(histrange[1], projrange), right = FALSE)]
  cal[, full := {
    m <- .BY$month
    oh <- obs_hist[month == m]
    sh <- hist_sim[month == m]
    isimip3(oh$tmean_obs, sh$tmean, tmean, oh$year, sh$year, year)
  }, by = .(month, calperiod)]
  cal[, year5 := (year %/% 5L) * 5L]
  cal[, demo := {
    m <- .BY$month
    ref <- cal[month == m & year5 == counterfactual_ref_year5]
    isimip3(ref$full, full, full, ref$year, year, year)
  }, by = .(month, year5)]

  z <- cal[year %in% future_years, .(date, year, with_cc = full, without_cc = demo)]
  z <- melt(z, id.vars = c("date", "year"), variable.name = "branch", value.name = "temperature")
  z[, period_start := 2020L + ((year - 2020L) %/% 5L) * 5L]
  ans <- z[, .(
    mean_temperature = mean(temperature),
    temperature_variance = var(temperature),
    extreme_heat_days_per_year = sum(temperature >= p97_5) / 5
  ), by = .(branch, period_start)]
  ans[, gcm := gcm]
  ans
}

by_gcm <- rbindlist(lapply(gcmlist, one_gcm))
setcolorder(by_gcm, c("gcm", "branch", "period_start", "mean_temperature", "temperature_variance", "extreme_heat_days_per_year"))
setorder(by_gcm, gcm, branch, period_start)

expected_rows <- length(gcmlist) * length(branch_levels) * length(seq(2020L, 2095L, 5L))
if (nrow(by_gcm) != expected_rows || any(!is.finite(unlist(by_gcm[, -c("gcm", "branch")]))) ||
    any(by_gcm$extreme_heat_days_per_year < 0 | by_gcm$extreme_heat_days_per_year > 365)) {
  stop("Temperature diagnostic grid or values are invalid.", call. = FALSE)
}

ensemble <- by_gcm[, .(
  mean_temperature = mean(mean_temperature),
  temperature_variance = mean(temperature_variance),
  extreme_heat_days_per_year = mean(extreme_heat_days_per_year),
  mean_temperature_gcm_sd = sd(mean_temperature),
  extreme_heat_days_gcm_sd = sd(extreme_heat_days_per_year)
), by = .(branch, period_start)]

fwrite(by_gcm, file.path(diag_dir, "temperature_pattern_by_gcm.csv"))
fwrite(ensemble, file.path(diag_dir, "temperature_pattern_ensemble.csv"))
fwrite(data.table(
  check_name = c("gcm_period_grid", "values_finite", "extreme_day_range"),
  status = "PASS",
  value = c(sprintf("%d rows", nrow(by_gcm)), "all finite", sprintf("%.1f to %.1f days/year",
    min(by_gcm$extreme_heat_days_per_year), max(by_gcm$extreme_heat_days_per_year)))
), file.path(diag_dir, "temperature_pattern_checks.csv"))

plot_dt <- melt(ensemble,
  id.vars = c("branch", "period_start"),
  measure.vars = c("mean_temperature", "temperature_variance", "extreme_heat_days_per_year"),
  variable.name = "measure", value.name = "value"
)
plot_dt[, measure := factor(measure,
  levels = c("mean_temperature", "temperature_variance", "extreme_heat_days_per_year"),
  labels = c("Mean temperature (deg C)", "Temperature variance (deg C squared)", "Extreme-heat days per year"))]
p <- ggplot(plot_dt, aes(period_start + 2, value, colour = branch)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  facet_wrap(~measure, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = branch_colors, labels = branch_labels, name = NULL) +
  labs(
    title = sprintf("%s: temperature-pattern diagnostic", city_name),
    subtitle = sprintf("%s, mean across 19 GCMs; extreme threshold = %.2f deg C", ssplabs[ssp_name], p97_5),
    x = "5-year period (midpoint)", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(diag_dir, "temperature_pattern_ensemble.png"), p, width = 10, height = 9, dpi = 160)

print(ensemble, digits = 4)
message("Saved temperature-pattern diagnostic to ", diag_dir)
