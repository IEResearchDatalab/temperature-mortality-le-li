#!/usr/bin/env Rscript

################################################################################
# Diagnostic requested by Simon Lloyd (25 Sep 2026): why does moderate heat
# contribute so little by the end of the century?
#
# This script checks both parts of the mechanism:
#   1. the shape of every city-age ERF between its MMT and fixed ERA5
#      1990-2019 p97.5 threshold; and
#   2. the observed split of projected heat-attributable deaths between the
#      moderate- and extreme-heat categories in Object 1.
################################################################################

source("pipeline/00_pkg_params.R")

obj1_file <- Sys.getenv("OBJECT1_FILE", "results/europe/collected/object1_dataset.parquet")
diag_dir <- Sys.getenv("DIAG_DIR", "results/europe/diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(obj1_file)) stop("Missing Object 1: ", obj1_file, call. = FALSE)
if (!file.exists("data/prep_data.RData")) stop("Missing data/prep_data.RData.", call. = FALSE)

load("data/prep_data.RData")
setDT(obs_data); setDT(thresholds)
coef_dt <- fread("data/coefs.csv")[agegroup %in% agelabs]
meta <- unique(fread("data/city_results.csv")[, .(URAU_CODE, region, country = CNTR_CODE)])

erf <- rbindlist(lapply(sort(unique(thresholds$URAU_CODE)), function(id) {
  obs <- obs_data[URAU_CODE == id, tmean_obs]
  tper <- quantile(obs, predper / 100, na.rm = TRUE)
  knots <- tper[paste0(varper, ".0%")]
  bound <- range(tper)
  p99 <- as.numeric(quantile(obs, 0.99, na.rm = TRUE))
  rbindlist(lapply(agelabs, function(a) {
    tr <- thresholds[URAU_CODE == id & agegroup == a]
    cf <- coef_dt[URAU_CODE == id & agegroup == a, .(b1, b2, b3, b4, b5)]
    if (nrow(tr) != 1L || nrow(cf) != 1L) stop("Incomplete ERF data for ", id, "/", a, call. = FALSE)
    mmt <- tr$mmt; p97 <- tr$p97_5
    mid <- if (mmt < p97) (mmt + p97) / 2 else NA_real_
    temps <- c(mmt, mid, p97, p99)
    b <- suppressWarnings(onebasis(temps, fun = varfun, degree = vardegree, knots = knots, Bound = bound))
    bc <- scale(b, center = b[1, ], scale = FALSE)
    rr <- pmax(exp(as.numeric(bc %*% t(as.matrix(cf)))), 1)
    data.table(
      URAU_CODE = id, agegroup = a, mmt = mmt, p97_5 = p97, p99 = p99,
      mmt_percentile = 100 * mean(obs <= mmt),
      moderate_band_exists = mmt < p97,
      moderate_band_width_c = pmax(p97 - mmt, 0),
      rr_mid_moderate = if (mmt < p97) rr[2] else NA_real_,
      rr_at_p97_5 = rr[3], rr_at_p99 = rr[4]
    )
  }))
}))
erf <- merge(erf, meta, by = "URAU_CODE", all.x = TRUE)
expected_erf_rows <- uniqueN(meta$URAU_CODE) * length(agelabs)
if (nrow(erf) != expected_erf_rows || uniqueN(erf$URAU_CODE) != uniqueN(meta$URAU_CODE)) {
  stop("ERF diagnostic city-age coverage is incomplete.", call. = FALSE)
}
if (any(!is.finite(erf$rr_at_p97_5)) || any(!is.finite(erf$rr_at_p99)) ||
    any(erf$rr_at_p97_5 < 1) || any(erf$rr_at_p99 < 1)) {
  stop("ERF diagnostic produced invalid RR values.", call. = FALSE)
}
fwrite(erf, file.path(diag_dir, "moderate_heat_erf_city_age.csv"))

summarise_erf <- function(x, label) {
  x <- copy(x)
  x[, geography := label]
  x[, .(
    n_city_age = .N,
    share_no_moderate_band = mean(!moderate_band_exists),
    median_mmt_percentile = median(mmt_percentile),
    median_moderate_band_width_c = median(moderate_band_width_c),
    median_rr_mid_moderate = median(rr_mid_moderate, na.rm = TRUE),
    median_rr_at_p97_5 = median(rr_at_p97_5),
    median_rr_at_p99 = median(rr_at_p99),
    share_rr_p97_below_1_01 = mean(rr_at_p97_5 < 1.01),
    share_rr_p97_below_1_05 = mean(rr_at_p97_5 < 1.05)
  ), by = .(geography, agegroup)]
}

erf_summary <- rbind(
  summarise_erf(erf, "Europe"),
  rbindlist(lapply(sort(unique(erf$region)), function(r) summarise_erf(erf[region == r], paste(r, "Europe")))),
  summarise_erf(erf[URAU_CODE == "ES001C"], "Madrid")
)
fwrite(erf_summary, file.path(diag_dir, "moderate_heat_erf_summary.csv"))

# Projected AN split in the first and last five-year periods. Population and
# rest deaths are irrelevant here; Object 1's heat causes are city-level ANs.
years_keep <- c(2020:2024, 2095:2099)
q <- open_dataset(obj1_file) %>%
  filter(ssp == !!as.integer(ssp_name), year %in% !!years_keep,
    cause %in% c("ModHeat", "ExtrHeat")) %>%
  group_by(city, region, scenario, year, cause) %>%
  summarise(AN = sum(deaths, na.rm = TRUE)) %>%
  collect() %>%
  as.data.table()
q[, period := fifelse(year <= 2024L, "2020-2024", "2095-2099")]
if (uniqueN(q$city) != uniqueN(meta$URAU_CODE) || any(!is.finite(q$AN)) || any(q$AN < 0)) {
  stop("Heat-AN diagnostic has incomplete city coverage or invalid AN values.", call. = FALSE)
}

region_an <- q[, .(AN = sum(AN)), by = .(region, scenario, period, year, cause)]
region_an[, geography := paste(region, "Europe")]
region_an <- region_an[, .(annual_mean_AN = mean(AN)), by = .(geography, scenario, period, cause)]
europe_an <- q[, .(AN = sum(AN)), by = .(scenario, period, year, cause)]
europe_an[, geography := "Europe"]
europe_an <- europe_an[, .(annual_mean_AN = mean(AN)), by = .(geography, scenario, period, cause)]
madrid_an <- q[city == "ES001C", .(AN = sum(AN)), by = .(scenario, period, year, cause)]
madrid_an[, geography := "Madrid"]
madrid_an <- madrid_an[, .(annual_mean_AN = mean(AN)), by = .(geography, scenario, period, cause)]
heat <- rbind(europe_an, region_an, madrid_an)
heat <- dcast(heat, geography + scenario + period ~ cause, value.var = "annual_mean_AN", fill = 0)
heat[, `:=`(
  total_heat_AN = ModHeat + ExtrHeat,
  moderate_share_of_heat_AN = fifelse(ModHeat + ExtrHeat > 0, ModHeat / (ModHeat + ExtrHeat), NA_real_)
)]
if (nrow(heat) != 6L * length(branch_levels) * 2L || any(!is.finite(heat$moderate_share_of_heat_AN))) {
  stop("Heat-AN summary grid is incomplete or invalid.", call. = FALSE)
}
fwrite(heat, file.path(diag_dir, "moderate_heat_attributable_summary.csv"))

# Simon's 2 Oct follow-up: the LE figure concerns change, not total burden.
# Keep both branch totals above, and make the with-minus-without change explicit.
heat_delta <- dcast(
  heat,
  geography + period ~ scenario,
  value.var = c("ModHeat", "ExtrHeat")
)
heat_delta[, `:=`(
  incremental_mod_heat_AN = ModHeat_with_cc - ModHeat_without_cc,
  incremental_extr_heat_AN = ExtrHeat_with_cc - ExtrHeat_without_cc
)]
heat_delta[, incremental_heat_AN := incremental_mod_heat_AN + incremental_extr_heat_AN]
heat_delta[, moderate_share_of_incremental_heat_AN := fifelse(
  incremental_heat_AN != 0,
  incremental_mod_heat_AN / incremental_heat_AN,
  NA_real_
)]
if (nrow(heat_delta) != 6L * 2L || any(!is.finite(heat_delta$moderate_share_of_incremental_heat_AN))) {
  stop("Incremental heat-AN summary grid is incomplete or invalid.", call. = FALSE)
}
fwrite(heat_delta, file.path(diag_dir, "moderate_heat_incremental_summary.csv"))

checks <- data.table(
  check_name = c("erf_city_age_coverage", "rr_finite_and_clamped", "object1_city_coverage", "heat_summary_grid"),
  status = "PASS",
  value = c(
    sprintf("%d rows, %d cities", nrow(erf), uniqueN(erf$URAU_CODE)),
    sprintf("RR range %.4f to %.4f", min(erf$rr_at_p97_5), max(erf$rr_at_p99)),
    sprintf("%d cities", uniqueN(q$city)),
    sprintf("%d rows", nrow(heat))
  )
)
fwrite(checks, file.path(diag_dir, "moderate_heat_checks.csv"))

p <- ggplot(heat[scenario == "with_cc"], aes(period, moderate_share_of_heat_AN, colour = geography, group = geography)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = c(
    "Europe" = "black", "Eastern Europe" = "#E6AB02", "Northern Europe" = "#1B9E77",
    "Southern Europe" = "#D95F02", "Western Europe" = "#7570B3", "Madrid" = "#E7298A"
  )) +
  labs(
    title = "Moderate heat as a share of projected heat-attributable deaths",
    subtitle = sprintf("%s, central ERF, fixed ERA5 1990-2019 thresholds", ssplabs[ssp_name]),
    x = NULL, y = "Moderate-heat share", colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(diag_dir, "moderate_heat_share.png"), p, width = 10, height = 6, dpi = 160)

heat_long <- melt(
  heat[period == "2095-2099"],
  id.vars = c("geography", "scenario", "period"),
  measure.vars = c("ModHeat", "ExtrHeat"),
  variable.name = "cause",
  value.name = "annual_mean_AN"
)
p_deaths <- ggplot(
  heat_long,
  aes(x = scenario, y = annual_mean_AN, fill = cause)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_wrap(~ geography, scales = "free_y", ncol = 3) +
  scale_x_discrete(labels = branch_labels) +
  scale_fill_manual(
    values = c(ModHeat = range_colors[["ModHeat"]], ExtrHeat = range_colors[["ExtrHeat"]]),
    labels = c(ModHeat = "Moderate heat", ExtrHeat = "Extreme heat")
  ) +
  labs(
    title = "Moderate- and extreme-heat attributable deaths",
    subtitle = sprintf("%s, annual mean in 2095-2099", ssplabs[ssp_name]),
    x = NULL, y = "Attributable deaths per year", fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(
  file.path(diag_dir, "moderate_extreme_heat_deaths_end_century.png"),
  p_deaths, width = 12, height = 7, dpi = 160
)

print(erf_summary, digits = 4)
print(heat, digits = 4)
print(heat_delta, digits = 4)
message("Saved moderate-heat diagnostic to ", diag_dir)
