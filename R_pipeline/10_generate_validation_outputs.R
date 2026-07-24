#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sf)
})

message("\n[10] Generating validation figures and tables from baseline outputs...")

dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)

expected_ranges <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")

leap_days <- function(year) {
  ifelse((year %% 400L == 0L) | (year %% 4L == 0L & year %% 100L != 0L), 366L, 365L)
}

baseline_files <- list.files("temp_results_baseline", pattern = "\\.rds$", full.names = TRUE)
if (!length(baseline_files)) {
  stop("No baseline files found in temp_results_baseline/", call. = FALSE)
}

city_meta <- fread("data/city_results.csv")
masselot <- fread("references/2025-masselot-zenodo/results/cityage.csv")

city_meta_key <- unique(city_meta[, .(
  URAU_CODE,
  LABEL,
  CNTR_CODE,
  cntr_name,
  region,
  lon,
  lat,
  agegroup,
  agepop,
  death,
  stdrate_cold_est,
  stdrate_heat_est,
  excess_cold_est,
  excess_heat_est
)])

message("Loading baseline point-estimate results for ", length(baseline_files), " cities...")

baseline_long <- rbindlist(lapply(baseline_files, function(f) {
  city_id <- sub("\\.rds$", "", basename(f))
  d <- as.data.table(readRDS(f))
  d <- d[gcm == "ERA5" & sim == 0L]
  if (!nrow(d)) return(NULL)
  d[, URAU_CODE := city_id]
  d
}), use.names = TRUE, fill = TRUE)

if (!nrow(baseline_long)) {
  stop("No ERA5 point-estimate rows found in temp_results_baseline/", call. = FALSE)
}

message("Collapsing annual outputs to exact baseline-period attributable numbers...")

baseline_complete <- baseline_long[, .(an = sum(an)), by = .(URAU_CODE, agegroup, year, range)]
full_index <- baseline_complete[, CJ(
  URAU_CODE = unique(URAU_CODE),
  agegroup = unique(agegroup),
  year = sort(unique(year)),
  range = expected_ranges,
  unique = TRUE
), by = .(URAU_CODE, agegroup)]
baseline_complete <- baseline_complete[full_index, on = .(URAU_CODE, agegroup, year, range)]
baseline_complete[is.na(an), an := 0]
baseline_complete[, year_days := leap_days(year)]

baseline_range <- baseline_complete[, .(
  an_est = sum(an * year_days) / sum(year_days)
), by = .(URAU_CODE, agegroup, range)]

baseline_wide <- dcast(baseline_range, URAU_CODE + agegroup ~ range, value.var = "an_est", fill = 0)
baseline_wide[, `:=`(
  our_cold = ExtrCold + ModCold,
  our_heat = ModHeat + ExtrHeat,
  our_total = ExtrCold + ModCold + ModHeat + ExtrHeat
)]

merged <- merge(
  baseline_wide,
  masselot[, .(URAU_CODE, agegroup,
               masselot_cold = excess_cold_est,
               masselot_heat = excess_heat_est,
               masselot_total = excess_total_est)],
  by = c("URAU_CODE", "agegroup"),
  all.x = TRUE
)

merged <- merge(merged, city_meta_key, by = c("URAU_CODE", "agegroup"), all.x = TRUE)
merged[, `:=`(
  rate_ExtrCold = ExtrCold / agepop * 100000,
  rate_ModCold = ModCold / agepop * 100000,
  rate_ModHeat = ModHeat / agepop * 100000,
  rate_ExtrHeat = ExtrHeat / agepop * 100000,
  rate_total_extreme = (ExtrCold + ExtrHeat) / agepop * 100000,
  rate_total_moderate = (ModCold + ModHeat) / agepop * 100000,
  rate_cold = our_cold / agepop * 100000,
  rate_heat = our_heat / agepop * 100000,
  pct_total_mortality_extreme = (ExtrCold + ExtrHeat) / death * 100,
  extreme_to_moderate_ratio = fifelse((ModCold + ModHeat) == 0, NA_real_,
                                      (ExtrCold + ExtrHeat) / (ModCold + ModHeat))
)]

merged[, `:=`(
  diff_cold = our_cold - masselot_cold,
  diff_heat = our_heat - masselot_heat,
  error_cold_pct = fifelse(abs(masselot_cold) < 1e-12, 0, abs(our_cold - masselot_cold) / abs(masselot_cold) * 100),
  error_heat_pct = fifelse(abs(masselot_heat) < 1e-12, 0, abs(our_heat - masselot_heat) / abs(masselot_heat) * 100)
)]

fwrite(merged, "results/validation/validation_city_age_4range.csv")

message("Creating Figure 1: additivity check...")

p_cold <- ggplot(merged, aes(x = masselot_cold, y = our_cold)) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linewidth = 0.8) +
  geom_point(alpha = 0.25, size = 1.2, color = "steelblue4") +
  labs(
    title = "Panel A. Cold Additivity Check",
    subtitle = sprintf("R² = %.6f | n = %d city-age combinations", cor(merged$masselot_cold, merged$our_cold)^2, nrow(merged)),
    x = "Masselot 2023 total cold attributable deaths/year",
    y = "ExtrCold + ModCold attributable deaths/year"
  ) +
  theme_minimal(base_size = 11)

p_heat <- ggplot(merged, aes(x = masselot_heat, y = our_heat)) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linewidth = 0.8) +
  geom_point(alpha = 0.25, size = 1.2, color = "darkorange3") +
  labs(
    title = "Panel B. Heat Additivity Check",
    subtitle = sprintf("R² = %.6f | n = %d city-age combinations", cor(merged$masselot_heat, merged$our_heat)^2, nrow(merged)),
    x = "Masselot 2023 total heat attributable deaths/year",
    y = "ModHeat + ExtrHeat attributable deaths/year"
  ) +
  theme_minimal(base_size = 11)

ggsave("results/validation/fig1_additivity_cold.pdf", p_cold, width = 6.5, height = 5)
ggsave("results/validation/fig1_additivity_heat.pdf", p_heat, width = 6.5, height = 5)

message("Creating Figure 2 and Table 1: regional summaries...")

regional_summary <- merged[, .(
  cities = uniqueN(URAU_CODE),
  person_years = sum(agepop),
  ExtrCold = sum(ExtrCold),
  ModCold = sum(ModCold),
  ModHeat = sum(ModHeat),
  ExtrHeat = sum(ExtrHeat),
  our_cold = sum(our_cold),
  our_heat = sum(our_heat),
  masselot_cold = sum(masselot_cold),
  masselot_heat = sum(masselot_heat)
), by = region][order(factor(region, levels = c("Northern", "Western", "Eastern", "Southern")))]

regional_summary[, `:=`(
  rate_ExtrCold = ExtrCold / person_years * 100000,
  rate_ModCold = ModCold / person_years * 100000,
  rate_ModHeat = ModHeat / person_years * 100000,
  rate_ExtrHeat = ExtrHeat / person_years * 100000,
  rate_total_cold = our_cold / person_years * 100000,
  rate_total_heat = our_heat / person_years * 100000,
  rate_masselot_cold = masselot_cold / person_years * 100000,
  rate_masselot_heat = masselot_heat / person_years * 100000,
  error_cold_pct = fifelse(abs(masselot_cold) < 1e-12, 0, abs(our_cold - masselot_cold) / abs(masselot_cold) * 100),
  error_heat_pct = fifelse(abs(masselot_heat) < 1e-12, 0, abs(our_heat - masselot_heat) / abs(masselot_heat) * 100)
)]

regional_total <- regional_summary[, .(
  region = "Total",
  cities = sum(cities),
  person_years = sum(person_years),
  ExtrCold = sum(ExtrCold),
  ModCold = sum(ModCold),
  ModHeat = sum(ModHeat),
  ExtrHeat = sum(ExtrHeat),
  our_cold = sum(our_cold),
  our_heat = sum(our_heat),
  masselot_cold = sum(masselot_cold),
  masselot_heat = sum(masselot_heat)
)]

regional_total[, `:=`(
  rate_ExtrCold = ExtrCold / person_years * 100000,
  rate_ModCold = ModCold / person_years * 100000,
  rate_ModHeat = ModHeat / person_years * 100000,
  rate_ExtrHeat = ExtrHeat / person_years * 100000,
  rate_total_cold = our_cold / person_years * 100000,
  rate_total_heat = our_heat / person_years * 100000,
  rate_masselot_cold = masselot_cold / person_years * 100000,
  rate_masselot_heat = masselot_heat / person_years * 100000,
  error_cold_pct = fifelse(abs(masselot_cold) < 1e-12, 0, abs(our_cold - masselot_cold) / abs(masselot_cold) * 100),
  error_heat_pct = fifelse(abs(masselot_heat) < 1e-12, 0, abs(our_heat - masselot_heat) / abs(masselot_heat) * 100)
)]

table1 <- rbindlist(list(regional_summary, regional_total), use.names = TRUE, fill = TRUE)
table1_export <- table1[, .(
  Region = region,
  Cities = cities,
  ExtrCold = round(rate_ExtrCold, 2),
  ModCold = round(rate_ModCold, 2),
  ModHeat = round(rate_ModHeat, 2),
  ExtrHeat = round(rate_ExtrHeat, 2),
  TotalCold = round(rate_total_cold, 2),
  TotalHeat = round(rate_total_heat, 2),
  MasselotCold = round(rate_masselot_cold, 2),
  MasselotHeat = round(rate_masselot_heat, 2),
  ErrorColdPct = round(error_cold_pct, 4),
  ErrorHeatPct = round(error_heat_pct, 4)
)]
fwrite(table1_export, "results/validation/table1_validation_summary.csv")

regional_long <- melt(
  regional_summary[, .(region, rate_ExtrCold, rate_ModCold, rate_ModHeat, rate_ExtrHeat)],
  id.vars = "region",
  variable.name = "range",
  value.name = "rate_per_100k"
)
regional_long[, range := factor(
  sub("^rate_", "", range),
  levels = c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
)]

p_region <- ggplot(regional_long, aes(x = factor(region, levels = c("Northern", "Western", "Eastern", "Southern")),
                                      y = rate_per_100k,
                                      fill = range)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    ExtrCold = "#2c7fb8",
    ModCold = "#7fcdbb",
    ModHeat = "#fdae61",
    ExtrHeat = "#d7191c"
  )) +
  labs(
    title = "Regional Contributions by Temperature Range",
    x = NULL,
    y = "Attributable deaths per 100,000 person-years",
    fill = "Range"
  ) +
  theme_minimal(base_size = 11)

ggsave("results/validation/fig2_regional_bars.pdf", p_region, width = 8, height = 5.5)

message("Creating Figure 3 and Table 3: age gradients...")

age_summary <- merged[, .(
  person_years = sum(agepop),
  ExtrCold = sum(ExtrCold),
  ModCold = sum(ModCold),
  ModHeat = sum(ModHeat),
  ExtrHeat = sum(ExtrHeat)
), by = agegroup]
age_summary[, age_order := match(agegroup, c("20-44", "45-64", "65-74", "75-84", "85+"))]
setorder(age_summary, age_order)
age_summary[, `:=`(
  rate_ExtrCold = ExtrCold / person_years * 100000,
  rate_ModCold = ModCold / person_years * 100000,
  rate_ModHeat = ModHeat / person_years * 100000,
  rate_ExtrHeat = ExtrHeat / person_years * 100000,
  cold_heat_ratio = (ExtrCold + ModCold) / (ModHeat + ExtrHeat)
)]

table3_export <- age_summary[, .(
  AgeGroup = agegroup,
  ExtrCold = round(rate_ExtrCold, 2),
  ModCold = round(rate_ModCold, 2),
  ModHeat = round(rate_ModHeat, 2),
  ExtrHeat = round(rate_ExtrHeat, 2),
  ColdHeatRatio = round(cold_heat_ratio, 2)
)]
fwrite(table3_export, "results/validation/table3_age_gradient_4range.csv")

age_long <- melt(
  age_summary[, .(agegroup, rate_ExtrCold, rate_ModCold, rate_ModHeat, rate_ExtrHeat)],
  id.vars = "agegroup",
  variable.name = "range",
  value.name = "rate_per_100k"
)
age_long[, range := factor(
  sub("^rate_", "", range),
  levels = c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
)]

p_age <- ggplot(age_long, aes(x = factor(agegroup, levels = c("20-44", "45-64", "65-74", "75-84", "85+")),
                              y = rate_per_100k,
                              color = range,
                              group = range)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    ExtrCold = "#2c7fb8",
    ModCold = "#7fcdbb",
    ModHeat = "#fdae61",
    ExtrHeat = "#d7191c"
  )) +
  labs(
    title = "Age Gradient by Temperature Range",
    x = "Age group",
    y = "Attributable deaths per 100,000 person-years",
    color = "Range"
  ) +
  theme_minimal(base_size = 11)

ggsave("results/validation/fig3_age_gradient.pdf", p_age, width = 7.5, height = 5.5)

message("Creating Figure 4 and Table 2: city-level extreme burden...")

city_summary <- merged[, .(
  LABEL = first(LABEL),
  cntr_name = first(cntr_name),
  region = first(region),
  lon = first(lon),
  lat = first(lat),
  pop_total = sum(agepop),
  death_total = sum(death),
  ExtrCold = sum(ExtrCold),
  ModCold = sum(ModCold),
  ModHeat = sum(ModHeat),
  ExtrHeat = sum(ExtrHeat)
), by = URAU_CODE]

city_summary[, `:=`(
  total_extreme = ExtrCold + ExtrHeat,
  total_moderate = ModCold + ModHeat,
  total_extreme_rate = (ExtrCold + ExtrHeat) / pop_total * 100000,
  total_temp_rate = (ExtrCold + ExtrHeat + ModCold + ModHeat) / pop_total * 100000,
  extreme_share_pct = (ExtrCold + ExtrHeat) / (ExtrCold + ExtrHeat + ModCold + ModHeat) * 100,
  pct_total_mortality = (ExtrCold + ExtrHeat) / death_total * 100,
  extreme_moderate_ratio = fifelse((ModCold + ModHeat) == 0, NA_real_,
                                   (ExtrCold + ExtrHeat) / (ModCold + ModHeat))
)]

table2_export <- city_summary[order(-total_extreme_rate)][1:20, .(
  City = LABEL,
  Country = cntr_name,
  Region = region,
  ExtrCold = round(ExtrCold / pop_total * 100000, 2),
  ExtrHeat = round(ExtrHeat / pop_total * 100000, 2),
  TotalExtreme = round(total_extreme_rate, 2),
  PctOfTotalMortality = round(pct_total_mortality, 3),
  ExtremeModerateRatio = round(extreme_moderate_ratio, 3)
)]
fwrite(table2_export, "results/validation/table2_top20_cities_extreme_burden.csv")

map_bbox <- st_bbox(c(xmin = -12, xmax = 35, ymin = 35, ymax = 72), crs = st_crs(4326))
europe_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
europe_map <- st_crop(europe_map, map_bbox)

p_map <- ggplot() +
  geom_sf(data = europe_map, fill = "grey96", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = city_summary,
    aes(x = lon, y = lat, color = extreme_moderate_ratio, size = total_extreme_rate),
    alpha = 0.9
  ) +
  scale_color_gradientn(
    colours = c("#2166ac", "#f7f7f7", "#b2182b"),
    name = "Extreme / Moderate"
  ) +
  scale_size_continuous(name = "Extreme burden\nper 100,000", range = c(0.7, 4)) +
  coord_sf(xlim = c(-12, 35), ylim = c(35, 72), expand = FALSE) +
  labs(
    title = "Extreme-to-Moderate Temperature Burden by City",
    subtitle = "Cities are overlaid on a Europe basemap; color shows composition, size shows extreme burden",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.2),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave("results/validation/fig4_extreme_moderate_ratio_map.pdf", p_map, width = 8.5, height = 6.5)

message("Creating Figure 5: why disaggregation matters at the city level...")

p_extreme_share <- ggplot(
  city_summary,
  aes(x = total_temp_rate, y = extreme_share_pct, color = region)
) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = c(
    Northern = "#1f78b4",
    Western = "#33a02c",
    Eastern = "#e31a1c",
    Southern = "#ff7f00"
  )) +
  labs(
    title = "Why Disaggregation Matters: Same Total Burden, Different Extreme Share",
    subtitle = "Cities with similar total temperature mortality can have very different fractions driven by extremes",
    x = "Total temperature-attributable deaths per 100,000 person-years",
    y = "Extreme share of total temperature burden (%)",
    color = "Region"
  ) +
  theme_minimal(base_size = 11)

ggsave("results/validation/fig5_city_extreme_share_vs_total.pdf", p_extreme_share, width = 8, height = 5.5)

message("Creating Figure 6: extreme shares by age group...")

age_share <- age_summary[, .(
  agegroup,
  cold_extreme_share = ExtrCold / (ExtrCold + ModCold) * 100,
  heat_extreme_share = ExtrHeat / (ModHeat + ExtrHeat) * 100
)]

age_share_long <- melt(
  age_share,
  id.vars = "agegroup",
  variable.name = "metric",
  value.name = "share_pct"
)
age_share_long[, metric := factor(
  metric,
  levels = c("cold_extreme_share", "heat_extreme_share"),
  labels = c("Extreme share within cold", "Extreme share within heat")
)]

p_age_share <- ggplot(
  age_share_long,
  aes(x = factor(agegroup, levels = c("20-44", "45-64", "65-74", "75-84", "85+")),
      y = share_pct,
      color = metric,
      group = metric)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "Extreme share within cold" = "#2c7fb8",
    "Extreme share within heat" = "#d7191c"
  )) +
  labs(
    title = "Extreme Shares by Age Group",
    subtitle = "Disaggregation shows how much of the cold and heat burden is driven by rare extremes",
    x = "Age group",
    y = "Extreme share within category (%)",
    color = NULL
  ) +
  theme_minimal(base_size = 11)

ggsave("results/validation/fig6_age_extreme_shares.pdf", p_age_share, width = 7.5, height = 5.5)

message("Writing summary metrics...")

summary_lines <- c(
  "# Validation Outputs Summary",
  "",
  sprintf("- Cities processed: %d", uniqueN(merged$URAU_CODE)),
  sprintf("- City-age combinations: %d", nrow(merged)),
  sprintf("- Cold R^2: %.6f", cor(merged$masselot_cold, merged$our_cold)^2),
  sprintf("- Heat R^2: %.6f", cor(merged$masselot_heat, merged$our_heat)^2),
  sprintf("- Mean cold error (%%): %.6f", mean(merged$error_cold_pct)),
  sprintf("- Mean heat error (%%): %.6f", mean(merged$error_heat_pct)),
  "",
  "## Files",
  "- fig1_additivity_cold.pdf",
  "- fig1_additivity_heat.pdf",
  "- fig2_regional_bars.pdf",
  "- fig3_age_gradient.pdf",
  "- fig4_extreme_moderate_ratio_map.pdf",
  "- fig5_city_extreme_share_vs_total.pdf",
  "- fig6_age_extreme_shares.pdf",
  "- table1_validation_summary.csv",
  "- table2_top20_cities_extreme_burden.csv",
  "- table3_age_gradient_4range.csv",
  "- validation_city_age_4range.csv"
)
writeLines(summary_lines, "results/validation/README.md")

message("✅ Validation outputs generated in results/validation/")