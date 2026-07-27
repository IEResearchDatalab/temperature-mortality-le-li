#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

message("\n[13] Generating future projection validation outputs...")

input_dir <- trimws(Sys.getenv("INPUT_DIR", unset = "temp_results"))
output_dir <- trimws(Sys.getenv("OUTPUT_DIR", unset = "results/validation_future"))

expected_ranges <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
age_levels <- c("20-44", "45-64", "65-74", "75-84", "85+")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

projection_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
if (!length(projection_files)) {
  stop(sprintf("No projection files found in %s", input_dir), call. = FALSE)
}

city_meta <- fread("data/city_results.csv")
city_key <- unique(city_meta[, .(URAU_CODE, LABEL, CNTR_CODE, cntr_name, region, lat, lon, agegroup, agepop, death)])

infer_period <- function(year) {
  fifelse(
    year <= 2019L,
    sprintf("%ds", floor(year / 10) * 10),
    sprintf("%ds", floor(year / 10) * 10)
  )
}

message("Loading projection files from ", input_dir, " ...")

proj_long <- rbindlist(lapply(projection_files, function(f) {
  city_id <- sub("\\.rds$", "", basename(f))
  d <- as.data.table(readRDS(f))
  if (!nrow(d)) return(NULL)
  d[, URAU_CODE := city_id]
  d
}), use.names = TRUE, fill = TRUE)

if (!nrow(proj_long)) {
  stop("Projection files were present but no rows were loaded.", call. = FALSE)
}

proj_long[, `:=`(
  sim = as.integer(sim),
  period = infer_period(year),
  range = factor(range, levels = expected_ranges),
  agegroup = factor(agegroup, levels = age_levels)
)]

proj_long <- merge(proj_long, city_key, by = c("URAU_CODE", "agegroup"), all.x = TRUE)

if (anyNA(proj_long$agepop)) {
  stop("Some projection rows could not be matched to city metadata.", call. = FALSE)
}

proj_long[, thermal_class := fifelse(range %in% c("ExtrCold", "ModCold"), "Cold", "Heat")]
proj_long[, intensity_class := fifelse(range %in% c("ExtrCold", "ExtrHeat"), "Extreme", "Moderate")]

summary_overall <- proj_long[, .(
  an_mean = mean(an),
  an_p2.5 = quantile(an, 0.025),
  an_p97.5 = quantile(an, 0.975),
  rate_per_100k_mean = mean(an / agepop * 100000)
), by = .(gcm, ssp, period, year, thermal_class, intensity_class, range)]

summary_decade <- proj_long[, .(
  an_mean = mean(an),
  an_p2.5 = quantile(an, 0.025),
  an_p97.5 = quantile(an, 0.975),
  rate_per_100k_mean = mean(an / agepop * 100000)
), by = .(gcm, ssp, period, thermal_class, intensity_class, range)]

summary_city_decade <- proj_long[, .(
  an_mean = mean(an),
  an_p2.5 = quantile(an, 0.025),
  an_p97.5 = quantile(an, 0.975),
  rate_per_100k_mean = mean(an / agepop * 100000)
), by = .(URAU_CODE, LABEL, region, gcm, ssp, period, thermal_class, intensity_class, range)]

summary_city_thermal_decade <- proj_long[, .(
  an_mean = mean(an),
  an_p2.5 = quantile(an, 0.025),
  an_p97.5 = quantile(an, 0.975),
  rate_per_100k_mean = mean(an / agepop * 100000)
), by = .(URAU_CODE, LABEL, region, gcm, ssp, period, thermal_class)]

summary_age_decade <- proj_long[, .(
  an_mean = mean(an),
  an_p2.5 = quantile(an, 0.025),
  an_p97.5 = quantile(an, 0.975),
  rate_per_100k_mean = mean(an / agepop * 100000)
), by = .(agegroup, gcm, ssp, period, thermal_class, intensity_class, range)]

net_summary <- proj_long[, .(
  an_mean = mean(an),
  an_p2.5 = quantile(an, 0.025),
  an_p97.5 = quantile(an, 0.975)
), by = .(gcm, ssp, year, period, thermal_class)]
net_wide <- dcast(net_summary, gcm + ssp + year + period ~ thermal_class, value.var = c("an_mean", "an_p2.5", "an_p97.5"), fill = 0)
net_wide[, `:=`(
  net_mean = an_mean_Heat - an_mean_Cold,
  net_p2.5 = an_p2.5_Heat - an_p97.5_Cold,
  net_p97.5 = an_p97.5_Heat - an_p2.5_Cold
)]

period_order <- c(sprintf("%ds", seq(1980, 2090, by = 10)))
summary_decade[, period := factor(period, levels = period_order)]
summary_city_decade[, period := factor(period, levels = period_order)]
summary_age_decade[, period := factor(period, levels = period_order)]
net_wide[, period := factor(period, levels = period_order)]

fwrite(summary_overall, file.path(output_dir, "future_validation_yearly_summary.csv"))
fwrite(summary_decade, file.path(output_dir, "future_validation_decade_summary.csv"))
fwrite(summary_city_decade, file.path(output_dir, "future_validation_city_decade_summary.csv"))
fwrite(summary_city_thermal_decade, file.path(output_dir, "future_validation_city_thermal_decade_summary.csv"))
fwrite(summary_age_decade, file.path(output_dir, "future_validation_age_decade_summary.csv"))
fwrite(net_wide, file.path(output_dir, "future_validation_net_heat_minus_cold.csv"))

table_period <- summary_decade[period %in% c("1990s", "2010s", "2020s", "2050s", "2090s"), .(
  gcm,
  ssp,
  period,
  range,
  an_mean = round(an_mean, 3),
  an_p2.5 = round(an_p2.5, 3),
  an_p97.5 = round(an_p97.5, 3),
  rate_per_100k_mean = round(rate_per_100k_mean, 3)
)]
fwrite(table_period, file.path(output_dir, "table1_key_decade_range_summary.csv"))

city_change <- dcast(
  summary_city_thermal_decade[period %in% c("2020s", "2090s"), .(URAU_CODE, LABEL, region, gcm, ssp, period, thermal_class, an_mean)],
  URAU_CODE + LABEL + region + gcm + ssp + thermal_class ~ period,
  value.var = "an_mean"
)
if (all(c("2020s", "2090s") %in% names(city_change))) {
  city_change[, delta_2090s_vs_2020s := `2090s` - `2020s`]
  fwrite(city_change, file.path(output_dir, "table2_city_change_2090s_vs_2020s.csv"))
}

message("Creating Figure 1: time series by range...")

p_range <- ggplot(
  summary_overall,
  aes(x = year, y = an_mean, color = range, fill = range)
) +
  geom_ribbon(aes(ymin = an_p2.5, ymax = an_p97.5), alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 0.8) +
  facet_grid(ssp ~ gcm, scales = "free_y") +
  scale_color_manual(values = c(
    ExtrCold = "#2166ac",
    ModCold = "#67a9cf",
    ModHeat = "#f4a582",
    ExtrHeat = "#b2182b"
  )) +
  scale_fill_manual(values = c(
    ExtrCold = "#2166ac",
    ModCold = "#67a9cf",
    ModHeat = "#f4a582",
    ExtrHeat = "#b2182b"
  )) +
  labs(
    title = "Projected attributable deaths by temperature range",
    x = "Year",
    y = "Annual attributable deaths",
    color = "Range",
    fill = "Range"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(output_dir, "fig1_range_timeseries.png"), p_range, width = 10, height = 6.2, dpi = 220)
ggsave(file.path(output_dir, "fig1_range_timeseries.pdf"), p_range, width = 10, height = 6.2)

message("Creating Figure 2: net heat minus cold...")

p_net <- ggplot(net_wide, aes(x = year, y = net_mean)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_ribbon(aes(ymin = net_p2.5, ymax = net_p97.5), fill = "grey70", alpha = 0.35) +
  geom_line(color = "black", linewidth = 0.9) +
  facet_grid(ssp ~ gcm, scales = "free_y") +
  labs(
    title = "Net burden shift: heat minus cold",
    subtitle = "Positive values indicate heat attributable deaths exceed cold attributable deaths",
    x = "Year",
    y = "Annual attributable deaths"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(output_dir, "fig2_net_heat_minus_cold.png"), p_net, width = 10, height = 5.8, dpi = 220)
ggsave(file.path(output_dir, "fig2_net_heat_minus_cold.pdf"), p_net, width = 10, height = 5.8)

message("Creating Figure 3: age profile in the 2090s...")

age_2090s <- summary_age_decade[period == "2090s"]
p_age <- ggplot(age_2090s, aes(x = agegroup, y = rate_per_100k_mean, color = range, group = range)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_grid(ssp ~ gcm) +
  scale_color_manual(values = c(
    ExtrCold = "#2166ac",
    ModCold = "#67a9cf",
    ModHeat = "#f4a582",
    ExtrHeat = "#b2182b"
  )) +
  labs(
    title = "2090s age gradient by temperature range",
    x = "Age group",
    y = "Attributable deaths per 100,000",
    color = "Range"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(output_dir, "fig3_age_gradient_2090s.png"), p_age, width = 10, height = 5.8, dpi = 220)
ggsave(file.path(output_dir, "fig3_age_gradient_2090s.pdf"), p_age, width = 10, height = 5.8)

if (exists("city_change") && "delta_2090s_vs_2020s" %in% names(city_change)) {
  message("Creating Figure 4: city-level change from 2020s to 2090s...")

  city_change_plot <- city_change[thermal_class %in% c("Cold", "Heat")]
  p_city <- ggplot(city_change_plot, aes(x = reorder(LABEL, delta_2090s_vs_2020s), y = delta_2090s_vs_2020s, fill = thermal_class)) +
    geom_col() +
    facet_grid(ssp ~ gcm + thermal_class, scales = "free_x", space = "free_x") +
    coord_flip() +
    scale_fill_manual(values = c(Cold = "#67a9cf", Heat = "#ef8a62")) +
    labs(
      title = "City-level change from the 2020s to the 2090s",
      x = NULL,
      y = "Change in annual attributable deaths",
      fill = NULL
    ) +
    theme_minimal(base_size = 10)

  ggsave(file.path(output_dir, "fig4_city_change_2090s_vs_2020s.png"), p_city, width = 11, height = 7.5, dpi = 220)
  ggsave(file.path(output_dir, "fig4_city_change_2090s_vs_2020s.pdf"), p_city, width = 11, height = 7.5)
}

readme_lines <- c(
  "# Future Projection Validation Outputs",
  "",
  sprintf("- Input directory: `%s`", input_dir),
  sprintf("- Cities included: %d", uniqueN(proj_long$URAU_CODE)),
  sprintf("- GCMs included: %s", paste(sort(unique(proj_long$gcm)), collapse = ", ")),
  sprintf("- SSPs included: %s", paste(sort(unique(proj_long$ssp)), collapse = ", ")),
  sprintf("- Year span: %d-%d", min(proj_long$year), max(proj_long$year)),
  "- Output interpretation: these are temperature-only projections from `R_pipeline/03_attribution.R` with baseline population and deaths held constant.",
  "- Main tables:",
  "  - `table1_key_decade_range_summary.csv`: decade-level attributable burden by temperature range.",
  "  - `table2_city_change_2090s_vs_2020s.csv`: city-level change in heat/cold burden between early and late-century decades when both decades are present.",
  "- Main figures:",
  "  - `fig1_range_timeseries`: annual trajectory by range.",
  "  - `fig2_net_heat_minus_cold`: when heat burden overtakes cold burden.",
  "  - `fig3_age_gradient_2090s`: late-century age profile.",
  "  - `fig4_city_change_2090s_vs_2020s`: which cities gain or lose the most burden."
)
writeLines(readme_lines, file.path(output_dir, "README.md"))

message("Future projection validation outputs written to ", output_dir)