#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

message("\n[14] Generating lightweight future preview from completed city outputs...")

input_dir <- trimws(Sys.getenv("INPUT_DIR", unset = "temp_results_future"))
output_dir <- trimws(Sys.getenv("OUTPUT_DIR", unset = "results/validation_future_partial_quick"))
city_filter <- trimws(Sys.getenv("CITY_FILTER", unset = ""))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
city_plot_dir <- file.path(output_dir, "by_city")
dir.create(city_plot_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
if (nzchar(city_filter)) {
  wanted_cities <- trimws(strsplit(city_filter, ",", fixed = TRUE)[[1]])
  files <- files[sub("\\.rds$", "", basename(files)) %in% wanted_cities]
}
if (!length(files)) {
  stop(sprintf("No projection files found in %s", input_dir), call. = FALSE)
}

message("Found ", length(files), " completed city files.")

expected_ranges <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
expected_ssps <- c(1, 2, 3)

annual_city_list <- vector("list", length(files))

for (index in seq_along(files)) {
  file_path <- files[[index]]
  city_id <- sub("\\.rds$", "", basename(file_path))

  dt <- as.data.table(readRDS(file_path))
  dt_annual <- dt[, .(an = sum(an)), by = .(year, range, ssp, gcm, sim)]
  dt_city <- dt_annual[, .(an = mean(an)), by = .(year, range, ssp)]
  dt_city <- dt_city[CJ(year = sort(unique(year)), range = expected_ranges, ssp = expected_ssps, unique = TRUE), on = .(year, range, ssp)]
  dt_city[, URAU_CODE := city_id]

  annual_city_list[[index]] <- dt_city

  city_category <- dt_city[, .(
    annual_deaths = sum(an)
  ), by = .(
    year,
    ssp,
    category = fifelse(range %in% c("ExtrCold", "ModCold"), "Cold", "Heat")
  )]

  p_city_cold_heat <- ggplot(city_category, aes(x = year, y = annual_deaths, color = category, linewidth = factor(ssp))) +
    geom_line() +
    scale_color_manual(values = c(Cold = "#2b8cbe", Heat = "#de2d26")) +
    scale_linewidth_manual(values = c("1" = 0.6, "2" = 1.0, "3" = 1.5)) +
    labs(
      title = sprintf("%s: cold versus heat burden over time", city_id),
      x = "Year",
      y = "Annual attributable deaths",
      color = NULL,
      linewidth = "SSP"
    ) +
    theme_minimal(base_size = 10)

  ggsave(
    file.path(city_plot_dir, sprintf("%s_cold_vs_heat.png", city_id)),
    p_city_cold_heat,
    width = 8,
    height = 4.8,
    dpi = 180
  )

  p_city_ranges <- ggplot(dt_city, aes(x = year, y = an, color = range, linewidth = factor(ssp))) +
    geom_line() +
    scale_color_manual(values = c(
      ExtrCold = "#2166ac",
      ModCold = "#67a9cf",
      ModHeat = "#f4a582",
      ExtrHeat = "#b2182b"
    )) +
    scale_linewidth_manual(values = c("1" = 0.6, "2" = 1.0, "3" = 1.5)) +
    labs(
      title = sprintf("%s: attributable deaths by temperature range", city_id),
      x = "Year",
      y = "Annual attributable deaths",
      color = "Range",
      linewidth = "SSP"
    ) +
    theme_minimal(base_size = 10)

  ggsave(
    file.path(city_plot_dir, sprintf("%s_temperature_ranges.png", city_id)),
    p_city_ranges,
    width = 8.6,
    height = 5.2,
    dpi = 180
  )

  if (index %% 10L == 0L || index == length(files)) {
    message("  processed ", index, " / ", length(files), " files")
  }
}

annual_city <- rbindlist(annual_city_list, use.names = TRUE)
annual_city[, decade := floor(year / 10) * 10]

annual_summary <- annual_city[, .(
  completed_cities = uniqueN(URAU_CODE),
  annual_deaths = sum(an)
), by = .(year, range, ssp)]

decade_summary <- annual_city[, .(
  completed_cities = uniqueN(URAU_CODE),
  mean_annual_deaths = mean(an)
), by = .(ssp, decade, range)]
setorder(decade_summary, ssp, decade, range)

category_summary <- annual_summary[, .(
  annual_deaths = sum(annual_deaths)
), by = .(
  year,
  ssp,
  category = fifelse(range %in% c("ExtrCold", "ModCold"), "Cold", "Heat")
)]

fwrite(annual_summary, file.path(output_dir, "annual_summary_by_range.csv"))
fwrite(decade_summary, file.path(output_dir, "mean_annual_deaths_by_decade_range.csv"))
fwrite(category_summary, file.path(output_dir, "annual_cold_vs_heat_summary.csv"))

p_cold_heat <- ggplot(category_summary, aes(x = year, y = annual_deaths, color = category, linewidth = factor(ssp))) +
  geom_line() +
  scale_color_manual(values = c(Cold = "#2b8cbe", Heat = "#de2d26")) +
  scale_linewidth_manual(values = c("1" = 0.7, "2" = 1.05, "3" = 1.45)) +
  labs(
    title = "Completed-city preview: cold versus heat burden over time",
    subtitle = sprintf("Based on %d completed city outputs currently in %s", uniqueN(annual_city$URAU_CODE), input_dir),
    x = "Year",
    y = "Annual attributable deaths across completed cities",
    color = NULL,
    linewidth = "SSP"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(output_dir, "cold_vs_heat_over_time.png"), p_cold_heat, width = 8.5, height = 5.2, dpi = 220)

p_ranges <- ggplot(annual_summary, aes(x = year, y = annual_deaths, color = range, linewidth = factor(ssp))) +
  geom_line() +
  scale_color_manual(values = c(
    ExtrCold = "#2166ac",
    ModCold = "#67a9cf",
    ModHeat = "#f4a582",
    ExtrHeat = "#b2182b"
  )) +
  scale_linewidth_manual(values = c("1" = 0.6, "2" = 0.95, "3" = 1.35)) +
  labs(
    title = "Completed-city preview: attributable deaths by temperature range",
    subtitle = sprintf("Based on %d completed city outputs currently in %s", uniqueN(annual_city$URAU_CODE), input_dir),
    x = "Year",
    y = "Annual attributable deaths across completed cities",
    color = "Range",
    linewidth = "SSP"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(output_dir, "temperature_ranges_over_time.png"), p_ranges, width = 9.2, height = 7.2, dpi = 220)

readme_lines <- c(
  "# Lightweight Future Preview",
  "",
  sprintf("- Input directory: `%s`", input_dir),
  sprintf("- Completed city files included: %d", uniqueN(annual_city$URAU_CODE)),
  sprintf("- Year span: %d-%d", min(annual_city$year), max(annual_city$year)),
  "- These summaries are interim and only reflect cities whose `.rds` outputs already exist.",
  "- Each city is first collapsed across age groups, then averaged across GCMs and simulations within each SSP/year/range before being summed across completed cities.",
  "- Main outputs:",
  "  - `cold_vs_heat_over_time.png`",
  "  - `temperature_ranges_over_time.png`",
  "  - `mean_annual_deaths_by_decade_range.csv`",
  "  - `by_city/*_cold_vs_heat.png`",
  "  - `by_city/*_temperature_ranges.png`"
)
writeLines(readme_lines, file.path(output_dir, "README.md"))

message("Lightweight preview written to ", output_dir)