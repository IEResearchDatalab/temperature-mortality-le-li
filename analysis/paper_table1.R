#!/usr/bin/env Rscript

################################################################################
# Paper Table 1: SSP3 end-century burden and LE65 erosion ratios.
# Presentation only: reads completed pooled Part 03/04 outputs and Object 1.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(arrow)
  library(dplyr)
})

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
figure_root <- Sys.getenv("FIGURE_ROOT", "results/paper_figures/20261013")
out_dir <- file.path(figure_root, "table1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

range_levels <- c("ExtrCold", "ModCold", "ModHeat", "NormExtrHeat", "ExtrHeat")
heat_levels <- c("ModHeat", "NormExtrHeat", "ExtrHeat")
period_start <- function(x) as.integer(substr(x, 1L, 4L))

geo_root <- file.path(batch_root, "geographies", "ssp3")
geo_index <- rbindlist(list(
  data.table(level = "europe", id = "EUROPE", order_group = 1L),
  data.table(level = "region", id = basename(list.dirs(file.path(geo_root, "region"), recursive = FALSE)), order_group = 2L),
  data.table(level = "country", id = basename(list.dirs(file.path(geo_root, "country"), recursive = FALSE)), order_group = 3L)
), use.names = TRUE)
if (nrow(geo_index) != 35L) stop("Expected 35 pooled geographies, found ", nrow(geo_index), call. = FALSE)

summarise_geo <- function(level, id, order_group) {
  d <- file.path(geo_root, level, id)
  master_path <- file.path(d, "03_master_table.csv")
  levels_path <- file.path(d, "04_le_li_levels.csv")
  within_path <- file.path(d, "04_within_branch_period_decomposition.csv")
  if (!all(file.exists(c(master_path, levels_path, within_path)))) {
    stop("Incomplete Table 1 geography directory: ", d, call. = FALSE)
  }

  master <- fread(master_path)
  label <- unique(master$label)
  if (length(label) != 1L) stop("Non-unique geography label in ", master_path, call. = FALSE)
  end <- master[year %between% c(2095L, 2099L) & range %chin% range_levels]
  deaths <- end[, .(deaths = sum(an)), by = .(branch, range)]
  pop <- unique(end[, .(branch, year, age, pop)])[, .(person_years = sum(pop)), by = branch]
  totals <- deaths[, .(deaths = sum(deaths)), by = branch]
  totals <- merge(totals, pop, by = "branch")
  totals[, rate := 1e5 * deaths / person_years]
  with_total <- totals[branch == "with_cc"]
  without_total <- totals[branch == "without_cc"]
  with_heat <- deaths[branch == "with_cc" & range %chin% heat_levels, sum(deaths)]

  levels <- fread(levels_path)[gcm == "ENSEMBLE"]
  le_period <- levels[year %in% c(2020:2024, 2095:2099), .(
    LE65 = mean(LE65)
  ), by = .(branch, period = fifelse(year <= 2024L, "start", "end"))]
  le_wide <- dcast(le_period, branch ~ period, value.var = "LE65")
  le_wide[, gain := end - start]
  gain_with <- le_wide[branch == "with_cc", gain]
  gain_without <- le_wide[branch == "without_cc", gain]

  within <- fread(within_path)[gcm == "ENSEMBLE"]
  within[, `:=`(from_start = period_start(period_from), to_start = period_start(period_to))]
  heat_gain <- within[
    from_start >= 2020L & to_start <= 2095L & cause %chin% heat_levels,
    .(heat_gain = sum(le_contribution)), by = branch
  ]
  heat_gain_with <- heat_gain[branch == "with_cc", heat_gain]
  heat_gain_without <- heat_gain[branch == "without_cc", heat_gain]
  small_denominator <- gain_without < (1 / 12)

  data.table(
    geography = label, level = level, id = id, order_group = order_group,
    deaths_per_100k = with_total$rate,
    climate_attributable_share_pct = 100 * (with_total$deaths - without_total$deaths) / with_total$deaths,
    modheat_share_pct = 100 * deaths[branch == "with_cc" & range == "ModHeat", deaths] / with_heat,
    normextrheat_share_pct = 100 * deaths[branch == "with_cc" & range == "NormExtrHeat", deaths] / with_heat,
    le65_gain_without_cc_months = gain_without * 12,
    heat_only_erosion_ratio = fifelse(
      small_denominator, NA_real_, (heat_gain_without - heat_gain_with) / gain_without
    ),
    all_temperature_erosion_ratio = fifelse(
      small_denominator, NA_real_, (gain_without - gain_with) / gain_without
    ),
    small_denominator = small_denominator
  )
}

tab <- rbindlist(Map(summarise_geo, geo_index$level, geo_index$id, geo_index$order_group))
tab[order_group == 2L & !grepl("Europe$", geography), geography := paste0(geography, " Europe")]
setorder(tab, order_group, geography)
if (nrow(tab) != 35L || any(!is.finite(tab$deaths_per_100k))) {
  stop("Table 1 is incomplete or has invalid death rates.", call. = FALSE)
}
fwrite(tab, file.path(out_dir, "table1_ssp3_geographies.csv"))

# Independent Object 1 recomputation for Europe and Spain.
object1_path <- file.path(batch_root, "collected", "object1_dataset.parquet")
if (!file.exists(object1_path)) stop("Missing Object 1: ", object1_path, call. = FALSE)
obj <- open_dataset(object1_path) |>
  filter(ssp == 3L, year >= 2095L, year <= 2099L,
         scenario == "with_cc", cause %in% !!range_levels) |>
  collect() |>
  as.data.table()
independent_rate <- function(x) {
  1e5 * x[, sum(deaths)] / unique(x[, .(city, year, age, pop)])[, sum(pop)]
}
checks <- data.table(
  geography = c("Europe", "Spain"),
  object1_rate = c(independent_rate(obj), independent_rate(obj[country == "ES"]))
)
checks[, pooled_rate := c(
  tab[level == "europe", deaths_per_100k],
  tab[level == "country" & id == "ES", deaths_per_100k]
)]
checks[, abs_difference := abs(object1_rate - pooled_rate)]
fwrite(checks, file.path(out_dir, "table1_independent_checks.csv"))

# Formatted landscape table.
display <- tab[, .(
  Geography = geography,
  `Deaths per 100,000` = sprintf("%.1f", deaths_per_100k),
  `Climate-attributable share` = sprintf("%.1f%%", climate_attributable_share_pct),
  `Moderate-heat share` = sprintf("%.1f%%", modheat_share_pct),
  `Normal-extreme share` = sprintf("%.1f%%", normextrheat_share_pct),
  `Heat-only LE65 erosion` = fifelse(is.na(heat_only_erosion_ratio), "n/a", sprintf("%.3f", heat_only_erosion_ratio)),
  `All-temperature LE65 erosion` = fifelse(is.na(all_temperature_erosion_ratio), "n/a", sprintf("%.3f", all_temperature_erosion_ratio))
)]
long <- melt(display, id.vars = "Geography", variable.name = "column", value.name = "value")
long[, `:=`(
  row = match(Geography, display$Geography),
  col = match(column, names(display)[-1])
)]
x_positions <- c(4.3, 6.9, 9.8, 12.5, 15.1, 18.0)
long[, x := x_positions[col]]
geo_text <- data.table(Geography = display$Geography, row = seq_len(nrow(display)), x = 0.05)
headers <- data.table(
  label = c("Geography", names(display)[-1]),
  x = c(0.05, x_positions), y = 0,
  hjust = c(0, rep(0.5, length(x_positions)))
)
row_tiles <- data.table(
  row = seq_len(nrow(display)),
  fill = rep(c("#FFFFFF", "#F3F6F8"), length.out = nrow(display))
)
small_names <- tab[small_denominator == TRUE, geography]
footnote <- if (length(small_names)) {
  paste0("n/a: gain in LE65 without climate change is below one month for ",
         paste(small_names, collapse = ", "), ".")
} else {
  "No geography has a below-one-month denominator."
}

p <- ggplot() +
  geom_rect(data = row_tiles,
            aes(xmin = -0.2, xmax = 19.5, ymin = row - 0.48, ymax = row + 0.48, fill = fill),
            colour = "#D9D9D9", linewidth = 0.15, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(data = geo_text, aes(x, row, label = Geography), hjust = 0, size = 2.45) +
  geom_text(data = long, aes(x, row, label = value), size = 2.35) +
  geom_rect(aes(xmin = -0.2, xmax = 19.5, ymin = -0.65, ymax = 0.55),
            fill = "#244A66", colour = "#244A66") +
  geom_text(data = headers, aes(x, y, label = label, hjust = hjust), colour = "white",
            fontface = "bold", size = 2.35) +
  scale_y_reverse(limits = c(nrow(display) + 1.8, -1.2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.2, 19.5), clip = "off") +
  labs(
    title = "Table 1  Temperature burden and LE65 erosion in 2095-99 under SSP3-7.0",
    subtitle = "Europe, four regions and 30 countries; people aged 65+",
    caption = paste(
      "Deaths are central estimates per 100,000 person-years with climate change. Heat shares use all heat deaths.",
      "Erosion ratio = climate-related loss of the LE65 gain / LE65 gain without climate change; the heat-only numerator includes the three heat ranges.",
      footnote,
      "Static ERA5-Land 2000-2014 baseline; anchored moving P97.5 uses the 15 preceding years. Central estimate, no uncertainty intervals."
    )
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 14, colour = "black"),
    plot.subtitle = element_text(size = 10, colour = "#555555"),
    plot.caption = element_text(size = 7.5, hjust = 0, margin = margin(t = 10))
  )

ggsave(file.path(out_dir, "table1_ssp3_geographies.pdf"), p,
       width = 16, height = 11.5, units = "in", device = cairo_pdf)
message("Saved Table 1 and independent checks to ", out_dir)
