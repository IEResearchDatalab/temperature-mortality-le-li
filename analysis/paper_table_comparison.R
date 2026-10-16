#!/usr/bin/env Rscript

################################################################################
# Cross-scenario paper table: end-century burden and LE65 erosion ratios.
# Presentation only: reads completed pooled Part 03/04 outputs.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
figure_root <- Sys.getenv("FIGURE_ROOT", "results/paper_figures/20261013")
ssp_set <- as.integer(strsplit(Sys.getenv("SSP_SET", "1,2,3"), ",", fixed = TRUE)[[1]])
out_dir <- file.path(figure_root, "table1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!length(ssp_set) || anyNA(ssp_set) || any(!ssp_set %in% 1:3)) {
  stop("SSP_SET must be a comma-separated subset of 1,2,3.", call. = FALSE)
}

ssp_names <- c(`1` = "SSP1-2.6", `2` = "SSP2-4.5", `3` = "SSP3-7.0")
ssp_order <- c(1L, 2L, 3L)
range_levels <- c("ExtrCold", "ModCold", "ModHeat", "NormExtrHeat", "ExtrHeat")
heat_levels <- c("ModHeat", "NormExtrHeat", "ExtrHeat")
period_start <- function(x) as.integer(substr(x, 1L, 4L))

geography_index <- function(ssp) {
  root <- file.path(batch_root, "geographies", paste0("ssp", ssp))
  ans <- rbindlist(list(
    data.table(level = "europe", id = "EUROPE", order_group = 1L),
    data.table(level = "region", id = basename(list.dirs(file.path(root, "region"), recursive = FALSE)), order_group = 2L),
    data.table(level = "country", id = basename(list.dirs(file.path(root, "country"), recursive = FALSE)), order_group = 3L)
  ))
  if (nrow(ans) != 35L) {
    stop("Expected 35 pooled geographies for SSP", ssp, ", found ", nrow(ans), call. = FALSE)
  }
  ans
}

summarise_geo <- function(ssp, level, id, order_group) {
  d <- file.path(batch_root, "geographies", paste0("ssp", ssp), level, id)
  paths <- file.path(d, c(
    "03_master_table.csv", "04_le_li_levels.csv",
    "04_within_branch_period_decomposition.csv"
  ))
  if (!all(file.exists(paths))) stop("Incomplete pooled geography: ", d, call. = FALSE)

  master <- fread(paths[1])
  label <- unique(master$label)
  if (length(label) != 1L) stop("Non-unique geography label: ", d, call. = FALSE)
  end <- master[year %between% c(2095L, 2099L) & range %chin% range_levels]
  deaths <- end[, .(deaths = sum(an)), by = .(branch, range)]
  pop <- unique(end[, .(branch, year, age, pop)])[
    , .(person_years = sum(pop)), by = branch
  ]
  totals <- merge(deaths[, .(deaths = sum(deaths)), by = branch], pop, by = "branch")
  totals[, rate := 1e5 * deaths / person_years]
  with_total <- totals[branch == "with_cc"]
  without_total <- totals[branch == "without_cc"]
  with_heat <- deaths[branch == "with_cc" & range %chin% heat_levels, sum(deaths)]

  levels <- fread(paths[2])[gcm == "ENSEMBLE"]
  le_period <- levels[year %in% c(2020:2024, 2095:2099),
                      .(LE65 = mean(LE65)),
                      by = .(branch, period = fifelse(year <= 2024L, "start", "end"))]
  le_wide <- dcast(le_period, branch ~ period, value.var = "LE65")
  le_wide[, gain := end - start]
  gain_with <- le_wide[branch == "with_cc", gain]
  gain_without <- le_wide[branch == "without_cc", gain]

  within <- fread(paths[3])[gcm == "ENSEMBLE"]
  within[, `:=`(
    from_start = period_start(period_from),
    to_start = period_start(period_to)
  )]
  heat_gain <- within[
    from_start >= 2020L & to_start <= 2095L & cause %chin% heat_levels,
    .(heat_gain = sum(le_contribution)), by = branch
  ]
  heat_gain_with <- heat_gain[branch == "with_cc", heat_gain]
  heat_gain_without <- heat_gain[branch == "without_cc", heat_gain]
  small_denominator <- gain_without < (1 / 12)

  data.table(
    ssp = ssp,
    scenario = ssp_names[as.character(ssp)],
    geography = label,
    level = level,
    id = id,
    order_group = order_group,
    deaths_per_100k = with_total$rate,
    climate_attributable_share_pct =
      100 * (with_total$deaths - without_total$deaths) / with_total$deaths,
    modheat_share_pct =
      100 * deaths[branch == "with_cc" & range == "ModHeat", deaths] / with_heat,
    normextrheat_share_pct =
      100 * deaths[branch == "with_cc" & range == "NormExtrHeat", deaths] / with_heat,
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

all_tables <- rbindlist(lapply(ssp_set, function(ssp) {
  idx <- geography_index(ssp)
  rbindlist(Map(
    function(level, id, order_group) summarise_geo(ssp, level, id, order_group),
    idx$level, idx$id, idx$order_group
  ))
}))
all_tables[order_group == 2L & !grepl("Europe$", geography),
           geography := paste0(geography, " Europe")]
all_tables[, scenario := factor(scenario, levels = ssp_names[as.character(ssp_order)])]
setorder(all_tables, scenario, order_group, geography)

if (any(!is.finite(all_tables$deaths_per_100k))) {
  stop("Invalid death rates in comparison table.", call. = FALSE)
}

scenario_tag <- paste0("ssp", paste(sort(ssp_set), collapse = "_ssp"))
fwrite(all_tables, file.path(out_dir, paste0("table1_", scenario_tag, "_geographies.csv")))
for (ssp in ssp_set) {
  current_ssp <- ssp
  fwrite(all_tables[ssp == current_ssp],
         file.path(out_dir, paste0("table1_ssp", ssp, "_geographies.csv")))
}

make_scenario_pdf <- function(tab, ssp) {
  display <- tab[, .(
    Geography = geography,
    `Deaths / 100k` = sprintf("%.1f", deaths_per_100k),
    `Climate share` = sprintf("%.1f%%", climate_attributable_share_pct),
    `Moderate heat` = sprintf("%.1f%%", modheat_share_pct),
    `Normal extreme` = sprintf("%.1f%%", normextrheat_share_pct),
    `Heat LE65 erosion` = fifelse(is.na(heat_only_erosion_ratio), "n/a", sprintf("%.3f", heat_only_erosion_ratio)),
    `All-temp LE65 erosion` = fifelse(is.na(all_temperature_erosion_ratio), "n/a", sprintf("%.3f", all_temperature_erosion_ratio))
  )]
  long <- melt(display, id.vars = "Geography", variable.name = "column", value.name = "value")
  long[, `:=`(row = match(Geography, display$Geography), col = match(column, names(display)[-1]))]
  x_positions <- seq(4.4, 19.0, length.out = ncol(display) - 1L)
  long[, x := x_positions[col]]
  geo_text <- data.table(Geography = display$Geography, row = seq_len(nrow(display)), x = 0.05)
  headers <- data.table(label = c("Geography", names(display)[-1]),
                        x = c(0.05, x_positions), y = 0,
                        hjust = c(0, rep(0.5, length(x_positions))))
  row_tiles <- data.table(row = seq_len(nrow(display)),
                          fill = rep(c("#FFFFFF", "#F3F6F8"), length.out = nrow(display)))
  patch_note <- if (ssp == 1L) {
    "SSP1 uses the documented 0.9995 replacement for survival ratios rounded to exactly one."
  } else if (ssp == 2L) {
    "SSP2 includes 12 Swiss cities under the documented mortality-shape warning."
  } else ""

  p <- ggplot() +
    geom_rect(data = row_tiles,
              aes(xmin = -0.2, xmax = 20.5, ymin = row - 0.48, ymax = row + 0.48, fill = fill),
              colour = "#D9D9D9", linewidth = 0.15, show.legend = FALSE) +
    scale_fill_identity() +
    geom_text(data = geo_text, aes(x, row, label = Geography), hjust = 0, size = 2.45) +
    geom_text(data = long, aes(x, row, label = value), size = 2.3) +
    geom_rect(aes(xmin = -0.2, xmax = 20.5, ymin = -0.65, ymax = 0.55),
              fill = "#244A66", colour = "#244A66") +
    geom_text(data = headers, aes(x, y, label = label, hjust = hjust),
              colour = "white", fontface = "bold", size = 2.15) +
    scale_y_reverse(limits = c(nrow(display) + 1.8, -1.2), expand = c(0, 0)) +
    coord_cartesian(xlim = c(-0.2, 20.5), clip = "off") +
    labs(
      title = paste0("Table 1 — Temperature burden and LE65 erosion in 2095-99 under ", ssp_names[as.character(ssp)]),
      subtitle = "Europe, four regions and 30 countries; people aged 65+",
      caption = paste(
        "Central estimates; no uncertainty intervals.", patch_note,
        "Static ERA5-Land 2000-2014 baseline; anchored moving P97.5 uses the 15 preceding years."
      )
    ) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, colour = "#555555"),
          plot.caption = element_text(size = 7.5, hjust = 0, margin = margin(t = 10)))

  ggsave(file.path(out_dir, paste0("table1_ssp", ssp, "_geographies.pdf")), p,
         width = 17, height = 11.5, units = "in", device = cairo_pdf)
}

for (ssp in ssp_set) {
  current_ssp <- ssp
  make_scenario_pdf(all_tables[ssp == current_ssp], ssp)
}

if (all(1:3 %in% ssp_set)) {
  main <- all_tables[order_group <= 2L]
  main[, metric_order := match(ssp, ssp_order)]
  setorder(main, order_group, geography, metric_order)
  display <- main[, .(
    Geography = geography,
    Scenario = as.character(scenario),
    `Deaths / 100k` = sprintf("%.1f", deaths_per_100k),
    `Climate share` = sprintf("%.1f%%", climate_attributable_share_pct),
    `Normal extreme` = sprintf("%.1f%%", normextrheat_share_pct),
    `All-temp erosion` = fifelse(is.na(all_temperature_erosion_ratio), "n/a", sprintf("%.3f", all_temperature_erosion_ratio))
  )]
  display[, row := .I]
  long <- melt(display, id.vars = c("Geography", "Scenario", "row"),
               variable.name = "column", value.name = "value")
  x_positions <- c(0.05, 4.2, 7.0, 10.2, 13.4, 16.5)
  long[, x := x_positions[match(column, names(display)[3:6]) + 2L]]
  labels <- display[, .(Geography, Scenario, row)]
  row_tiles <- labels[, .(row, fill = rep(c("#FFFFFF", "#F3F6F8", "#E8EEF2"), 5L))]
  headers <- data.table(
    label = c("Geography", "Scenario", names(display)[3:6]),
    x = x_positions, y = 0,
    hjust = c(0, 0, rep(0.5, 4L))
  )
  p <- ggplot() +
    geom_rect(data = row_tiles,
              aes(xmin = -0.2, xmax = 18.5, ymin = row - 0.48, ymax = row + 0.48, fill = fill),
              colour = "#D9D9D9", linewidth = 0.15, show.legend = FALSE) +
    scale_fill_identity() +
    geom_text(data = labels, aes(0.05, row, label = Geography), hjust = 0, size = 2.8) +
    geom_text(data = labels, aes(4.2, row, label = Scenario), hjust = 0, size = 2.8) +
    geom_text(data = long, aes(x, row, label = value), size = 2.8) +
    geom_rect(aes(xmin = -0.2, xmax = 18.5, ymin = -0.65, ymax = 0.55),
              fill = "#244A66", colour = "#244A66") +
    geom_text(data = headers, aes(x, y, label = label, hjust = hjust),
              colour = "white", fontface = "bold", size = 2.6) +
    scale_y_reverse(limits = c(nrow(display) + 2.2, -1.2), expand = c(0, 0)) +
    coord_cartesian(xlim = c(-0.2, 18.5), clip = "off") +
    labs(
      title = "SSP1-2.6, SSP2-4.5 and SSP3-7.0 comparison in 2095-99",
      subtitle = "Europe and four regions; central estimates for people aged 65+",
      caption = paste(
        "SSP1 uses the documented 0.9995 survival-ratio replacement; SSP2 retains its mortality-shape warning.",
        "Central estimates, no uncertainty intervals. Full 35-geography values are supplied in the source CSV."
      )
    ) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, colour = "#555555"),
          plot.caption = element_text(size = 7.5, hjust = 0, margin = margin(t = 10)))
  ggsave(file.path(out_dir, "table1_ssp1_ssp2_ssp3_comparison.pdf"), p,
         width = 17, height = 8.5, units = "in", device = cairo_pdf)
}

message("Saved scenario-specific and cross-scenario tables to ", out_dir)
