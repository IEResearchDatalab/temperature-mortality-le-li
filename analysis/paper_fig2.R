#!/usr/bin/env Rscript

################################################################################
# Figure 2: contributions to the gain in LE65 by age and temperature range.
# Requires: first run pipeline Parts 00-04 for all SSPs and pool the requested
# geography with Parts 03-04. For every SSP, read
# <batch>/geographies/ssp<SSP>/<level>/<id>/03_master_table.csv and
# 04_within_branch_period_decomposition.csv. The script sources
# analysis/decade_decomposition.R, analysis/paper_helpers.R and
# analysis/output_names.R automatically.
# Produces: FIGURE_ROOT/{figures,figures_pdf,data}/figure02_* and decade checks.
# How to run: BATCH_ROOT=results/europe/heatsplit_20261013 \
#   FIGURE_ROOT=results/paper_output_new GEO_LEVEL=europe GEO_ID=EUROPE \
#   SSP_SET=1,2,3 Rscript analysis/paper_fig2.R
# Notes: SSP1 replacement and SSP1/SSP2 shape-check bypasses are exploratory;
#   see the caveats in README.md.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("analysis/output_names.R")
source("analysis/decade_decomposition.R")
source("analysis/paper_helpers.R")

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
figure_root <- Sys.getenv("FIGURE_ROOT", "results/figures_simon_20261123")
geo_level <- Sys.getenv("GEO_LEVEL", "europe")
geo_id <- Sys.getenv("GEO_ID", "EUROPE")
fpath <- function(filename) figure_file(figure_root, filename, geo_id)
data_dir <- file.path(figure_root, "data")

range_levels <- c("ExtrCold", "ModCold", "ModHeat", "NormExtrHeat", "ExtrHeat")
range_labels <- c(
  ExtrCold = "Extreme cold", ModCold = "Moderate cold",
  ModHeat = "Moderate heat", NormExtrHeat = "Normal extreme heat",
  ExtrHeat = "Extreme heat"
)
range_colours <- c(
  ExtrCold = "#24558A", ModCold = "#86BCE2",
  ModHeat = "#F3C623", NormExtrHeat = "#F06B18", ExtrHeat = "#AD1A2B"
)
ssp_labels <- c(`1` = "SSP1–2.6", `2` = "SSP2–4.5", `3` = "SSP3–7.0")
ssp_order <- c("1", "2", "3")
age_levels <- c("65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100+")
ssp_set <- as.integer(strsplit(Sys.getenv("SSP_SET", "1,2,3"), ",", fixed = TRUE)[[1]])
if (!length(ssp_set) || anyNA(ssp_set) || any(!ssp_set %in% 1:3)) {
  stop("SSP_SET must be a comma-separated subset of 1,2,3.", call. = FALSE)
}
exploratory_style_value <- Sys.getenv("EXPLORATORY_SCENARIO_STYLE", "0")
if (!exploratory_style_value %in% c("0", "1")) {
  stop("EXPLORATORY_SCENARIO_STYLE must be 0 or 1.", call. = FALSE)
}
exploratory_style <- exploratory_style_value == "1"

check_stack_order <- function() {
  test <- data.table(
    x = rep(c("positive", "negative"), each = length(range_levels)),
    cause = factor(rep(range_levels, 2L), levels = range_levels),
    value = c(seq_len(length(range_levels)), -seq_len(length(range_levels)))
  )
  pos <- ggplot_build(
    ggplot(test[value > 0], aes(x, value, fill = cause)) +
      geom_col(position = position_stack(reverse = FALSE)) +
      scale_fill_manual(values = range_colours)
  )$data[[1]]
  neg <- ggplot_build(
    ggplot(test[value < 0], aes(x, value, fill = cause)) +
      geom_col(position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = range_colours)
  )$data[[1]]
  positive_total <- sum(test[value > 0, value])
  negative_total <- sum(test[value < 0, value])
  if (pos$ymin[which(pos$fill == range_colours[["ExtrHeat"]])] != 0 ||
      pos$ymax[which(pos$fill == range_colours[["ExtrCold"]])] != positive_total ||
      neg$ymax[which(neg$fill == range_colours[["ExtrCold"]])] != 0 ||
      neg$ymin[which(neg$fill == range_colours[["ExtrHeat"]])] != negative_total) {
    stop("Figure 2 positive/negative stack order does not run hottest to coldest away from zero.", call. = FALSE)
  }
}
check_stack_order()

read_ssp <- function(ssp) {
  path <- file.path(batch_root, "geographies", paste0("ssp", ssp), geo_level,
                    geo_id, "04_within_branch_period_decomposition.csv")
  if (!file.exists(path)) stop("Missing within-branch decomposition: ", path, call. = FALSE)
  x <- fread(path)[gcm == "ENSEMBLE" & branch == "with_cc" & cause %chin% range_levels]
  if (!nrow(x) || length(setdiff(range_levels, unique(x$cause)))) {
    stop("Incomplete with-CC ensemble cause domain in ", path, call. = FALSE)
  }
  x[, `:=`(
    ssp = ssp,
    ssp_label = ssp_labels[as.character(ssp)],
    from_start = period_start(period_from),
    to_start = period_start(period_to),
    age_band = age_band(age)
  )]
  x
}

period_data <- rbindlist(lapply(ssp_set, read_ssp))
period_data[, ssp_label := factor(ssp_label, levels = ssp_labels[ssp_order])]
geo_labels <- unique(period_data$label)
if (length(geo_labels) != 1L) stop("Geography label is not unique.", call. = FALSE)
geo_label <- geo_labels[[1]]

decade_result <- compute_decade_decomposition(
  batch_root = batch_root,
  geo_level = geo_level,
  geo_id = geo_id,
  ssp_set = ssp_set,
  output_dir = data_dir,
  n_horiuchi = 400L
)
decade_data <- decade_result$data[cause %chin% range_levels]
decade_data[, `:=`(
  ssp_label = factor(ssp_labels[as.character(ssp)], levels = ssp_labels[ssp_order]),
  decade_label = factor(
    sprintf("%d–%02d", decade_to, (decade_to + 9L) %% 100L),
    levels = sprintf("%d–%02d", seq(2030L, 2090L, 10L), seq(39L, 99L, 10L))
  ),
  cause = factor(cause, levels = range_levels)
)]
when <- decade_data[, .(value = sum(le_contribution_months)),
                    by = .(ssp_label, decade_label, cause)]
when_net <- when[, .(value = sum(value)), by = .(ssp_label, decade_label)]
net_check <- merge(
  when_net,
  decade_data[, .(
    displayed_sum = sum(le_contribution_months),
    temperature_ranges = uniqueN(cause)
  ), by = .(ssp_label, decade_label)],
  by = c("ssp_label", "decade_label")
)
net_check[, series := "Cold and heat together"]
if (any(abs(net_check$value - net_check$displayed_sum) > 1e-12) ||
    any(net_check$temperature_ranges != length(range_levels))) {
  stop("Figure 2 diamonds do not equal the five displayed temperature ranges.", call. = FALSE)
}

plot_stacked <- function(data, xvar, net, x_label, y_label, show_strips, limits = NULL) {
  positive <- data[value > 0]
  negative <- data[value < 0]
  bar_mapping <- if (exploratory_style) {
    aes(x = .data[[xvar]], y = value, fill = cause, alpha = ssp_label)
  } else {
    aes(x = .data[[xvar]], y = value, fill = cause)
  }
  base <- ggplot() +
    geom_col(data = positive, mapping = bar_mapping,
             width = 0.82, position = position_stack(reverse = FALSE)) +
    geom_col(data = negative, mapping = bar_mapping,
             width = 0.82, position = position_stack(reverse = TRUE)) +
    geom_point(data = net, aes(x = .data[[xvar]], y = value, shape = "Cold and heat together"),
               inherit.aes = FALSE, size = 2.2, colour = "black") +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    facet_grid(ssp_label ~ ., scales = "fixed", switch = "y") +
    scale_fill_manual(values = range_colours, labels = range_labels, drop = FALSE) +
    scale_shape_manual(
      name = NULL, values = c("Cold and heat together" = 18),
      breaks = "Cold and heat together"
    ) +
    scale_y_continuous(
      breaks = function(limits) {
        breaks <- scales::breaks_extended(n = 6)(limits)
        if (limits[[1]] < 0 && !any(breaks < 0)) {
          breaks <- c(breaks, signif(limits[[1]], 1))
        }
        sort(unique(breaks))
      },
      labels = scales::label_number(style_negative = "minus")
    ) +
    labs(x = x_label, y = y_label, fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = if (show_strips) element_text(face = "bold", angle = 0, size = 11)
      else element_blank(),
      legend.position = "bottom"
    )
  if (!is.null(limits)) base <- base + coord_cartesian(ylim = limits)
  if (exploratory_style) {
    base <- base + scale_alpha_manual(
      values = c("SSP1–2.6" = 0.48, "SSP2–4.5" = 0.48, "SSP3–7.0" = 1),
      guide = "none"
    )
  }
  base
}

p_when <- plot_stacked(
  when, "decade_label", when_net, "10-year period", "Months of LE65 per decade", TRUE
) +
  labs(title = "a. By 10-year step") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8))

make_age_panel <- function(start_period, end_period, title) {
  years <- end_period - start_period
  z <- period_data[
    from_start >= start_period & to_start <= end_period,
    .(value = sum(le_contribution) * 12 * 10 / years),
    by = .(ssp_label, age_band, cause)
  ]
  z[, `:=`(
    age_band = factor(age_band, levels = age_levels),
    cause = factor(cause, levels = range_levels)
  )]
  net <- z[, .(value = sum(value)), by = .(ssp_label, age_band)]
  list(data = z, net = net, title = title)
}

first <- make_age_panel(2020L, 2055L, "b. By age group, 2020-24 to 2055-59")
second <- make_age_panel(2055L, 2095L, "c. By age group, 2055-59 to 2095-99")
age_limits <- range(c(0, first$data$value, first$net$value, second$data$value, second$net$value))
age_limits <- age_limits + c(-1, 1) * diff(age_limits) * 0.06

p_first <- plot_stacked(
  first$data, "age_band", first$net, "Age group", "Months of LE65 per decade",
  FALSE, age_limits
) +
  labs(title = first$title) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
        legend.position = "none")
p_second <- plot_stacked(
  second$data, "age_band", second$net, "Age group", "Months of LE65 per decade",
  FALSE, age_limits
) +
  labs(title = second$title) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
        legend.position = "none")

p <- (p_when | p_first | p_second) +
  plot_layout(widths = c(1.8, 1, 1), guides = "collect") +
  plot_annotation(
    caption = if (exploratory_style) {
      "SSP1 and SSP2 are exploratory; lighter bars mark these scenarios."
    } else {
      NULL
    },
    theme = theme(legend.position = "bottom")
  )

stem <- "figure02_le65_contributions"
ggsave(fpath(paste0(stem, ".png")), p, width = 18, height = 11, units = "in", dpi = 300)
ggsave(fpath(paste0(stem, ".pdf")), p, width = 18, height = 11, units = "in", device = cairo_pdf)
fwrite(when, fpath(paste0(stem, "_panel_a_data.csv")))
fwrite(first$data, fpath(paste0(stem, "_panel_b_data.csv")))
fwrite(second$data, fpath(paste0(stem, "_panel_c_data.csv")))
fwrite(net_check, fpath(paste0(stem, "_diamond_checks.csv")))
caveat_lines <- c(
  "SSP2 includes 12 Swiss cities with within-age-band mortality-shape warnings. The maximum decline is 11.09%."
)
if (exploratory_style) {
  caveat_lines <- c(
    caveat_lines,
    "Exploratory style enabled: SSP1 and SSP2 bars use 48% opacity; SSP3 bars use full opacity."
  )
}
writeLines(caveat_lines, fpath(paste0(stem, "_caveats.txt")))
message("Saved Figure 2 for ", geo_id, " under ", figure_root)
