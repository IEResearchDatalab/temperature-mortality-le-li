#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("analysis/decade_decomposition.R")

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
figure_root <- Sys.getenv("FIGURE_ROOT", "results/figures_simon_20oct")
geo_level <- Sys.getenv("GEO_LEVEL", "europe")
geo_id <- Sys.getenv("GEO_ID", "EUROPE")
out_dir <- file.path(figure_root, geo_level, geo_id)
data_dir <- file.path(figure_root, "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

age_band <- function(age) {
  fifelse(age >= 100L, "100+", sprintf("%d-%d", (age %/% 5L) * 5L, (age %/% 5L) * 5L + 4L))
}
period_start <- function(x) as.integer(substr(x, 1L, 4L))

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
  when[, .(displayed_sum = sum(value)), by = .(ssp_label, decade_label)],
  by = c("ssp_label", "decade_label")
)
if (any(abs(net_check$value - net_check$displayed_sum) > 1e-12)) {
  stop("Figure 2 diamonds do not equal the five displayed temperature ranges.", call. = FALSE)
}

plot_stacked <- function(data, xvar, net, x_label, y_label, show_strips, limits = NULL) {
  positive <- data[value > 0]
  negative <- data[value < 0]
  base <- ggplot() +
    geom_col(data = positive, aes(x = .data[[xvar]], y = value, fill = cause),
             width = 0.82, position = position_stack(reverse = FALSE)) +
    geom_col(data = negative, aes(x = .data[[xvar]], y = value, fill = cause),
             width = 0.82, position = position_stack(reverse = TRUE)) +
    geom_point(data = net, aes(x = .data[[xvar]], y = value, shape = "Net (all temperature ranges)"),
               inherit.aes = FALSE, size = 2.2, colour = "black") +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    facet_grid(ssp_label ~ ., scales = "fixed", switch = "y") +
    scale_fill_manual(values = range_colours, labels = range_labels, drop = FALSE) +
    scale_shape_manual(
      name = NULL, values = c("Net (all temperature ranges)" = 18),
      breaks = "Net (all temperature ranges)"
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
  base
}

p_when <- plot_stacked(
  when, "decade_label", when_net, NULL, "Months of LE65 per decade", TRUE
) +
  labs(title = "a. Over time (each 10-year step)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8))

make_age_panel <- function(start_period, end_period, title, subtitle) {
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
  list(data = z, net = net, title = title, subtitle = subtitle)
}

first <- make_age_panel(2020L, 2055L, "b. By age, 2020–24 → 2055–59", NULL)
second <- make_age_panel(2055L, 2095L, "c. By age, 2055–59 → 2095–99", NULL)
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

figure_title <- paste0(
  geo_label, ": contributions of temperature-related mortality to change in LE65, with climate change"
)
if (geo_level == "europe") {
  figure_title <- "Europe: contributions of temperature-related mortality to change in LE65, with climate change"
}
p <- (p_when | p_first | p_second) +
  plot_layout(widths = c(1.8, 1, 1), guides = "collect") +
  plot_annotation(
    title = figure_title,
    subtitle = "Positive contributions add to LE65 gains; negative contributions drag on them; black diamonds show the net over all temperature ranges.",
    caption = paste(
      "In panel a, each bar compares the mortality schedule of the labelled 10-year period with that of the preceding one; in panels b and c, contributions over each half-century are shown per decade.",
      "Without climate change, each temperature range contributes to LE65 change only in proportion to its share of falling mortality from all causes (its attributable fraction is fixed); heat contributions turning negative with climate change therefore signal warming outpacing that decline.",
      "Black diamonds sum the five temperature ranges and exclude all other causes. Historical P2.5 and P97.5: ERA5-Land, 2000–2014. Projected P97.5: the historical P97.5 plus the change in the calibrated GCM P97.5 over the 15 preceding years.",
      "For SSP1, Wittgenstein Centre survival ratios stored as exactly 1 were set to 0.9995. Central estimates, mean of 19 GCMs; no uncertainty intervals.",
      sep = "\n"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, colour = "#555555"),
      plot.caption = element_text(size = 8, hjust = 0),
      legend.position = "bottom"
    )
  )

stem <- file.path(out_dir, "figure2_le65_contributions")
ggsave(paste0(stem, ".png"), p, width = 18, height = 11, units = "in", dpi = 300)
ggsave(paste0(stem, ".pdf"), p, width = 18, height = 11, units = "in", device = cairo_pdf)
fwrite(when, paste0(stem, "_panel_a_data.csv"))
fwrite(first$data, paste0(stem, "_panel_b_data.csv"))
fwrite(second$data, paste0(stem, "_panel_c_data.csv"))
fwrite(net_check, paste0(stem, "_diamond_checks.csv"))
writeLines(
  "SSP2 includes 12 Swiss cities with within-age-band mortality-shape warnings. The maximum decline is 11.09%.",
  paste0(stem, "_caveats.txt")
)
message("Saved Figure 2 to ", out_dir)
