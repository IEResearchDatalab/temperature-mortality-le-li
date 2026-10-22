#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
})

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
figure_root <- Sys.getenv("FIGURE_ROOT", "results/figures_simon_20oct")
geo_level <- Sys.getenv("GEO_LEVEL", "europe")
geo_id <- Sys.getenv("GEO_ID", "EUROPE")
out_dir <- file.path(figure_root, geo_level, geo_id)
data_dir <- file.path(figure_root, "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

range_levels <- c("ExtrCold", "ModCold", "ModHeat", "NormExtrHeat", "ExtrHeat")
range_labels <- c(
  ExtrCold = "Extreme cold", ModCold = "Moderate cold",
  ModHeat = "Moderate heat", NormExtrHeat = "Normal extreme heat",
  ExtrHeat = "Extreme heat"
)
cold_causes <- c("ExtrCold", "ModCold")
heat_causes <- c("ModHeat", "NormExtrHeat", "ExtrHeat")
group_colours <- c(cold = "#24558A", heat = "#AD1A2B")
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
cause_group <- function(x) fifelse(x %chin% cold_causes, "cold",
                                   fifelse(x %chin% heat_causes, "heat", "rest"))

read_ssp <- function(ssp) {
  path <- file.path(batch_root, "geographies", paste0("ssp", ssp), geo_level,
                    geo_id, "04_within_branch_period_decomposition.csv")
  if (!file.exists(path)) stop("Missing within-branch decomposition: ", path, call. = FALSE)
  x <- fread(path)[gcm == "ENSEMBLE" & branch == "with_cc"]
  if (!nrow(x) || length(setdiff(range_levels, unique(x$cause)))) {
    stop("Incomplete with-CC ensemble data in ", path, call. = FALSE)
  }
  x[, `:=`(
    ssp = ssp,
    ssp_label = ssp_labels[as.character(ssp)],
    from_start = period_start(period_from),
    to_start = period_start(period_to),
    age_band = age_band(age),
    cause_group = cause_group(cause)
  )]
  x
}

d <- rbindlist(lapply(ssp_set, read_ssp))
d[, ssp_label := factor(ssp_label, levels = ssp_labels[ssp_order])]
geo_labels <- unique(d$label)
if (length(geo_labels) != 1L) stop("Geography label is not unique.", call. = FALSE)
geo_label <- geo_labels[[1]]

intervals <- data.table(
  interval = c("a. By age, 2020–24 → 2055–59", "b. By age, 2055–59 → 2095–99"),
  start = c(2020L, 2055L), end = c(2055L, 2095L)
)

make_bars <- function(causes, collapse = TRUE) {
  source <- d[cause %chin% causes]
  if (collapse) source[, display_cause := cause_group(cause)]
  else source[, display_cause := cause]
  bars <- rbindlist(lapply(seq_len(nrow(intervals)), function(i) {
    yrs <- intervals$end[i] - intervals$start[i]
    source[
      from_start >= intervals$start[i] & to_start <= intervals$end[i],
      .(value = sum(li_contribution) * 365 * 10 / yrs),
      by = .(ssp, ssp_label, age_band, display_cause)
    ][, interval := intervals$interval[i]]
  }))
  bars[, `:=`(
    age_band = factor(age_band, levels = age_levels),
    age_index = match(age_band, age_levels),
    interval = factor(interval, levels = intervals$interval),
    display_cause = if (collapse) factor(display_cause, levels = c("cold", "heat"))
    else factor(display_cause, levels = range_levels)
  )]
  bars
}

bars <- make_bars(c(cold_causes, heat_causes))
crossovers <- bars[display_cause == "cold", {
  ordered <- .SD[order(age_index)]
  first_nonnegative <- ordered[which(value >= 0)[1], age_index]
  .(crossover = if (length(first_nonnegative)) first_nonnegative else Inf)
}, by = .(ssp, ssp_label, interval)]
crossovers[, crossover_label := ifelse(
  is.finite(crossover),
  paste0("Crossover ≈ ", c("65", "70", "75", "80", "85", "90", "95", "100+")[crossover]),
  "No non-negative band"
)]
crossovers[, crossover_side_label := "← below crossover age | above →"]
net_labels <- bars[, .(
  cold_net = sum(value[display_cause == "cold"]),
  heat_net = sum(value[display_cause == "heat"])
), by = .(ssp, ssp_label, interval)]
net_labels[, label := paste0(
  "Net, all ages (days per decade): cold ",
  ifelse(cold_net < 0, "−", "+"), sprintf("%.1f", abs(cold_net)),
  ", heat ", ifelse(heat_net < 0, "−", "+"), sprintf("%.1f", abs(heat_net))
)]
net_labels <- merge(net_labels, crossovers, by = c("ssp", "ssp_label", "interval"))

bar_extent <- bars[, .(
  positive = sum(pmax(value, 0)),
  negative = sum(pmin(value, 0))
), by = .(ssp, ssp_label, interval, age_band)]
bar_panel_extent <- bar_extent[, .(
  positive = max(positive),
  negative = min(negative)
), by = .(ssp, ssp_label, interval)]
net_labels <- merge(net_labels, bar_panel_extent, by = c("ssp", "ssp_label", "interval"))
bar_limits <- range(c(0, bars$value, bar_extent$positive, bar_extent$negative))
bar_padding <- diff(bar_limits) * 0.24
bar_limits <- bar_limits + c(-bar_padding, bar_padding)
net_labels[, `:=`(x = 0.55, y = bar_limits[1] + diff(bar_limits) * 0.025)]
net_labels[, `:=`(
  crossover_x = fifelse(is.finite(crossover), pmin(crossover + 0.1, 7.9), 8),
  crossover_y = positive + diff(bar_limits) * 0.035,
  crossover_side_y = positive + diff(bar_limits) * 0.12
)]
net_label_plot <- copy(net_labels)
net_label_plot[, `:=`(
  label = paste0(
    "Net, all ages\n(days per decade)\ncold ",
    ifelse(cold_net < 0, "−", "+"), sprintf("%.1f", abs(cold_net)),
    ", heat ", ifelse(heat_net < 0, "−", "+"), sprintf("%.1f", abs(heat_net))
  ),
  x = fifelse(is.finite(crossover) & (crossover - 0.5) <= 4.5, 8.35, 0.65),
  hjust = fifelse(is.finite(crossover) & (crossover - 0.5) <= 4.5, 1, 0),
  y = bar_limits[1] + diff(bar_limits) * 0.05
)]
main_axis_checks <- bar_panel_extent[, .(
  data_min = min(negative),
  data_max = max(positive)
), by = .(ssp, ssp_label, interval)]
main_axis_checks[, `:=`(axis_min = bar_limits[1], axis_max = bar_limits[2])]
main_axis_checks[, status := ifelse(
  data_min >= axis_min & data_max <= axis_max,
  "PASS", "FAIL"
)]
if (any(main_axis_checks$status != "PASS") ||
    any(net_labels$y < bar_limits[1] | net_labels$y > bar_limits[2]) ||
    any(net_labels$crossover_side_y > bar_limits[2], na.rm = TRUE)) {
  stop("Figure 3 age bars or labels exceed the shared y-axis limits.", call. = FALSE)
}
message(
  "Main age-bar shared y-limits for ", geo_level, "/", geo_id, ": ",
  paste(format(bar_limits, digits = 6), collapse = " to "),
  "; data extents: ",
  paste(format(range(c(bar_panel_extent$negative, bar_panel_extent$positive)),
               digits = 6), collapse = " to ")
)
print(main_axis_checks)

p_bars <- ggplot(bars, aes(age_index, value, fill = display_cause)) +
  geom_rect(
    data = crossovers[is.finite(crossover)],
    aes(xmin = crossover - 0.5, xmax = 8.5, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#F0F0F0", colour = NA
  ) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_vline(
    data = crossovers[is.finite(crossover)],
    aes(xintercept = crossover - 0.5), inherit.aes = FALSE,
    colour = "#777777", linewidth = 0.35, linetype = "dotted"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_text(data = net_label_plot, aes(x = x, y = y, label = label, hjust = hjust),
            inherit.aes = FALSE, vjust = 0, size = 2.5, colour = "#444444") +
  geom_text(data = net_labels[is.finite(crossover)],
            aes(x = crossover - 0.5, y = crossover_side_y,
                label = crossover_side_label),
            inherit.aes = FALSE, hjust = 0.5, vjust = 0, size = 2.6,
            colour = "#666666") +
  geom_text(data = net_labels[is.finite(crossover)],
            aes(x = crossover_x, y = crossover_y, label = crossover_label),
            inherit.aes = FALSE, hjust = 0, vjust = 0, size = 2.5, colour = "#444444") +
  facet_grid(ssp_label ~ interval, scales = "fixed", switch = "y") +
  scale_x_continuous(breaks = seq_along(age_levels), labels = age_levels) +
  scale_fill_manual(
    values = group_colours,
    labels = c(cold = "Cold (extreme + moderate)",
               heat = "Heat (moderate + normal extreme + extreme)")
  ) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 6),
                     labels = scales::label_number(style_negative = "minus")) +
  coord_cartesian(ylim = bar_limits) +
  labs(x = "Age group", y = "Days of LI65+ per decade", fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(face = "bold", angle = 0, size = 10),
    strip.text.x = element_text(face = "bold", size = 9),
    legend.position = "none"
  )

plane_source <- d[cause_group %chin% c("cold", "heat"), .(
  le = sum(le_contribution), li = sum(li_contribution)
), by = .(ssp, ssp_label, period_to, to_start, cause_group)]
setorder(plane_source, ssp, cause_group, to_start)
plane_source[, `:=`(
  le_months = cumsum(le) * 12,
  li_days = cumsum(li) * 365
), by = .(ssp, cause_group)]

expected_endpoints <- data.table(
  ssp = rep(1:3, each = 4L),
  cause_group = rep(c("cold", "heat"), 6L),
  to_start = rep(c(2055L, 2055L, 2095L, 2095L), 3L),
  expected_le = c(5.454, 0.922, 14.722, 4.038,
                  3.767, -0.131, 8.715, 1.161,
                  2.252, -0.856, 4.608, -2.401),
  expected_li = c(-23.486, -8.893, -31.786, -3.029,
                  -8.248, -6.282, -43.913, -19.990,
                  1.561, -8.531, -3.178, -23.129)
)
endpoint_check <- plane_source[to_start %in% c(2055L, 2095L)]
if (geo_level == "europe" && geo_id == "EUROPE") {
  endpoint_check <- merge(
    endpoint_check, expected_endpoints,
    by = c("ssp", "cause_group", "to_start"), all = TRUE
  )
  endpoint_check[, `:=`(
    le_error = le_months - expected_le,
    li_error = li_days - expected_li
  )]
  if (nrow(endpoint_check) != nrow(expected_endpoints) ||
      anyNA(endpoint_check$le_error) || anyNA(endpoint_check$li_error) ||
      any(abs(endpoint_check$le_error) > 0.001) ||
      any(abs(endpoint_check$li_error) > 0.001)) {
    stop("Plane endpoint check does not match the audited source values.", call. = FALSE)
  }
  endpoint_check[, audit_status := "PASS"]
} else {
  endpoint_check[, `:=`(
    expected_le = NA_real_, expected_li = NA_real_,
    le_error = NA_real_, li_error = NA_real_,
    audit_status = "NOT_APPLICABLE"
  )]
}

origin <- data.table(
  ssp = rep(ssp_set, each = 2L),
  ssp_label = factor(rep(ssp_labels[as.character(ssp_set)], each = 2L),
                     levels = ssp_labels[ssp_order]),
  cause_group = rep(c("cold", "heat"), length(ssp_set)),
  li_days = 0, le_months = 0, point_order = 0L,
  period_label = "2020–24"
)
mid_end <- plane_source[to_start %in% c(2055L, 2095L)][
  , .(ssp, ssp_label, cause_group, li_days, le_months,
      point_order = fifelse(to_start == 2055L, 1L, 2L),
      period_label = fifelse(to_start == 2055L, "2055–59", "2095–99"))
]
plane_main <- rbind(origin, mid_end)
setorder(plane_main, ssp, cause_group, point_order)
arrow_segments <- merge(
  plane_main[point_order == 1L, .(
    ssp, ssp_label, cause_group, x = li_days, y = le_months
  )],
  plane_main[point_order == 2L, .(
    ssp, cause_group, xend = li_days, yend = le_months
  )],
  by = c("ssp", "cause_group")
)
endpoint_labels <- copy(plane_main[point_order == 2L])
endpoint_x_span <- diff(range(c(0, endpoint_labels$li_days)))
endpoint_y_span <- diff(range(c(0, endpoint_labels$le_months)))
endpoint_x_pad <- max(1.5, endpoint_x_span * 0.035)
endpoint_y_pad <- max(0.6, endpoint_y_span * 0.04)
endpoint_x_bounds <- range(c(0, endpoint_labels$li_days))
endpoint_labels[, label_side := fifelse(cause_group == "cold", 1, -1)]
endpoint_labels[li_days <= endpoint_x_bounds[1] + 0.12 * endpoint_x_span, label_side := 1]
endpoint_labels[li_days >= endpoint_x_bounds[2] - 0.12 * endpoint_x_span, label_side := -1]
endpoint_labels[, `:=`(
  label = paste(ssp_labels[as.character(ssp)], cause_group),
  label_x = li_days + label_side * endpoint_x_pad,
  label_y = le_months + fifelse(le_months >= 0, endpoint_y_pad, -endpoint_y_pad),
  label_hjust = fifelse(label_side > 0, 0, 1)
)]

x_range <- range(c(0, plane_main$li_days, endpoint_labels$label_x))
y_range <- range(c(0, plane_main$le_months, endpoint_labels$label_y))
x_span <- diff(x_range)
y_span <- diff(y_range)
plane_heading <- if (geo_level == "europe") {
  "c. Where cold and heat take older Europeans"
} else {
  paste0("c. Where cold and heat take older people in ", geo_label)
}
p_plane <- ggplot(plane_main, aes(li_days, le_months,
                                  colour = cause_group,
                                  group = interaction(ssp, cause_group))) +
  annotate("rect", xmin = x_range[1], xmax = 0,
           ymin = 0, ymax = y_range[2], fill = "#EEF3F8", colour = NA) +
  annotate("rect", xmin = x_range[1], xmax = 0,
           ymin = y_range[1], ymax = 0, fill = "#FBEDEE", colour = NA) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_path(linewidth = 0.75) +
  geom_point(data = plane_main[point_order == 1L], size = 2) +
  geom_segment(
    data = arrow_segments,
    aes(x = x, y = y, xend = xend, yend = yend, colour = cause_group),
    inherit.aes = FALSE, linewidth = 0.8, show.legend = FALSE,
    arrow = arrow(length = unit(0.13, "inches"), type = "closed")
  ) +
  geom_segment(
    data = endpoint_labels,
    aes(x = li_days, y = le_months, xend = label_x, yend = label_y,
        colour = cause_group),
    inherit.aes = FALSE, linewidth = 0.2, alpha = 0.55, show.legend = FALSE
  ) +
  geom_text(
    data = endpoint_labels,
    aes(
      x = label_x, y = label_y, label = label, colour = cause_group,
      hjust = label_hjust
    ),
    inherit.aes = FALSE, size = 2.7, show.legend = FALSE
  ) +
  annotate("text", x = x_range[1] + 0.08 * x_span,
           y = y_range[2] - 0.07 * y_span,
           label = "Longer, more equal lives", hjust = 0,
           size = 2.7, fontface = "italic", colour = "#24558A") +
  annotate("text", x = x_range[1] + 0.08 * x_span,
           y = y_range[1] + 0.08 * y_span,
           label = "Shorter, more equal lives", hjust = 0,
           size = 2.7, fontface = "italic", colour = "#9B1C31") +
  scale_colour_manual(
    values = group_colours,
    labels = c(cold = "Cold (extreme + moderate)",
               heat = "Heat (moderate + normal extreme + extreme)")
  ) +
  scale_x_continuous(
    labels = scales::label_number(style_negative = "minus"),
    expand = expansion(mult = c(0.08, 0.15))
  ) +
  scale_y_continuous(
    labels = scales::label_number(style_negative = "minus"),
    expand = expansion(mult = c(0.08, 0.1))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = plane_heading,
    subtitle = "Cumulative since 2020–24; dot: 2055–59; arrowhead: 2095–99",
    x = "Contribution to change in LI65+ since 2020–24 (days)",
    y = "Contribution to change in LE65 since 2020–24 (months)",
    colour = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 10, margin = margin(b = 1)),
    plot.subtitle = element_text(size = 7, margin = margin(b = 2)),
    plot.margin = margin(0, 2, 0, 2)
  ) +
  guides(colour = guide_legend(override.aes = list(
    linetype = 0, shape = 22, fill = unname(group_colours), size = 4
  )))

figure_title <- paste0(geo_label, ": contributions of cold and heat to change in LI65+, with climate change")
if (geo_level == "europe") {
  figure_title <- "Europe: contributions of cold and heat to change in LI65+, with climate change"
}
caption_paragraphs <- c(
  "In panels a and b, contributions over each half-century are shown per decade; grey shading marks ages above the crossover age, where cold's contributions change sign. In panel c, each line shows how changes in cold (or heat) mortality alone have moved LE65 and LI65+ since 2020–24; the total change, which also includes all other causes, is given in the text.",
  "Without climate change, each temperature range's contributions follow the age pattern of falling mortality from all causes (its attributable fraction is fixed); heat contributions that reverse this pattern signal warming outpacing that decline. A five-range version is in the Supplementary Material. Historical P2.5 and P97.5: ERA5-Land, 2000–2014. Projected P97.5: the historical P97.5 plus the change in the calibrated GCM P97.5 over the 15 preceding years.",
  "For SSP1, Wittgenstein Centre survival ratios stored as exactly 1 were set to 0.9995. Central estimates, mean of 19 GCMs; no uncertainty intervals."
)
five_range_caption <- caption_paragraphs
five_range_caption[1] <- paste(
  "Panels a and b show the same contributions as Figure 3 split into five",
  "temperature ranges."
)
step_plane_caption <- caption_paragraphs
step_plane_caption[1] <- paste(
  "The plane shows the five-year cumulative paths of changes in cold (or heat)",
  "mortality alone to LE65 and LI65+ since 2020–24. Dots mark 2055–59 and",
  "arrowheads mark 2095–99; the total change, which also includes all other",
  "causes, is given in the text."
)
figure_caption <- paste(
  vapply(caption_paragraphs, function(x) paste(strwrap(x, width = 190), collapse = "\n"),
         character(1)),
  collapse = "\n"
)
format_caption <- function(paragraphs) {
  paste(
    vapply(paragraphs, function(x) paste(strwrap(x, width = 190), collapse = "\n"),
           character(1)),
    collapse = "\n"
  )
}
five_range_figure_caption <- format_caption(five_range_caption)
step_plane_figure_caption <- format_caption(step_plane_caption)
p <- (p_bars | p_plane) +
  plot_layout(widths = c(1.65, 1.2), guides = "collect") +
  plot_annotation(
    title = figure_title,
    subtitle = "Positive contributions raise LI65+ (lifespans more unequal); negative contributions lower it. Grey shading: ages above the crossover age, where an added death lowers LI65+.",
    caption = figure_caption,
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, colour = "#555555"),
      plot.caption = element_text(size = 8, hjust = 0),
      legend.position = "bottom"
    )
  )

stem <- file.path(out_dir, "figure3_li65_contributions")
ggsave(paste0(stem, ".png"), p, width = 18, height = 11, units = "in", dpi = 300)
ggsave(paste0(stem, ".pdf"), p, width = 18, height = 11, units = "in", device = cairo_pdf)
fwrite(bars, paste0(stem, "_bars_data.csv"))
fwrite(net_labels, paste0(stem, "_net_crossover_data.csv"))
fwrite(main_axis_checks, paste0(stem, "_axis_checks.csv"))
fwrite(plane_main, paste0(stem, "_plane_data.csv"))
fwrite(endpoint_check, paste0(stem, "_plane_endpoint_checks.csv"))
writeLines(
  "SSP2 includes 12 Swiss cities with within-age-band mortality-shape warnings. The maximum decline is 11.09%.",
  paste0(stem, "_caveats.txt")
)

plane_steps <- copy(plane_source)
plane_steps[, `:=`(
  x = li_days,
  y = le_months,
  ssp_label = factor(ssp_labels[as.character(ssp)], levels = ssp_labels[ssp_order])
)]
step_x_range <- range(c(0, plane_steps$x))
step_y_range <- range(c(0, plane_steps$y))
step_ends <- plane_steps[to_start == max(to_start)]
step_previous <- plane_steps[to_start == max(to_start) - 5L]
step_arrows <- merge(
  step_previous[, .(ssp, ssp_label, cause_group, x = x, y = y)],
  step_ends[, .(ssp, ssp_label, cause_group, xend = x, yend = y)],
  by = c("ssp", "ssp_label", "cause_group")
)
step_labels <- copy(step_ends)
step_labels[, label := paste(ssp_labels[as.character(ssp)], cause_group)]
step_x_span <- diff(range(c(0, step_labels$x)))
step_y_span <- diff(range(c(0, step_labels$y)))
step_x_pad <- max(1.5, step_x_span * 0.035)
step_y_pad <- max(0.6, step_y_span * 0.04)
step_x_bounds <- range(c(0, step_labels$x))
step_labels[, label_side := fifelse(cause_group == "cold", 1, -1)]
step_labels[x <= step_x_bounds[1] + 0.12 * step_x_span, label_side := 1]
step_labels[x >= step_x_bounds[2] - 0.12 * step_x_span, label_side := -1]
step_labels[, `:=`(
  label_x = x + label_side * step_x_pad,
  label_y = y + fifelse(y >= 0, step_y_pad, -step_y_pad),
  label_hjust = fifelse(label_side > 0, 0, 1)
)]
step_x_range <- range(c(step_x_range, step_labels$label_x))
step_y_range <- range(c(step_y_range, step_labels$label_y))
step_x_range <- step_x_range + c(-0.08, 0.08) * diff(step_x_range)
step_y_range <- step_y_range + c(-0.05, 0.05) * diff(step_y_range)
step_plot <- ggplot(plane_steps, aes(x, y, colour = cause_group,
                                     group = interaction(ssp, cause_group))) +
  annotate("rect", xmin = step_x_range[1], xmax = 0, ymin = 0, ymax = step_y_range[2],
           fill = "#EEF3F8", colour = NA) +
  annotate("rect", xmin = step_x_range[1], xmax = 0, ymin = step_y_range[1], ymax = 0,
           fill = "#FBEDEE", colour = NA) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_path(linewidth = 0.65) +
  geom_point(data = plane_steps[to_start %in% c(2055L, 2095L)], size = 1.5) +
  geom_segment(
    data = step_arrows,
    aes(x = x, y = y, xend = xend, yend = yend, colour = cause_group),
    inherit.aes = FALSE, linewidth = 0.7, show.legend = FALSE,
    arrow = arrow(length = unit(0.11, "inches"), type = "closed")
  ) +
  geom_segment(
    data = step_labels,
    aes(x = x, y = y, xend = label_x, yend = label_y, colour = cause_group),
    inherit.aes = FALSE, linewidth = 0.2, alpha = 0.55, show.legend = FALSE
  ) +
  geom_text(
    data = step_labels,
    aes(
      x = label_x, y = label_y, label = label, colour = cause_group,
      hjust = label_hjust
    ),
    inherit.aes = FALSE, size = 2.6, show.legend = FALSE
  ) +
  facet_grid(ssp_label ~ ., switch = "y") +
  scale_colour_manual(
    values = group_colours,
    labels = c(cold = "Cold (extreme + moderate)",
               heat = "Heat (moderate + normal extreme + extreme)")
  ) +
  scale_x_continuous(labels = scales::label_number(style_negative = "minus")) +
  scale_y_continuous(labels = scales::label_number(style_negative = "minus")) +
  labs(
    title = paste0(geo_label, ": five-year cumulative cold and heat contributions"),
    x = "Contribution to change in LI65+ since 2020–24 (days)",
    y = "Contribution to change in LE65 since 2020–24 (months)",
    colour = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(face = "bold", angle = 0)) +
  coord_cartesian(xlim = step_x_range, ylim = step_y_range, clip = "off") +
  plot_annotation(caption = step_plane_figure_caption) &
  theme(plot.caption = element_text(size = 7, hjust = 0))
step_stem <- file.path(out_dir, "figure3_supplement_plane_steps")
ggsave(paste0(step_stem, ".png"), step_plot, width = 10, height = 8, units = "in", dpi = 300)
ggsave(paste0(step_stem, ".pdf"), step_plot, width = 10, height = 8, units = "in", device = cairo_pdf)
fwrite(plane_steps, paste0(step_stem, "_data.csv"))

five_range <- make_bars(range_levels, collapse = FALSE)
five_net <- five_range[, .(value = sum(value)), by = .(ssp, ssp_label, age_band, age_index, interval)]
five_extent <- five_range[, .(
  positive = sum(pmax(value, 0)),
  negative = sum(pmin(value, 0))
), by = .(ssp, ssp_label, interval, age_band)]
five_limits <- range(c(0, five_range$value, five_extent$positive,
                       five_extent$negative, five_net$value))
five_limits <- five_limits + c(-1, 1) * diff(five_limits) * 0.12
five_axis_checks <- five_extent[, .(
  data_min = min(negative),
  data_max = max(positive)
), by = .(ssp, ssp_label, interval)]
five_axis_checks <- merge(
  five_axis_checks,
  five_net[, .(diamond_min = min(value), diamond_max = max(value)),
           by = .(ssp, ssp_label, interval)],
  by = c("ssp", "ssp_label", "interval")
)
five_axis_checks[, `:=`(axis_min = five_limits[1], axis_max = five_limits[2])]
five_axis_checks[, status := ifelse(
  data_min >= axis_min & data_max <= axis_max &
    diamond_min >= axis_min & diamond_max <= axis_max,
  "PASS", "FAIL"
)]
if (any(five_axis_checks$status != "PASS")) {
  stop("Five-range contributions exceed the shared y-axis limits.", call. = FALSE)
}
five_plot <- ggplot(five_range, aes(age_index, value, fill = display_cause)) +
  geom_rect(
    data = crossovers[is.finite(crossover)],
    aes(xmin = crossover - 0.5, xmax = 8.5, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#F0F0F0", colour = NA
  ) +
  geom_col(
    data = five_range[value > 0],
    position = position_stack(reverse = FALSE), width = 0.82
  ) +
  geom_col(
    data = five_range[value < 0],
    position = position_stack(reverse = TRUE), width = 0.82
  ) +
  geom_point(data = five_net, aes(age_index, value), inherit.aes = FALSE,
             shape = 18, size = 1.8, colour = "black") +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_grid(ssp_label ~ interval, switch = "y") +
  scale_x_continuous(breaks = seq_along(age_levels), labels = age_levels) +
  scale_fill_manual(values = range_colours, labels = range_labels, drop = FALSE) +
  scale_y_continuous(labels = scales::label_number(style_negative = "minus")) +
  coord_cartesian(ylim = five_limits) +
  labs(
    title = paste0(geo_label, ": five-range contributions to change in LI65+"),
    x = "Age group", y = "Days of LI65+ per decade", fill = NULL
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(face = "bold", angle = 0),
    strip.text.x = element_text(face = "bold"),
    legend.position = "bottom"
  )
five_plot <- five_plot + plot_annotation(caption = five_range_figure_caption) &
  theme(plot.caption = element_text(size = 7, hjust = 0))
message(
  "Five-range shared y-limits for ", geo_level, "/", geo_id, ": ",
  paste(format(five_limits, digits = 6), collapse = " to ")
)
print(five_axis_checks)
five_stem <- file.path(out_dir, "figure3_supplement_five_ranges")
ggsave(paste0(five_stem, ".png"), five_plot, width = 15, height = 9, units = "in", dpi = 300)
ggsave(paste0(five_stem, ".pdf"), five_plot, width = 15, height = 9, units = "in", device = cairo_pdf)
fwrite(five_range, paste0(five_stem, "_bars_data.csv"))
fwrite(five_net, paste0(five_stem, "_net_data.csv"))
fwrite(five_axis_checks, paste0(five_stem, "_axis_checks.csv"))
writeLines(
  "SSP2 includes 12 Swiss cities with within-age-band mortality-shape warnings. The maximum decline is 11.09%.",
  paste0(five_stem, "_caveats.txt")
)
writeLines(
  "SSP2 includes 12 Swiss cities with within-age-band mortality-shape warnings. The maximum decline is 11.09%.",
  paste0(step_stem, "_caveats.txt")
)
message("Saved Figure 3 and supplements to ", out_dir)
