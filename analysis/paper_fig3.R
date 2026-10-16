#!/usr/bin/env Rscript

################################################################################
# Paper Figure 3: cold/heat contributions to LI65+ and the LE/LI plane.
# Presentation only: aggregates completed Part 04 Horiuchi contributions.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
figure_root <- Sys.getenv("FIGURE_ROOT", "results/paper_figures/20261013")
geo_level <- Sys.getenv("GEO_LEVEL", "europe")
geo_id <- Sys.getenv("GEO_ID", "EUROPE")
out_dir <- file.path(figure_root, geo_level, geo_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

range_levels <- c("ExtrCold", "ModCold", "ModHeat", "NormExtrHeat", "ExtrHeat")
cold_causes <- c("ExtrCold", "ModCold")
heat_causes <- c("ModHeat", "NormExtrHeat", "ExtrHeat")
group_colours <- c(cold = "#24558A", heat = "#AD1A2B")
ssp_labels <- c(
  `1` = "SSP1-2.6",
  `2` = "SSP2-4.5",
  `3` = "SSP3-7.0"
)
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
  if (!nrow(x)) stop("No with-CC ensemble decomposition in ", path, call. = FALSE)
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
d[, ssp_label := factor(ssp_label, levels = ssp_labels[c("1", "2", "3")])]
geo_label <- as.character(unique(d[["label"]]))
if (length(geo_label) != 1L) stop("Geography label is not unique.", call. = FALSE)

intervals <- data.table(
  interval = c("2020-24 → 2055-59", "2055-59 → 2095-99"),
  start = c(2020L, 2055L), end = c(2055L, 2095L)
)
bars <- rbindlist(lapply(seq_len(nrow(intervals)), function(i) {
  yrs <- intervals$end[i] - intervals$start[i]
  d[
    from_start >= intervals$start[i] & to_start <= intervals$end[i] &
      cause_group %chin% c("cold", "heat"),
    .(value = sum(li_contribution) * 365 * 10 / yrs),
    by = .(ssp_label, age_band, cause_group)
  ][, interval := intervals$interval[i]]
}))
bars[, `:=`(
  age_band = factor(age_band, levels = age_levels),
  age_index = match(age_band, age_levels),
  interval = factor(interval, levels = intervals$interval),
  cause_group = factor(cause_group, levels = c("cold", "heat"))
)]

crossovers <- bars[cause_group == "cold", {
  z <- value[order(age_index)]
  idx <- age_index[order(age_index)]
  first_nonnegative <- idx[which(z >= 0)[1]]
  .(crossover = if (length(first_nonnegative)) first_nonnegative else Inf)
}, by = .(ssp_label, interval)]
net_labels <- bars[, .(
  label = sprintf("Net, all ages: cold %+.1f   heat %+.1f",
                  sum(value[cause_group == "cold"]), sum(value[cause_group == "heat"])),
  x = max(age_index),
  y = min(c(0, value)) - 0.08 * diff(range(c(0, value)))
), by = .(ssp_label, interval)]

p_bars <- ggplot(bars, aes(age_index, value, fill = cause_group)) +
  geom_rect(
    data = crossovers[is.finite(crossover)],
    aes(xmin = crossover - 0.5, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#F0F0F0", colour = NA
  ) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_text(data = net_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.6, colour = "#444444") +
  facet_grid(ssp_label ~ interval, scales = "free_y") +
  scale_x_continuous(breaks = seq_along(age_levels), labels = age_levels) +
  scale_fill_manual(values = group_colours, labels = c(cold = "Cold", heat = "Heat")) +
  labs(x = "Age group", y = "Days of LI65+ per decade", fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"), legend.position = "bottom"
  )

# Cumulative LE/LI plane since 2020-24.
plane_steps <- d[cause_group %chin% c("cold", "heat"), .(
  le = sum(le_contribution), li = sum(li_contribution)
), by = .(ssp, ssp_label, period_to, to_start, cause_group)]
setorder(plane_steps, ssp, cause_group, to_start)
plane_steps[, `:=`(
  x = cumsum(le) * 12,
  y = cumsum(li) * 365
), by = .(ssp, cause_group)]
plane_steps[, cause_group := factor(cause_group, levels = c("cold", "heat"))]
midpoints <- plane_steps[to_start == 2055L]
endpoints <- plane_steps[to_start == 2095L]
endsegments <- merge(
  plane_steps[to_start == 2090L, .(ssp, cause_group, x, y)],
  endpoints[, .(ssp, cause_group, xend = x, yend = y, ssp_label)],
  by = c("ssp", "cause_group")
)

rest_end <- d[cause_group == "rest", .(
  x = sum(le_contribution) * 12,
  y = sum(li_contribution) * 365
), by = ssp]
rest_vector <- rest_end[, .(x = mean(x), y = mean(y))]
xy_scale <- max(abs(c(plane_steps$x, plane_steps$y)), na.rm = TRUE)
rest_norm <- sqrt(rest_vector$x^2 + rest_vector$y^2)
rest_vector[, `:=`(x = x / rest_norm * xy_scale * 0.70,
                   y = y / rest_norm * xy_scale * 0.70)]

x_rng <- range(c(0, plane_steps$x, rest_vector$x), na.rm = TRUE)
y_rng <- range(c(0, plane_steps$y, rest_vector$y), na.rm = TRUE)
label_offset <- 0.050 * diff(x_rng)
label_vertical_offset <- 0.025 * diff(y_rng)
endpoints[, `:=`(
  label_x = x + fifelse(x < 0, -label_offset, label_offset),
  label_y = y + fifelse(ssp == 1L, label_vertical_offset,
                       fifelse(ssp == 3L, -label_vertical_offset, 0)),
  label_hjust = fifelse(x < 0, 1, 0),
  label_vjust = 0.5
)]
p_plane <- ggplot(plane_steps, aes(x, y, colour = cause_group,
                                  linetype = ssp_label,
                                  group = interaction(ssp, cause_group))) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FBEDEE", colour = NA) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#EEF3F8", colour = NA) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_path(linewidth = 0.75) +
  geom_point(data = midpoints, size = 2) +
  geom_segment(
    data = endsegments,
    aes(x = x, y = y, xend = xend, yend = yend,
        colour = cause_group, linetype = ssp_label),
    inherit.aes = FALSE, linewidth = 0.8,
    arrow = arrow(length = unit(0.13, "inches"), type = "closed")
  ) +
  geom_text(data = if (length(ssp_set) > 1L) endpoints else endpoints[0],
            aes(x = label_x, y = label_y, label = paste0("SSP", ssp),
                hjust = label_hjust, vjust = label_vjust),
            size = 2.6, show.legend = FALSE) +
  geom_segment(data = rest_vector,
               aes(x = 0, y = 0, xend = x, yend = y), inherit.aes = FALSE,
               colour = "#888888", linewidth = 0.6, linetype = "dashed",
               arrow = arrow(length = unit(0.11, "inches"))) +
  annotate("text", x = rest_vector$x, y = rest_vector$y,
           label = "direction of social progress\n(all other causes)",
           colour = "#777777", size = 2.6,
           hjust = ifelse(rest_vector$x >= 0, 1, 0)) +
  annotate("text", x = x_rng[1] + 0.18 * diff(x_rng), y = y_rng[1],
           label = "Shorter, more\nequal lives",
           hjust = 0, vjust = 0, size = 2.7, fontface = "italic", colour = "#9B1C31") +
  annotate("text", x = x_rng[2], y = y_rng[1], label = "Longer, more\nequal lives",
           hjust = 1, vjust = 0, size = 2.7, fontface = "italic", colour = "#24558A") +
  scale_colour_manual(values = group_colours, labels = c(cold = "Cold", heat = "Heat")) +
  scale_linetype_manual(
    values = setNames(c("solid", "dashed", "dotdash"), ssp_labels[c("1", "2", "3")]),
    breaks = ssp_labels[as.character(ssp_set)],
    labels = setNames(c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0"),
                      ssp_labels[c("1", "2", "3")])[ssp_labels[as.character(ssp_set)]]
  ) +
  labs(
    title = "Where cold and heat take older people",
    subtitle = "Dot: 2055-59; arrowhead: 2095-99",
    x = "Contribution to change in LE65 since 2020-24 (months)",
    y = "Contribution to change in LI65+ since 2020-24 (days)",
    colour = NULL, linetype = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))

p <- (p_bars | p_plane) +
  plot_layout(widths = c(2.15, 1.25), guides = "collect") +
  plot_annotation(
    title = paste0(geo_label, ": contributions of cold and heat to change in lifespan inequality above 65"),
    subtitle = "Bars show age-specific contributions; the plane shows cumulative cold, heat and social-progress directions",
    caption = paste(
      "Bars are contributions between 5-year-period mean mortality schedules, expressed as days of LI65+ per decade; shading begins at the data-derived cold crossover age.",
      "Static ERA5-Land baseline 2000-2014; moving P97.5 uses the preceding 15 years and is anchored as observed baseline plus calibrated GCM change.",
      "For SSP1, WCDE survival ratios rounded to exactly one were replaced with 0.9995.",
      "Central estimate, no uncertainty intervals."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, colour = "#555555"),
      plot.caption = element_text(size = 8, hjust = 0),
      legend.position = "bottom"
    )
  )

stem <- file.path(out_dir, "figure3_li65_contributions")
ggsave(paste0(stem, ".png"), p, width = 16, height = 11, units = "in", dpi = 300)
ggsave(paste0(stem, ".pdf"), p, width = 16, height = 11, units = "in", device = cairo_pdf)
fwrite(bars, paste0(stem, "_bars_data.csv"))
fwrite(plane_steps, paste0(stem, "_plane_data.csv"))
message("Saved Figure 3 to ", out_dir)
