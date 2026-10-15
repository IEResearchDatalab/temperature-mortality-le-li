#!/usr/bin/env Rscript

################################################################################
# Paper Figure 2: contributions to changes in LE65.
# Presentation only: sums completed Part 04 Horiuchi contributions over causes,
# ages and consecutive 5-year-period steps. It does not rerun decomposition.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
figure_root <- Sys.getenv("FIGURE_ROOT", "results/paper_figures/20261013")
geo_level <- Sys.getenv("GEO_LEVEL", "europe")
geo_id <- Sys.getenv("GEO_ID", "EUROPE")
out_dir <- file.path(figure_root, geo_level, geo_id)
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
ssp_labels <- c(
  `1` = "SSP1-2.6 — survival ratios patched (WCDE rounding)",
  `2` = "SSP2-4.5",
  `3` = "SSP3-7.0"
)
ssp_order <- c("1", "3", "2")
age_levels <- c("65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100+")
ssp_set <- as.integer(strsplit(Sys.getenv("SSP_SET", "1,2,3"), ",", fixed = TRUE)[[1]])
if (!length(ssp_set) || anyNA(ssp_set) || any(!ssp_set %in% 1:3)) {
  stop("SSP_SET must be a comma-separated subset of 1,2,3.", call. = FALSE)
}

age_band <- function(age) {
  fifelse(age >= 100L, "100+", sprintf("%d-%d", (age %/% 5L) * 5L, (age %/% 5L) * 5L + 4L))
}
period_start <- function(x) as.integer(substr(x, 1L, 4L))

read_ssp <- function(ssp) {
  path <- file.path(batch_root, "geographies", paste0("ssp", ssp), geo_level,
                    geo_id, "04_within_branch_period_decomposition.csv")
  if (!file.exists(path)) stop("Missing within-branch decomposition: ", path, call. = FALSE)
  x <- fread(path)[gcm == "ENSEMBLE" & branch == "with_cc" & cause %chin% range_levels]
  if (!nrow(x)) stop("No with-CC ensemble decomposition rows in ", path, call. = FALSE)
  if (length(setdiff(range_levels, unique(x$cause)))) {
    stop("Incomplete cause domain in ", path, call. = FALSE)
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

d <- rbindlist(lapply(ssp_set, read_ssp))
d[, ssp_label := factor(ssp_label, levels = ssp_labels[ssp_order])]
geo_label <- unique(d$label)
if (length(geo_label) != 1L) stop("Geography label is not unique.", call. = FALSE)

# Panel a: each 5-year step, expressed as months per decade.
when <- d[, .(value = sum(le_contribution) * 12 * 2),
          by = .(ssp_label, period_to, cause)]
when[, cause := factor(cause, levels = range_levels)]
when_net <- when[, .(value = sum(value)), by = .(ssp_label, period_to)]

p_when <- ggplot(when, aes(period_to, value, fill = cause)) +
  geom_col(width = 0.82) +
  geom_point(data = when_net, aes(period_to, value), inherit.aes = FALSE,
             shape = 18, size = 2.2, colour = "black") +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_grid(ssp_label ~ ., scales = "fixed") +
  scale_fill_manual(values = range_colours, labels = range_labels, drop = FALSE) +
  labs(title = "a  When", subtitle = "Every 5-year step", x = NULL,
       y = "Months of LE65 per decade") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
    panel.grid.minor = element_blank(), strip.text.y = element_text(face = "bold"),
    legend.position = "none", plot.title = element_text(face = "bold")
  )

make_age_panel <- function(start_period, end_period, title, subtitle) {
  years <- end_period - start_period
  z <- d[from_start >= start_period & to_start <= end_period,
         .(value = sum(le_contribution) * 12 * 10 / years),
         by = .(ssp_label, age_band, cause)]
  z[, `:=`(
    age_band = factor(age_band, levels = age_levels),
    cause = factor(cause, levels = range_levels)
  )]
  net <- z[, .(value = sum(value)), by = .(ssp_label, age_band)]
  list(
    plot = ggplot(z, aes(age_band, value, fill = cause)) +
      geom_col(width = 0.82) +
      geom_point(data = net, aes(age_band, value), inherit.aes = FALSE,
                 shape = 18, size = 2.2, colour = "black") +
      geom_hline(yintercept = 0, linewidth = 0.35) +
      facet_grid(ssp_label ~ ., scales = "fixed") +
      scale_fill_manual(values = range_colours, labels = range_labels, drop = FALSE) +
      labs(title = title, subtitle = subtitle, x = "Age group",
           y = "Months of LE65 per decade", fill = NULL) +
      theme_minimal(base_size = 9) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
        panel.grid.minor = element_blank(), strip.text.y = element_blank(),
        legend.position = "bottom", plot.title = element_text(face = "bold")
      ),
    data = z,
    net = net
  )
}

first <- make_age_panel(2020L, 2055L, "b  For whom: first half", "2020-24 → 2055-59")
second <- make_age_panel(2055L, 2095L, "c  For whom: second half", "2055-59 → 2095-99")

figure_title <- paste0(geo_label, ": temperature-related contributions to change in LE65 with climate change")
if (nchar(figure_title) > 92L) {
  figure_title <- sub(": ", ":\n", figure_title, fixed = TRUE)
}

p <- (p_when | first$plot | second$plot) +
  plot_layout(widths = c(2.25, 1, 1), guides = "collect") +
  plot_annotation(
    title = figure_title,
    subtitle = "Positive contributions add to LE65 gains; negative contributions erode them; black diamonds are net temperature contributions",
    caption = paste(
      "In panel a, each bar compares the labelled 5-year-period mean mortality schedule with the immediately preceding 5-year period; contributions are scaled to months of LE65 per decade.",
      "Black diamonds sum the five displayed temperature categories and exclude all other causes (social progress).",
      "Static ERA5-Land baseline 2000-2014; moving P97.5 uses the preceding 15 years and is anchored as observed baseline plus calibrated GCM change.",
      "Central estimate, no uncertainty intervals.",
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
ggsave(paste0(stem, ".png"), p, width = 16, height = 11, units = "in", dpi = 300)
ggsave(paste0(stem, ".pdf"), p, width = 16, height = 11, units = "in", device = cairo_pdf)
fwrite(when, paste0(stem, "_panel_a_data.csv"))
fwrite(first$data, paste0(stem, "_panel_b_data.csv"))
fwrite(second$data, paste0(stem, "_panel_c_data.csv"))
message("Saved Figure 2 to ", out_dir)
