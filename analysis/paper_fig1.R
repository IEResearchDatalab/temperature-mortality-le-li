#!/usr/bin/env Rscript

################################################################################
# Paper Figure 1: absolute temperature-attributable mortality risk by age.
# Presentation only: reads completed pooled Part 03 tables and performs no
# temperature attribution, life-table or decomposition calculation.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
  `1` = "SSP1-2.6",
  `2` = "SSP2-4.5",
  `3` = "SSP3-7.0"
)
age_levels <- c("65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95+")
period_levels <- c("2055-2059", "2095-2099")
ssp_set <- as.integer(strsplit(Sys.getenv("SSP_SET", "1,2,3"), ",", fixed = TRUE)[[1]])
if (!length(ssp_set) || anyNA(ssp_set) || any(!ssp_set %in% 1:3)) {
  stop("SSP_SET must be a comma-separated subset of 1,2,3.", call. = FALSE)
}

age_band <- function(age) {
  cut(age, breaks = c(64, 69, 74, 79, 84, 89, 94, Inf),
      labels = age_levels, right = TRUE)
}

read_ssp <- function(ssp) {
  path <- file.path(batch_root, "geographies", paste0("ssp", ssp),
                    geo_level, geo_id, "03_master_table.csv")
  if (!file.exists(path)) stop("Missing pooled master table: ", path, call. = FALSE)
  x <- fread(path)
  missing_causes <- setdiff(range_levels, unique(x$range))
  if (length(missing_causes)) {
    stop("Missing Figure 1 causes for SSP", ssp, ": ",
         paste(missing_causes, collapse = ", "), call. = FALSE)
  }
  label <- unique(x$label)
  if (length(label) != 1L) stop("Geography label is not unique in ", path, call. = FALSE)
  x <- x[year %in% c(2055:2059, 2095:2099) & range %chin% range_levels]
  x[, `:=`(
    age_band = age_band(age),
    period = fifelse(year <= 2059L, "2055-2059", "2095-2099")
  )]
  numerator <- x[, .(attributable_deaths = sum(an)),
                 by = .(branch, period, age_band, range)]
  denominator <- unique(x[, .(branch, year, age, age_band, period, pop)])[
    , .(person_years = sum(pop)), by = .(branch, period, age_band)]
  ans <- merge(numerator, denominator,
               by = c("branch", "period", "age_band"), all.x = TRUE)
  ans[, rate := 1e5 * attributable_deaths / person_years]
  ans[, `:=`(
    rate_signed = fifelse(range %chin% c("ExtrCold", "ModCold"), -rate, rate),
    ssp = ssp,
    geography = label
  )]
  ans
}

plot_data <- rbindlist(lapply(ssp_set, read_ssp))
if (any(!is.finite(plot_data$rate)) || any(plot_data$rate < 0)) {
  stop("Figure 1 rates are missing, non-finite or negative.", call. = FALSE)
}
plot_data[, `:=`(
  age_band = factor(age_band, levels = age_levels),
  period = factor(period, levels = period_levels),
  ssp_label = factor(ssp_labels[as.character(ssp)], levels = ssp_labels[c("1", "2", "3")]),
  range = factor(range, levels = range_levels),
  y = as.numeric(age_band) + fifelse(branch == "with_cc", 0.16, -0.16)
)]

label_data <- plot_data[, .(
  x = max(c(0, rate_signed)) * 0.04,
  y = 1 + fifelse(branch[1] == "with_cc", 0.16, -0.16),
  label = fifelse(branch[1] == "with_cc", "with climate change", "without climate change")
), by = .(ssp_label, period, branch)]

geo_label <- unique(plot_data$geography)
p <- ggplot(plot_data, aes(x = rate_signed, y = y, fill = range, alpha = branch)) +
  geom_col(width = 0.27, colour = NA, orientation = "y") +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#333333") +
  geom_text(
    data = label_data,
    aes(x = x, y = y, label = label), inherit.aes = FALSE,
    size = 2.35, hjust = 0, colour = "#333333"
  ) +
  facet_grid(ssp_label ~ period, scales = "free_x") +
  scale_y_continuous(breaks = seq_along(age_levels), labels = age_levels,
                     expand = expansion(mult = c(0.04, 0.12))) +
  scale_x_continuous(labels = function(z) format(abs(z), big.mark = ",", trim = TRUE)) +
  scale_fill_manual(values = range_colours, labels = range_labels, drop = FALSE) +
  scale_alpha_manual(values = c(with_cc = 1, without_cc = 0.32), guide = "none") +
  labs(
    title = paste0(geo_label, ": absolute risk of temperature-attributable death by age and temperature range"),
    subtitle = "Upper bar with climate change; lower pale bar without climate change; 5-year-period means",
    x = "Attributable deaths per 100,000 person-years (cold to the left; heat to the right)",
    y = "Age group", fill = NULL,
    caption = paste(
      "Static thresholds: ERA5-Land 2000-2014. The moving P97.5 is the observed baseline plus the",
      "change in the calibrated GCM P97.5 over the 15 preceding years; the current year is excluded.",
      "Without climate change the normal-extreme band is empty by construction.",
      "For SSP1, WCDE survival ratios rounded to exactly one were replaced with 0.9995.",
      "Central estimate, no uncertainty intervals."
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "#555555"),
    plot.caption = element_text(size = 8, hjust = 0),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

stem <- file.path(out_dir, "figure1_absolute_risk")
ggsave(paste0(stem, ".png"), p, width = 14, height = 14, units = "in", dpi = 300)
ggsave(paste0(stem, ".pdf"), p, width = 14, height = 14, units = "in", device = cairo_pdf)
fwrite(plot_data, paste0(stem, "_data.csv"))
message("Saved Figure 1 to ", out_dir)
