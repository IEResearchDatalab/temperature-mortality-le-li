################################################################################
#
# Validation: Step 03 — Mortality backbone
#   Diagnose the projected country-level mortality schedules.
#
# Plots:
#   03a_mx_projection.png  — Historical and projected mx for selected ages
#   03b_improvement.png    — Annual improvement rates by age
#   03c_heatmap.png        — Full mx matrix heatmap (age x year)
#
################################################################################

library(ggplot2); library(scales); library(data.table)

validation_dir <- "../results/validation_plots"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Validation Step 03: Mortality backbone\n")
cat("========================================\n\n")

backbone_file <- file.path(path_out, "mortality_backbone.rds")
if (!file.exists(backbone_file)) {
  cat("  WARNING: mortality_backbone.rds not found. Run 03_MortalityBackbone.R first.\n")
  quit(save = "no")
}

backbone <- readRDS(backbone_file)
lt_hist <- readRDS(file.path(path_out, "lt_eu_raw_filtered.rds"))

cat(sprintf("  Backbone rows: %s\n", format(nrow(backbone), big.mark = ",")))
cat(sprintf("  Countries: %s\n", paste(sort(unique(backbone$country)), collapse = ", ")))
cat(sprintf("  Ages: %d to %d\n", min(backbone$age), max(backbone$age)))
cat(sprintf("  Years: %d to %d\n", min(backbone$year), max(backbone$year)))
cat(sprintf("  Baseline year range: %d-%d\n",
    min(backbone$baseline_year, na.rm = TRUE),
    max(backbone$baseline_year, na.rm = TRUE)))

#---------------------------
# 03a: Historical + projected mx for Spain
#---------------------------

spain_hist <- lt_hist[country == "ES" & age %in% c(45, 65, 75, 85)]
spain_proj <- backbone[country == "ES" & age %in% c(45, 65, 75, 85)]

p <- ggplot() +
  geom_line(data = spain_hist, aes(x = year, y = mx, color = factor(age)),
    linewidth = 1, linetype = "solid") +
  geom_line(data = spain_proj, aes(x = year, y = mx_proj, color = factor(age)),
    linewidth = 1, linetype = "dashed") +
  scale_y_log10(labels = comma) +
  labs(title = "Spain: Historical and projected mortality rates",
       subtitle = "Solid = historical (Eurostat). Dashed = projected (log-linear trend).",
       caption = "Mortality backbone method: log-linear trend extrapolation for each age.
Baseline year = most recent historical year with complete data.",
       x = "Year", y = "Mortality rate (mx)", color = "Age") +
  geom_vline(xintercept = 2019, linetype = "dotted", alpha = 0.5) +
  annotate("text", x = 2019, y = Inf, label = "end of historical", vjust = 1.5,
    hjust = 1.1, size = 3, alpha = 0.7) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "03a_mx_projection.png"), p, width = 10, height = 6)
cat(sprintf("  03a: mx projection saved\n"))

#---------------------------
# 03b: Annual improvement rates by age (Spain)
#---------------------------

spain_rates <- backbone[country == "ES", .(b = unique(b)), by = age]
# Convert log-linear slope to annual % change: (exp(b) - 1) * 100
spain_rates[, pct_change := (exp(b) - 1) * 100]

p <- ggplot(spain_rates[age >= 40 & age <= 95], aes(x = age, y = pct_change)) +
  geom_line(linewidth = 1, color = "#2980b9") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Spain: Annual mortality improvement rate by age",
       subtitle = "Negative = mortality declining. Derived from log-linear trend 1990-2019.",
       caption = "Improvement rate = exp(b) - 1 where b is the slope of ln(mx) ~ year.
Typical range: -3% to +1% per year depending on age.",
       x = "Age", y = "Annual change in mx (%)") +
  theme_minimal()

ggsave(file.path(validation_dir, "03b_improvement.png"), p, width = 9, height = 5)
cat(sprintf("  03b: improvement rates saved\n"))

#---------------------------
# 03c: mx heatmap (Spain, age x year)
#---------------------------

spain_heat <- backbone[country == "ES" & year >= 2020 & year <= 2100 &
  age >= 65 & age <= 100]

p <- ggplot(spain_heat, aes(x = year, y = age, fill = mx_proj)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", labels = comma) +
  labs(title = "Spain: Projected mortality rates (age x year heatmap)",
       subtitle = "Log scale. Age 65+ only.",
       caption = "Each tile = mx for one age-year combination.
Dark = low mortality, light = high mortality.",
       x = "Year", y = "Age", fill = "mx") +
  theme_minimal()

ggsave(file.path(validation_dir, "03c_heatmap.png"), p, width = 10, height = 7)
cat(sprintf("  03c: heatmap saved\n"))

cat("\n--- Step 03 validation complete ---\n")