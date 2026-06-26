################################################################################
#
# Validation: Step 06 — Cause-specific bookkeeping
#   Diagnose the identity: all-cause = residual + temp components.
#
# Plots:
#   06a_death_decomp.png    — All-cause deaths decomposed by component over time
#   06b_residual_check.png  — Residual age profile (should be smooth, no negatives)
#   06c_identity_check.png  — Verification that identity holds
#
################################################################################

library(ggplot2); library(scales); library(data.table)

validation_dir <- "../results/validation_plots"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Validation Step 06: Cause-specific bookkeeping\n")
cat("========================================\n\n")

analysis_file <- file.path(path_out, "analysis_dataset_all_cities.csv")
if (!file.exists(analysis_file)) {
  cat("  WARNING: analysis_dataset_all_cities.csv not found. Run 06 first.\n")
  quit(save = "no")
}

analysis <- fread(analysis_file)
cat(sprintf("  Analysis dataset rows: %s\n", format(nrow(analysis), big.mark = ",")))
cat(sprintf("  Cities: %d | GCMs: %d | SSPs: %s\n",
    uniqueN(analysis$city_code), uniqueN(analysis$gcm),
    paste(unique(analysis$ssp), collapse = ", ")))

#---------------------------
# Identity check
#---------------------------

analysis[, identity_check := deaths - (residual_deaths + extreme_cold +
    moderate_cold + moderate_heat + extreme_heat)]
max_err <- max(abs(analysis$identity_check))
n_fail <- sum(abs(analysis$identity_check) > 0.01)
n_neg <- sum(analysis$residual_deaths < -0.001)

cat(sprintf("  Identity check:\n"))
cat(sprintf("    Max error: %.10f\n", max_err))
cat(sprintf("    Rows with error > 0.01: %d / %s\n", n_fail,
    format(nrow(analysis), big.mark = ",")))
cat(sprintf("    Negative residuals: %d\n", n_neg))

if (n_fail > 0) {
  cat("  WARNING: Identity failures detected. Investigate!\n")
}
if (n_neg > 0) {
  cat("  WARNING: Negative residuals found. Bookkeeping issue!\n")
}

#---------------------------
# 06a: Death decomposition (Madrid, one GCM-SSP)
#---------------------------

madrid <- analysis[city_code == "ES001C" & gcm == unique(analysis$gcm)[1] & ssp == "2"]
madrid_ts <- madrid[, .(
  deaths = sum(deaths),
  residual = sum(residual_deaths),
  extreme_cold = sum(extreme_cold),
  moderate_cold = sum(moderate_cold),
  moderate_heat = sum(moderate_heat),
  extreme_heat = sum(extreme_heat)
), by = year]

# Check identity holds at aggregate level
madrid_ts[, check := deaths - (residual + extreme_cold + moderate_cold +
    moderate_heat + extreme_heat)]
cat(sprintf("  Madrid aggregate identity max error: %.10f\n", max(abs(madrid_ts$check))))

p <- ggplot(madrid_ts, aes(x = year)) +
  geom_area(aes(y = residual, fill = "Non-temperature"), alpha = 0.8) +
  geom_area(aes(y = residual + extreme_cold, fill = "Extreme cold"), alpha = 0.8) +
  geom_area(aes(y = residual + extreme_cold + moderate_cold,
    fill = "Moderate cold"), alpha = 0.8) +
  geom_area(aes(y = residual + extreme_cold + moderate_cold + moderate_heat,
    fill = "Moderate heat"), alpha = 0.8) +
  geom_area(aes(y = deaths, fill = "Extreme heat"), alpha = 0.8) +
  labs(title = "Madrid SSP3-7.0: Death decomposition by component",
       subtitle = sprintf("All-cause = non-temperature + 4 temperature ranges (%s)",
         unique(analysis$gcm)[1]),
       caption = "Stacked areas show cumulative contribution of each component.
The top boundary of each area = all-cause deaths up to that component.
Non-temperature residual = all-cause minus all temperature-attributable deaths.",
       x = "Year", y = "Annual deaths") +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c(
    "Non-temperature" = "#2c3e50",
    "Extreme cold" = "#3498db",
    "Moderate cold" = "#85c1e9",
    "Moderate heat" = "#e74c3c",
    "Extreme heat" = "#c0392b"
  )) +
  guides(fill = guide_legend(title = "Component")) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "06a_death_decomp.png"), p, width = 12, height = 7)
cat(sprintf("  06a: Death decomposition saved\n"))

#---------------------------
# 06b: Residual age profile
#---------------------------

madrid_2050 <- madrid[year == 2050]
p <- ggplot(madrid_2050, aes(x = age)) +
  geom_line(aes(y = residual_deaths, color = "Non-temperature"), linewidth = 1) +
  geom_line(aes(y = extreme_cold, color = "Extreme cold"), linewidth = 1) +
  geom_line(aes(y = moderate_cold, color = "Moderate cold"), linewidth = 1) +
  geom_line(aes(y = moderate_heat, color = "Moderate heat"), linewidth = 1) +
  geom_line(aes(y = extreme_heat, color = "Extreme heat"), linewidth = 1) +
  labs(title = "Madrid 2050 SSP3-7.0: Deaths by age and component",
       subtitle = "Single-age resolution showing where temperature effects concentrate.",
       caption = "Non-temperature residual should be positive and smooth.
Temperature effects peak at older ages where baseline mortality is higher.",
       x = "Age", y = "Deaths") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c(
    "Non-temperature" = "#2c3e50",
    "Extreme cold" = "#3498db",
    "Moderate cold" = "#85c1e9",
    "Moderate heat" = "#e74c3c",
    "Extreme heat" = "#c0392b"
  )) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "06b_residual_check.png"), p, width = 10, height = 6)
cat(sprintf("  06b: Residual age profile saved\n"))

#---------------------------
# 06c: Identity check (age profile)
#---------------------------

madrid_2050[, identity_val := deaths - (residual_deaths + extreme_cold +
    moderate_cold + moderate_heat + extreme_heat)]
p <- ggplot(madrid_2050, aes(x = age, y = identity_val)) +
  geom_line(color = "#27ae60", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Madrid 2050: Bookkeeping identity check",
       subtitle = "deaths - (residual + extreme_cold + moderate_cold + moderate_heat + extreme_heat)",
       caption = "Should be zero at all ages. Non-zero values indicate bookkeeping failure.
Values within ±0.001 are floating-point acceptable.",
       x = "Age", y = "Identity error") +
  theme_minimal()

ggsave(file.path(validation_dir, "06c_identity_check.png"), p, width = 9, height = 5)
cat(sprintf("  06c: Identity check saved\n"))

cat("\n--- Step 06 validation complete ---\n")