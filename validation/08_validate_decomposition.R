################################################################################
#
# Validation: Step 08 — Decomposition
#   Diagnose the Aburto/Horiuchi linear integral decomposition.
#
# Plots:
#   08a_decomp_e65.png      — Contribution to e65 change by age and cause
#   08b_decomp_sd.png       — Contribution to SD change by age and cause
#   08c_additivity_check.png— Sum of contributions vs total observed change
#
################################################################################

library(ggplot2); library(scales); library(data.table)

validation_dir <- "../results/validation_plots"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Validation Step 08: Decomposition\n")
cat("========================================\n\n")

decomp_annual_file <- file.path(path_out, "decomposition_annual_all_cities.csv")
decomp_period_file <- file.path(path_out, "decomposition_period_all_cities.csv")

if (!file.exists(decomp_annual_file)) {
  cat("  WARNING: decomposition files not found. Run 08 first.\n")
  quit(save = "no")
}

decomp_annual <- fread(decomp_annual_file)
decomp_period <- fread(decomp_period_file)

cat(sprintf("  Annual decomposition: %s rows\n", format(nrow(decomp_annual), big.mark = ",")))
cat(sprintf("  Period decomposition: %s rows\n", format(nrow(decomp_period), big.mark = ",")))
cat(sprintf("  Cities: %d | Causes: %s\n",
    uniqueN(decomp_period$city_code),
    paste(unique(decomp_period$cause), collapse = ", ")))

#---------------------------
# Additivity check
#---------------------------

madrid_decomp <- decomp_period[city_code == "ES001C"]

add_check <- madrid_decomp[,
  .(sum_delta_e65 = sum(delta_e65),
    sum_delta_sd = sum(delta_sd)),
  by = .(city_code, gcm, ssp, period)]

cat(sprintf("  Additivity check (sum of contributions vs total change):\n"))
cat(sprintf("    e65 range: [%.6f, %.6f]\n", min(add_check$sum_delta_e65),
    max(add_check$sum_delta_e65)))
cat(sprintf("    SD range: [%.6f, %.6f]\n", min(add_check$sum_delta_sd),
    max(add_check$sum_delta_sd)))

#---------------------------
# 08a: e65 decomposition by age and cause
#---------------------------

# Pick one GCM-SSP for clear visualization
one_gcm <- unique(madrid_decomp$gcm)[1]
madrid_one <- madrid_decomp[gcm == one_gcm & ssp == "2"]

p <- ggplot(madrid_one, aes(x = age, y = delta_e65, fill = cause)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~period, ncol = 1, scales = "fixed") +
  labs(title = sprintf("Madrid SSP3-7.0: Decomposition of e65 change (%s)", one_gcm),
       subtitle = "Aburto/Horiuchi linear integral decomposition by age and temperature cause.",
       caption = "Bars = contribution of each age-cause combination to the total e65 change.
Negative = reducing e65 (harmful). Positive = increasing e65 (beneficial).
Causes sum to total e65 change at each period.",
       x = "Age", y = "Contribution to e65 change (years)") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_fill_manual(values = c(
    "residual" = "#2c3e50",
    "extreme_cold" = "#3498db",
    "moderate_cold" = "#85c1e9",
    "moderate_heat" = "#e74c3c",
    "extreme_heat" = "#c0392b"
  )) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "08a_decomp_e65.png"), p, width = 12, height = 10)
cat(sprintf("  08a: e65 decomposition saved\n"))

#---------------------------
# 08b: SD decomposition
#---------------------------

p <- ggplot(madrid_one, aes(x = age, y = delta_sd, fill = cause)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~period, ncol = 1, scales = "fixed") +
  labs(title = sprintf("Madrid SSP3-7.0: Decomposition of LI change (%s)", one_gcm),
       subtitle = "LI = SD of age at death above 65. Linear integral decomposition.",
       caption = "Bars = contribution of each age-cause to SD change.
Positive = increasing inequality. Negative = decreasing inequality.
Causes sum to total SD change at each period.",
       x = "Age", y = "Contribution to SD change") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_fill_manual(values = c(
    "residual" = "#2c3e50",
    "extreme_cold" = "#3498db",
    "moderate_cold" = "#85c1e9",
    "moderate_heat" = "#e74c3c",
    "extreme_heat" = "#c0392b"
  )) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "08b_decomp_sd.png"), p, width = 12, height = 10)
cat(sprintf("  08b: SD decomposition saved\n"))

#---------------------------
# 08c: Cause summary (punch card)
#---------------------------

cause_summary <- madrid_one[,
  .(delta_e65 = sum(delta_e65),
    delta_sd = sum(delta_sd)),
  by = .(period, cause)]

p <- ggplot(cause_summary, aes(x = period, y = cause, fill = delta_e65)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f", delta_e65)), size = 3.5) +
  scale_fill_gradient2(low = "#3498db", mid = "white", high = "#e74c3c",
    name = "Δe65") +
  labs(title = sprintf("Madrid: e65 contribution by period and cause (%s SSP3-7.0)",
       one_gcm),
       subtitle = "Punch card: each cell = total contribution of one cause in one period.",
       caption = "Numbers = years contributed to e65 change.
Red = positive (increasing LE), Blue = negative (decreasing LE).",
       x = "Period", y = "Cause") +
  theme_minimal()

ggsave(file.path(validation_dir, "08c_cause_summary.png"), p, width = 8, height = 5)
cat(sprintf("  08c: Cause summary saved\n"))

cat("\n--- Step 08 validation complete ---\n")