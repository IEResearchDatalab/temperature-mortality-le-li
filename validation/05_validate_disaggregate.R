################################################################################
#
# Validation: Step 05 — AN disaggregation (PCLM)
#   Diagnose the PCLM-based disaggregation from grouped to single ages.
#
# Plots:
#   05a_disagg_profile.png   — Single-age AN profile for selected year/SSP
#   05b_grouped_vs_single.png— Comparison before/after disaggregation
#   05c_disagg_by_range.png  — Single-age ANs split by temperature range
#
################################################################################

library(ggplot2); library(scales); library(data.table)

validation_dir <- "../results/validation_plots"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Validation Step 05: AN disaggregation\n")
cat("========================================\n\n")

single_file <- file.path(path_out, "ans_single_age_all_cities.csv")
if (!file.exists(single_file)) {
  cat("  WARNING: ans_single_age_all_cities.csv not found. Run 05 first.\n")
  quit(save = "no")
}

single_an <- fread(single_file)
cat(sprintf("  Single-age AN rows: %s\n", format(nrow(single_an), big.mark = ",")))
cat(sprintf("  Cities: %d | Ages: %d-%d | Temp ranges: %s\n",
    uniqueN(single_an$city_code), min(single_an$age), max(single_an$age),
    paste(unique(single_an$temp_range), collapse = ", ")))

#---------------------------
# 05a: Single-age AN profile (Madrid 2050, SSP 2, moderate heat)
#---------------------------

madrid_2050 <- single_an[city_code == "ES001C" & year == 2050 & ssp == "2"]
madrid_agg <- madrid_2050[, .(AN = mean(AN)), by = .(age, temp_range)]

p <- ggplot(madrid_agg, aes(x = age, y = AN, fill = temp_range)) +
  geom_bar(stat = "identity", width = 1) +
  labs(title = "Madrid 2050 SSP3-7.0: Single-age ANs by temperature range",
       subtitle = "PCLM disaggregation. Averages across GCMs.",
       caption = "Each bar = attributable deaths at one single age from one temperature range.
PCLM preserves grouped totals and keeps age distribution smooth.",
       x = "Age", y = "Attributable deaths", fill = "Temp range") +
  scale_y_continuous(labels = comma) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "05a_disagg_profile.png"), p, width = 10, height = 6)
cat(sprintf("  05a: Disaggregation profile saved\n"))

#---------------------------
# 05b: Grouped total preservation check
#---------------------------

# Re-load grouped ANs for comparison
grouped_an <- fread(file.path(path_out, "ans_annual_all_cities.csv"))
grouped_total <- grouped_an[, .(total_grouped = sum(an_est)),
  by = .(city_code, gcm, ssp, year, temp_range)]
single_total <- single_an[, .(total_single = sum(AN)),
  by = .(city_code, gcm, ssp, year, temp_range)]

comp <- merge(grouped_total, single_total,
  by = c("city_code", "gcm", "ssp", "year", "temp_range"))
comp[, diff := total_grouped - total_single]
comp[, pct_diff := 100 * abs(diff) / pmax(total_grouped, 1)]

cat(sprintf("  Total preservation check:\n"))
cat(sprintf("    Max abs diff: %.6f\n", max(abs(comp$diff))))
cat(sprintf("    Median pct diff: %.4f%%\n", median(comp$pct_diff, na.rm = TRUE)))
cat(sprintf("    Rows with >1%% diff: %d / %d\n",
    nrow(comp[comp$pct_diff > 1]), nrow(comp)))

p <- ggplot(comp, aes(x = pct_diff)) +
  geom_histogram(bins = 50, fill = "#3498db", alpha = 0.7) +
  labs(title = "PCLM total preservation: grouped vs single-age totals",
       subtitle = "Difference should be near zero if PCLM preserves group totals.",
       caption = "Difference = |grouped_total - single_total| / grouped_total × 100%.
Values >1% indicate PCLM convergence issues.",
       x = "Percent difference", y = "Number of combinations") +
  scale_x_log10(labels = comma) +
  theme_minimal()

ggsave(file.path(validation_dir, "05b_grouped_vs_single.png"), p, width = 9, height = 5)
cat(sprintf("  05b: Total preservation saved\n"))

cat("\n--- Step 05 validation complete ---\n")