################################################################################
# 
# Validation Script: Compare 4-Range Disaggregation Against Masselot Baseline
#
# Goal: Verify that ExtrCold + ModCold = excess_cold and
#                     ExtrHeat + ModHeat = excess_heat
#       from Masselot et al. 2025 baseline period (2000-2014)
#
################################################################################

library(data.table)
library(dplyr)
library(ggplot2)

message("\n=== MASSELOT REPRODUCIBILITY VALIDATION ===\n")

#----- Load Masselot baseline results (2000-2014 historical period)

masselot_baseline <- fread("data/city_results.csv")
masselot_ref <- masselot_baseline[, .(
  URAU_CODE, 
  agegroup, 
  masselot_cold = excess_cold_est,
  masselot_heat = excess_heat_est,
  masselot_total = excess_total_est
)]

message("Loaded Masselot baseline for ", length(unique(masselot_ref$URAU_CODE)), 
        " cities and ", length(unique(masselot_ref$agegroup)), " age groups")

#----- Load our temp_results and aggregate for 2000-2014

# Identify all city RDS files
rds_files <- list.files("temp_results", pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " city result files")

# Process each city
our_results_list <- list()

pb <- txtProgressBar(max = length(rds_files), style = 3)

for (i in seq_along(rds_files)) {
  city_file <- rds_files[i]
  city_id <- gsub(".rds", "", basename(city_file))
  
  d <- try(readRDS(city_file), silent = TRUE)
  if (inherits(d, "try-error")) next
  
  setDT(d)
  
  # Filter for 2000-2014 period (matching Masselot's histrange)
  d_hist <- d[year >= 2000 & year <= 2014]
  
  if (nrow(d_hist) == 0) next
  
  # Aggregate by agegroup, range, ssp, gcm across years and simulations
  # Take mean across simulations (sim column)
  d_agg <- d_hist[, .(an_mean = mean(an)), by = .(agegroup, range, ssp, gcm, year)]
  
  # Sum across all years in the period to get total AN for 2000-2014
  d_total <- d_agg[, .(an_total = sum(an_mean)), by = .(agegroup, range, ssp, gcm)]
  
  # Average across GCMs and SSPs to get ensemble mean (matching Masselot's approach)
  d_ensemble <- d_total[, .(an_ensemble = mean(an_total)), by = .(agegroup, range)]
  
  # Add city ID
  d_ensemble[, URAU_CODE := city_id]
  
  our_results_list[[length(our_results_list) + 1]] <- d_ensemble
  
  setTxtProgressBar(pb, i)
}

close(pb)

our_results <- rbindlist(our_results_list)

message("\nProcessed results for ", length(unique(our_results$URAU_CODE)), " cities")

#----- Reshape to compare Cold and Heat totals

# Sum the 4 ranges into Cold (ExtrCold + ModCold) and Heat (ModHeat + ExtrHeat)
our_totals <- our_results[, .(
  our_cold = sum(an_ensemble[range %in% c("ExtrCold", "ModCold")]),
  our_heat = sum(an_ensemble[range %in% c("ModHeat", "ExtrHeat")]),
  our_extrcold = sum(an_ensemble[range == "ExtrCold"]),
  our_modcold = sum(an_ensemble[range == "ModCold"]),
  our_modheat = sum(an_ensemble[range == "ModHeat"]),
  our_extrheat = sum(an_ensemble[range == "ExtrHeat"])
), by = .(URAU_CODE, agegroup)]

#----- Merge with Masselot baseline

comparison <- merge(our_totals, masselot_ref, by = c("URAU_CODE", "agegroup"), all = TRUE)

# Remove rows with missing data
comparison <- comparison[complete.cases(comparison)]

message("Merged data for ", nrow(comparison), " city-age combinations")

#----- Calculate differences and relative errors

comparison[, `:=`(
  diff_cold = our_cold - masselot_cold,
  diff_heat = our_heat - masselot_heat,
  rel_error_cold = (our_cold - masselot_cold) / abs(masselot_cold) * 100,
  rel_error_heat = (our_heat - masselot_heat) / abs(masselot_heat) * 100
)]

#----- Summary statistics

message("\n--- COLD-RELATED MORTALITY ---")
message(sprintf("Mean absolute difference: %.2f deaths", mean(abs(comparison$diff_cold), na.rm=TRUE)))
message(sprintf("Median absolute difference: %.2f deaths", median(abs(comparison$diff_cold), na.rm=TRUE)))
message(sprintf("Mean relative error: %.2f%%", mean(abs(comparison$rel_error_cold), na.rm=TRUE)))
message(sprintf("Median relative error: %.2f%%", median(abs(comparison$rel_error_cold), na.rm=TRUE)))

cold_concordance <- cor(comparison$our_cold, comparison$masselot_cold, use = "complete.obs")
message(sprintf("Pearson correlation: %.4f", cold_concordance))

message("\n--- HEAT-RELATED MORTALITY ---")
message(sprintf("Mean absolute difference: %.2f deaths", mean(abs(comparison$diff_heat), na.rm=TRUE)))
message(sprintf("Median absolute difference: %.2f deaths", median(abs(comparison$diff_heat), na.rm=TRUE)))
message(sprintf("Mean relative error: %.2f%%", mean(abs(comparison$rel_error_heat), na.rm=TRUE)))
message(sprintf("Median relative error: %.2f%%", median(abs(comparison$rel_error_heat), na.rm=TRUE)))

heat_concordance <- cor(comparison$our_heat, comparison$masselot_heat, use = "complete.obs")
message(sprintf("Pearson correlation: %.4f", heat_concordance))

#----- Identify largest discrepancies

message("\n--- TOP 10 DISCREPANCIES (Cold) ---")
top_cold <- comparison[order(-abs(diff_cold))][1:10]
print(top_cold[, .(URAU_CODE, agegroup, our_cold, masselot_cold, diff_cold, rel_error_cold)])

message("\n--- TOP 10 DISCREPANCIES (Heat) ---")
top_heat <- comparison[order(-abs(diff_heat))][1:10]
print(top_heat[, .(URAU_CODE, agegroup, our_heat, masselot_heat, diff_heat, rel_error_heat)])

#----- Save validation results

fwrite(comparison, "data/validation_masselot_comparison.csv")
message("\nFull comparison saved to: data/validation_masselot_comparison.csv")

#----- Create validation plots

dir.create("figures", showWarnings = FALSE)

# Scatter plot: Cold
p1 <- ggplot(comparison, aes(x = masselot_cold, y = our_cold)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_point(alpha = 0.3, size = 1) +
  labs(
    title = "Cold-Related Mortality: 4-Range Disaggregation vs Masselot Baseline",
    subtitle = sprintf("2000-2014 period | n=%d city-age groups | r=%.3f", 
                      nrow(comparison), cold_concordance),
    x = "Masselot Baseline (deaths)",
    y = "Our 4-Range Sum (ExtrCold + ModCold)"
  ) +
  theme_minimal() +
  coord_fixed()

# Scatter plot: Heat
p2 <- ggplot(comparison, aes(x = masselot_heat, y = our_heat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_point(alpha = 0.3, size = 1) +
  labs(
    title = "Heat-Related Mortality: 4-Range Disaggregation vs Masselot Baseline",
    subtitle = sprintf("2000-2014 period | n=%d city-age groups | r=%.3f", 
                      nrow(comparison), heat_concordance),
    x = "Masselot Baseline (deaths)",
    y = "Our 4-Range Sum (ModHeat + ExtrHeat)"
  ) +
  theme_minimal() +
  coord_fixed()

ggsave("figures/validation_cold_scatter.png", p1, width = 8, height = 7, dpi = 300)
ggsave("figures/validation_heat_scatter.png", p2, width = 8, height = 7, dpi = 300)

message("\nValidation plots saved to figures/")

#----- Breakdown by range component

message("\n--- COLD DISAGGREGATION (% of total) ---")
comparison[, cold_total_ours := our_extrcold + our_modcold]
comparison[, `:=`(
  pct_extrcold = our_extrcold / cold_total_ours * 100,
  pct_modcold = our_modcold / cold_total_ours * 100
)]

message(sprintf("Extreme Cold: %.1f%% (mean)", mean(comparison$pct_extrcold, na.rm=TRUE)))
message(sprintf("Moderate Cold: %.1f%% (mean)", mean(comparison$pct_modcold, na.rm=TRUE)))

message("\n--- HEAT DISAGGREGATION (% of total) ---")
comparison[, heat_total_ours := our_modheat + our_extrheat]
comparison[, `:=`(
  pct_modheat = our_modheat / heat_total_ours * 100,
  pct_extrheat = our_extrheat / heat_total_ours * 100
)]

message(sprintf("Moderate Heat: %.1f%% (mean)", mean(comparison$pct_modheat, na.rm=TRUE)))
message(sprintf("Extreme Heat: %.1f%% (mean)", mean(comparison$pct_extrheat, na.rm=TRUE)))

#----- Final verdict

message("\n=== VALIDATION SUMMARY ===")
if (cold_concordance > 0.95 && heat_concordance > 0.95) {
  message("✓ EXCELLENT: Strong agreement with Masselot baseline (r > 0.95)")
} else if (cold_concordance > 0.90 && heat_concordance > 0.90) {
  message("✓ GOOD: Good agreement with Masselot baseline (r > 0.90)")
} else if (cold_concordance > 0.80 && heat_concordance > 0.80) {
  message("⚠ MODERATE: Moderate agreement with Masselot baseline (r > 0.80)")
  message("  Consider investigating systematic differences")
} else {
  message("✗ POOR: Low agreement with Masselot baseline (r < 0.80)")
  message("  Major differences detected - review methodology")
}

message("\nValidation complete!")
