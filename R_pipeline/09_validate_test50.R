################################################################################
# 
# TEST VALIDATION: Compare 50 test cities with Masselot baseline
#
# This script validates whether changing hist_years to 2000:2014 
# improves agreement with Masselot et al. (2025)
#
################################################################################

library(data.table)
library(ggplot2)

message("\n[TEST] Starting city-level validation for 50 test cities...")
message("[TEST] Comparing temp_results_test50/ against Masselot baseline")

#----- Load Masselot baseline data

masselot_baseline <- fread("references/2025-masselot-zenodo/results/cityage.csv")
setnames(masselot_baseline, 
         c("excess_cold_est", "excess_heat_est"), 
         c("masselot_cold", "masselot_heat"))

message("[TEST] Loaded Masselot baseline: ", nrow(masselot_baseline), " city-age combinations")

#----- Process test results

test_files <- list.files("temp_results_test50", pattern = "\\.rds$", full.names = TRUE)
message("[TEST] Found ", length(test_files), " test result files")

if(length(test_files) == 0) {
  stop("[TEST] ERROR: No test results found. Run 03_attribution_test50.R first.")
}

comparison_list <- list()

for(rds_file in test_files) {
  city_id <- gsub(".rds$", "", basename(rds_file))
  d <- readRDS(rds_file)
  
  # Filter to historical period 2000-2014
  d_hist <- d[year >= 2000 & year <= 2014]
  
  if(nrow(d_hist) == 0) next
  
  # Average across simulations first
  d_sim <- d_hist[, .(an_mean = mean(an)), 
                  by = .(agegroup, range, ssp, gcm, year)]
  
  # Average across years (annual average, not cumulative)
  d_year <- d_sim[, .(an_annual = mean(an_mean)), 
                  by = .(agegroup, range, ssp, gcm)]
  
  # Average across SSPs and GCMs (ensemble mean for historical period)
  d_ens <- d_year[, .(an_ens = mean(an_annual)), 
                  by = .(agegroup, range)]
  
  # Sum to cold/heat
  d_wide <- dcast(d_ens, agegroup ~ range, value.var = "an_ens", fill = 0)
  d_wide[, `:=`(
    our_cold = ExtrCold + ModCold,
    our_heat = ModHeat + ExtrHeat,
    URAU_CODE = city_id
  )]
  
  comparison_list[[city_id]] <- d_wide[, .(URAU_CODE, agegroup, our_cold, our_heat)]
}

our_results <- rbindlist(comparison_list)
message("[TEST] Processed ", uniqueN(our_results$URAU_CODE), " cities")

#----- Merge with Masselot baseline

comparison <- merge(our_results, 
                    masselot_baseline[, .(URAU_CODE, agegroup, masselot_cold, masselot_heat)],
                    by = c("URAU_CODE", "agegroup"),
                    all = FALSE)

message("[TEST] Matched ", nrow(comparison), " city-age combinations")

#----- Calculate validation metrics

# Cold mortality
cold_cor <- cor(comparison$our_cold, comparison$masselot_cold, use = "complete.obs")
cold_rel_error <- abs(comparison$our_cold - comparison$masselot_cold) / abs(comparison$masselot_cold)
cold_median_error <- median(cold_rel_error[is.finite(cold_rel_error)], na.rm = TRUE) * 100

# Heat mortality
heat_cor <- cor(comparison$our_heat, comparison$masselot_heat, use = "complete.obs")
heat_rel_error <- abs(comparison$our_heat - comparison$masselot_heat) / abs(comparison$masselot_heat)
heat_median_error <- median(heat_rel_error[is.finite(heat_rel_error)], na.rm = TRUE) * 100

#----- Print results

cat("\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("TEST VALIDATION RESULTS (50 cities)\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("\n")
cat("COLD MORTALITY:\n")
cat("  Correlation (r):        ", sprintf("%.4f", cold_cor), "\n")
cat("  Median Relative Error:  ", sprintf("%.2f%%", cold_median_error), "\n")
cat("\n")
cat("HEAT MORTALITY:\n")
cat("  Correlation (r):        ", sprintf("%.4f", heat_cor), "\n")
cat("  Median Relative Error:  ", sprintf("%.2f%%", heat_median_error), "\n")
cat("\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("\n")

#----- Create output directory

dir.create("results/test50_validation", showWarnings = FALSE, recursive = TRUE)

#----- Save comparison data

fwrite(comparison, "results/test50_validation/test_comparison.csv")
message("[TEST] Saved comparison data to results/test50_validation/test_comparison.csv")

#----- Create scatter plots

# Cold mortality
p_cold <- ggplot(comparison, aes(x = masselot_cold, y = our_cold)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = sprintf("TEST: Cold Mortality (50 cities) - r = %.4f", cold_cor),
    subtitle = sprintf("Median relative error: %.2f%%", cold_median_error),
    x = "Masselot et al. (2025) - Cold AN",
    y = "Our Results - Cold AN"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("results/test50_validation/test_scatter_cold.png", p_cold, width = 8, height = 6, dpi = 300)

# Heat mortality
p_heat <- ggplot(comparison, aes(x = masselot_heat, y = our_heat)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = sprintf("TEST: Heat Mortality (50 cities) - r = %.4f", heat_cor),
    subtitle = sprintf("Median relative error: %.2f%%", heat_median_error),
    x = "Masselot et al. (2025) - Heat AN",
    y = "Our Results - Heat AN"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("results/test50_validation/test_scatter_heat.png", p_heat, width = 8, height = 6, dpi = 300)

message("[TEST] Saved scatter plots to results/test50_validation/")

#----- Write summary

summary_text <- sprintf("
TEST VALIDATION SUMMARY (50 cities)
====================================

Historical Period Used: 2000-2014 (matching Masselot)

RESULTS:
--------
Cold Mortality:
  - Correlation: r = %.4f
  - Median Relative Error: %.2f%%

Heat Mortality:
  - Correlation: r = %.4f
  - Median Relative Error: %.2f%%

Cities Tested: %d
City-Age Combinations: %d

COMPARISON TO PREVIOUS RESULTS (854 cities, hist_years=1990:2014):
------------------------------------------------------------------
Previous Cold: r = 0.9746, error = 44.69%%
Previous Heat: r = 0.9274, error = 47.81%%

STATUS: %s

", 
cold_cor, cold_median_error,
heat_cor, heat_median_error,
uniqueN(comparison$URAU_CODE),
nrow(comparison),
ifelse(cold_median_error < 40 & heat_median_error < 40, 
       "IMPROVEMENT DETECTED - Proceed with full rerun",
       "NO SIGNIFICANT IMPROVEMENT - Investigate further"))

cat(summary_text)
writeLines(summary_text, "results/test50_validation/TEST_SUMMARY.txt")
message("[TEST] Saved summary to results/test50_validation/TEST_SUMMARY.txt")

message("\n[TEST] Validation complete!")
