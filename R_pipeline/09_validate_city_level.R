################################################################################
# 
# CITY-LEVEL Validation: Direct City-by-City Comparison with Masselot
#
# Goal: Compare our 4-range disaggregation city-by-city against Masselot's
#       cityage-level baseline (2000-2014) WITHOUT any aggregation
#
# Approach: For each city-age combination:
#   - Our Cold = ExtrCold + ModCold
#   - Our Heat = ModHeat + ExtrHeat
#   - Compare directly to excess_cold_est and excess_heat_est from Masselot
#
################################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(doSNOW)
library(foreach)

# Create output directory
dir.create("results/masselot_validation_city", recursive = TRUE, showWarnings = FALSE)

message("\n=== CITY-LEVEL VALIDATION (NO AGGREGATION) ===\n")

#----- Load Masselot city-age baseline (2000-2014)

masselot_cityage <- fread("references/2025-masselot-zenodo/results/cityage.csv")

masselot_ref <- masselot_cityage[, .(
  URAU_CODE, 
  agegroup,
  CNTR_CODE,
  cntr_name,
  masselot_cold = excess_cold_est,
  masselot_heat = excess_heat_est,
  masselot_total = excess_total_est
)]

message("Loaded Masselot baseline:")
message("  - Cities: ", length(unique(masselot_ref$URAU_CODE)))
message("  - Age groups: ", paste(unique(masselot_ref$agegroup), collapse=", "))
message("  - Total combinations: ", nrow(masselot_ref))

#----- Load our temp_results for 2000-2014 period

rds_files <- list.files("temp_results", pattern = "\\.rds$", full.names = TRUE)
message("\nFound ", length(rds_files), " city result files")

# Function to process each city
process_city_file <- function(city_file) {
  city_id <- gsub(".rds", "", basename(city_file))
  
  d <- try(readRDS(city_file), silent = TRUE)
  if (inherits(d, "try-error")) return(NULL)
  
  setDT(d)
  
  # Filter for 2000-2014 period ONLY
  d_hist <- d[year >= 2000 & year <= 2014]
  
  if (nrow(d_hist) == 0) return(NULL)
  
  # Average across simulations first
  d_sim <- d_hist[, .(an_mean = mean(an)), 
                  by = .(agegroup, range, ssp, gcm, year)]
  
  # Average across years (annual average, not cumulative)
  d_year <- d_sim[, .(an_annual = mean(an_mean)), 
                  by = .(agegroup, range, ssp, gcm)]
  
  # Average across SSPs and GCMs (ensemble mean for historical period)
  # This matches Masselot's approach of providing a single baseline estimate
  d_ens <- d_year[, .(an_ens = mean(an_annual)), 
                  by = .(agegroup, range)]
  
  # Reshape to wide format
  d_wide <- dcast(d_ens, agegroup ~ range, value.var = "an_ens", fill = 0)
  
  # Calculate total cold and heat
  d_wide[, `:=`(
    our_cold = ExtrCold + ModCold,
    our_heat = ModHeat + ExtrHeat,
    our_total = ExtrCold + ModCold + ModHeat + ExtrHeat
  )]
  
  # Add city ID
  d_wide[, URAU_CODE := city_id]
  
  return(d_wide)
}

# Parallel processing
n_cores <- min(8, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

message("Processing cities in parallel using ", n_cores, " cores...")

pb <- txtProgressBar(max = length(rds_files), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

our_results_list <- foreach(f = rds_files, 
                            .packages = c("data.table"),
                            .options.snow = opts) %dopar% {
  process_city_file(f)
}

stopCluster(cl)
close(pb)

# Combine results
our_results_list <- our_results_list[!sapply(our_results_list, is.null)]
our_results <- rbindlist(our_results_list, fill = TRUE)

message("\nProcessed ", nrow(our_results), " city-age combinations")

#----- Merge with Masselot baseline

comparison <- merge(
  our_results[, .(URAU_CODE, agegroup, our_cold, our_heat, our_total, 
                  ExtrCold, ModCold, ModHeat, ExtrHeat)],
  masselot_ref,
  by = c("URAU_CODE", "agegroup"),
  all = FALSE  # Inner join - only cities in both datasets
)

message("\nMatched ", nrow(comparison), " city-age combinations with Masselot baseline")
message("Cities in comparison: ", length(unique(comparison$URAU_CODE)))

#----- Calculate differences and metrics

comparison[, `:=`(
  diff_cold = our_cold - masselot_cold,
  diff_heat = our_heat - masselot_heat,
  diff_total = our_total - masselot_total,
  rel_error_cold = abs(our_cold - masselot_cold) / abs(masselot_cold) * 100,
  rel_error_heat = abs(our_heat - masselot_heat) / abs(masselot_heat) * 100,
  rel_error_total = abs(our_total - masselot_total) / abs(masselot_total) * 100,
  # Extreme split percentages
  pct_extrcold = ExtrCold / (ExtrCold + ModCold) * 100,
  pct_extrheat = ExtrHeat / (ModHeat + ExtrHeat) * 100
)]

#----- Summary statistics

message("\n=== VALIDATION METRICS (CITY-LEVEL) ===\n")

message("--- COLD-RELATED MORTALITY ---")
message(sprintf("Mean absolute difference: %.2f deaths/year", 
                mean(abs(comparison$diff_cold), na.rm=TRUE)))
message(sprintf("Median absolute difference: %.2f deaths/year", 
                median(abs(comparison$diff_cold), na.rm=TRUE)))
message(sprintf("Mean relative error: %.2f%%", 
                mean(abs(comparison$rel_error_cold), na.rm=TRUE)))
message(sprintf("Median relative error: %.2f%%", 
                median(abs(comparison$rel_error_cold), na.rm=TRUE)))

cold_concordance <- cor(comparison$our_cold, comparison$masselot_cold, 
                        use = "complete.obs")
message(sprintf("Pearson correlation: %.4f", cold_concordance))

message("\n--- HEAT-RELATED MORTALITY ---")
message(sprintf("Mean absolute difference: %.2f deaths/year", 
                mean(abs(comparison$diff_heat), na.rm=TRUE)))
message(sprintf("Median absolute difference: %.2f deaths/year", 
                median(abs(comparison$diff_heat), na.rm=TRUE)))
message(sprintf("Mean relative error: %.2f%%", 
                mean(abs(comparison$rel_error_heat), na.rm=TRUE)))
message(sprintf("Median relative error: %.2f%%", 
                median(abs(comparison$rel_error_heat), na.rm=TRUE)))

heat_concordance <- cor(comparison$our_heat, comparison$masselot_heat, 
                        use = "complete.obs")
message(sprintf("Pearson correlation: %.4f", heat_concordance))

message("\n--- EXTREME vs MODERATE SPLIT ---")
message(sprintf("ExtrCold: %.1f%% of cold deaths", 
                mean(comparison$pct_extrcold, na.rm=TRUE)))
message(sprintf("ExtrHeat: %.1f%% of heat deaths", 
                mean(comparison$pct_extrheat, na.rm=TRUE)))

#----- By age group

message("\n--- BY AGE GROUP ---")
age_summary <- comparison[, .(
  n = .N,
  r_cold = cor(our_cold, masselot_cold, use="complete.obs"),
  r_heat = cor(our_heat, masselot_heat, use="complete.obs"),
  mean_diff_cold = mean(diff_cold, na.rm=TRUE),
  mean_diff_heat = mean(diff_heat, na.rm=TRUE),
  mean_err_cold_pct = mean(rel_error_cold, na.rm=TRUE),
  mean_err_heat_pct = mean(rel_error_heat, na.rm=TRUE)
), by = agegroup]

print(age_summary)

#----- By country

message("\n--- BY COUNTRY (Top 20) ---")
country_summary <- comparison[, .(
  n_cities = uniqueN(URAU_CODE),
  r_cold = cor(our_cold, masselot_cold, use="complete.obs"),
  r_heat = cor(our_heat, masselot_heat, use="complete.obs"),
  mean_diff_cold = mean(diff_cold, na.rm=TRUE),
  mean_diff_heat = mean(diff_heat, na.rm=TRUE),
  mean_err_cold_pct = mean(rel_error_cold, na.rm=TRUE),
  mean_err_heat_pct = mean(rel_error_heat, na.rm=TRUE)
), by = .(CNTR_CODE, cntr_name)]

country_summary <- country_summary[order(-n_cities)]
print(head(country_summary, 20))

#----- Identify largest discrepancies

message("\n--- TOP 20 DISCREPANCIES (Cold) ---")
top_cold <- comparison[order(-abs(diff_cold))][1:20]
print(top_cold[, .(URAU_CODE, cntr_name, agegroup, our_cold, masselot_cold, 
                   diff_cold, rel_error_cold)])

message("\n--- TOP 20 DISCREPANCIES (Heat) ---")
top_heat <- comparison[order(-abs(diff_heat))][1:20]
print(top_heat[, .(URAU_CODE, cntr_name, agegroup, our_heat, masselot_heat, 
                   diff_heat, rel_error_heat)])

#----- Save full comparison

fwrite(comparison, "results/masselot_validation_city/validation_city_comparison.csv")
message("\nFull comparison saved to: results/masselot_validation_city/validation_city_comparison.csv")

#----- Create validation plots

# Scatter plot: Cold (all city-age combinations)
p1 <- ggplot(comparison, aes(x = masselot_cold, y = our_cold)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_point(aes(color = agegroup), alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
  labs(
    title = "City-Level Cold Mortality: Our 4-Range Sum vs Masselot Baseline",
    subtitle = sprintf("2000-2014 period | n=%d city-age combinations | r=%.3f", 
                      nrow(comparison), cold_concordance),
    x = "Masselot Baseline (deaths/year)",
    y = "Our 4-Range Sum: ExtrCold + ModCold (deaths/year)",
    color = "Age Group"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Scatter plot: Heat (all city-age combinations)
p2 <- ggplot(comparison, aes(x = masselot_heat, y = our_heat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_point(aes(color = agegroup), alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
  labs(
    title = "City-Level Heat Mortality: Our 4-Range Sum vs Masselot Baseline",
    subtitle = sprintf("2000-2014 period | n=%d city-age combinations | r=%.3f", 
                      nrow(comparison), heat_concordance),
    x = "Masselot Baseline (deaths/year)",
    y = "Our 4-Range Sum: ModHeat + ExtrHeat (deaths/year)",
    color = "Age Group"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("results/masselot_validation_city/validation_city_cold_scatter.png", 
       p1, width = 10, height = 7, dpi = 300)
ggsave("results/masselot_validation_city/validation_city_heat_scatter.png", 
       p2, width = 10, height = 7, dpi = 300)

message("\nPlots saved to results/masselot_validation_city/")

# Bland-Altman plots
p3 <- ggplot(comparison, aes(x = (our_cold + masselot_cold)/2, 
                             y = diff_cold)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(aes(color = agegroup), alpha = 0.4, size = 1.5) +
  geom_smooth(se = TRUE, color = "blue") +
  labs(
    title = "Bland-Altman: Cold Mortality Agreement",
    subtitle = "Difference vs Average (city-age level)",
    x = "Average of Our Estimate and Masselot (deaths/year)",
    y = "Difference: Our - Masselot (deaths/year)",
    color = "Age Group"
  ) +
  theme_minimal()

p4 <- ggplot(comparison, aes(x = (our_heat + masselot_heat)/2, 
                             y = diff_heat)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(aes(color = agegroup), alpha = 0.4, size = 1.5) +
  geom_smooth(se = TRUE, color = "blue") +
  labs(
    title = "Bland-Altman: Heat Mortality Agreement",
    subtitle = "Difference vs Average (city-age level)",
    x = "Average of Our Estimate and Masselot (deaths/year)",
    y = "Difference: Our - Masselot (deaths/year)",
    color = "Age Group"
  ) +
  theme_minimal()

ggsave("results/masselot_validation_city/bland_altman_cold.png", 
       p3, width = 10, height = 7, dpi = 300)
ggsave("results/masselot_validation_city/bland_altman_heat.png", 
       p4, width = 10, height = 7, dpi = 300)

#----- Generate summary markdown

summary_md <- sprintf("# City-Level Validation Summary

**Date**: %s  
**Method**: Direct city-by-city comparison (NO aggregation)  
**Baseline period**: 2000-2014 (matching Masselot)

## Dataset Coverage

- **Cities matched**: %d
- **City-age combinations**: %d
- **Age groups**: %s

## Validation Metrics

### Cold-Related Mortality
- **Pearson correlation**: %.4f
- **Mean absolute difference**: %.2f deaths/year
- **Median absolute difference**: %.2f deaths/year
- **Mean relative error**: %.2f%%
- **Median relative error**: %.2f%%

### Heat-Related Mortality
- **Pearson correlation**: %.4f
- **Mean absolute difference**: %.2f deaths/year
- **Median absolute difference**: %.2f deaths/year
- **Mean relative error**: %.2f%%
- **Median relative error**: %.2f%%

### 4-Range Disaggregation Split

- **ExtrCold**: %.1f%% of cold deaths (below 2.5th percentile)
- **ModCold**: %.1f%% of cold deaths (2.5th percentile to MMT)
- **ModHeat**: %.1f%% of heat deaths (MMT to 97.5th percentile)
- **ExtrHeat**: %.1f%% of heat deaths (above 97.5th percentile)

## Interpretation

", 
Sys.Date(),
length(unique(comparison$URAU_CODE)),
nrow(comparison),
paste(unique(comparison$agegroup), collapse=", "),
cold_concordance,
mean(abs(comparison$diff_cold), na.rm=TRUE),
median(abs(comparison$diff_cold), na.rm=TRUE),
mean(abs(comparison$rel_error_cold), na.rm=TRUE),
median(abs(comparison$rel_error_cold), na.rm=TRUE),
heat_concordance,
mean(abs(comparison$diff_heat), na.rm=TRUE),
median(abs(comparison$diff_heat), na.rm=TRUE),
mean(abs(comparison$rel_error_heat), na.rm=TRUE),
median(abs(comparison$rel_error_heat), na.rm=TRUE),
mean(comparison$pct_extrcold, na.rm=TRUE),
100 - mean(comparison$pct_extrcold, na.rm=TRUE),
100 - mean(comparison$pct_extrheat, na.rm=TRUE),
mean(comparison$pct_extrheat, na.rm=TRUE)
)

if (cold_concordance > 0.95 && heat_concordance > 0.85) {
  summary_md <- paste0(summary_md, 
    "✅ **Excellent spatial agreement**: High correlations confirm methodology validity.\n\n")
} else if (cold_concordance > 0.8 && heat_concordance > 0.7) {
  summary_md <- paste0(summary_md,
    "⚠️ **Good spatial agreement**: Correlations acceptable but room for improvement.\n\n")
} else {
  summary_md <- paste0(summary_md,
    "❌ **Poor agreement**: Methodology may need revision.\n\n")
}

summary_md <- paste0(summary_md, "## Files Generated

- `validation_city_comparison.csv` - Full city-age level comparison
- `validation_city_cold_scatter.png` - Scatter plot for cold mortality
- `validation_city_heat_scatter.png` - Scatter plot for heat mortality
- `bland_altman_cold.png` - Agreement analysis for cold
- `bland_altman_heat.png` - Agreement analysis for heat

## Next Steps

")

if (mean(abs(comparison$rel_error_cold), na.rm=TRUE) > 50) {
  summary_md <- paste0(summary_md,
    "1. Investigate why absolute magnitudes differ despite good correlation\n",
    "2. Check if Masselot uses different aggregation period (annual vs cumulative)\n",
    "3. Verify MMT centering and DLNM specification match Masselot exactly\n",
    "4. Compare temperature data sources (GCM vs observed)\n")
} else {
  summary_md <- paste0(summary_md,
    "✅ Validation successful - proceed with projections and LE/LI decomposition\n")
}

writeLines(summary_md, "results/masselot_validation_city/CITY_VALIDATION_SUMMARY.md")

message("\n=== VALIDATION COMPLETE ===")
message("See results/masselot_validation_city/ for detailed outputs")
