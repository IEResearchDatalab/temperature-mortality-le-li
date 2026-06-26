################################################################################
#
# Validation: Step 00 — Setup and parameters
#   Verify global constants, paths, and package loading.
#
################################################################################

library(ggplot2)

validation_dir <- "../results/validation_plots"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Validation Step 00: Setup parameters\n")
cat("========================================\n\n")

#--- Parameter dump
cat(sprintf("  Paths:\n"))
cat(sprintf("    data:   %s\n", path_data))
cat(sprintf("    output: %s\n", path_out))
cat(sprintf("  Cities: %s\n", if (is.null(city_subset)) "ALL (854)" else city_subset))
cat(sprintf("  SSPs:   %s\n", paste(names(ssp_labels), ssp_labels, sep = "=", collapse = ", ")))
cat(sprintf("  Years:  hist %d-%d, proj %d-%d\n",
    hist_year_min, hist_year_max, proj_year_min, proj_year_max))
cat(sprintf("  Age groups: %s\n", paste(agelabs, collapse = ", ")))
cat(sprintf("  Single ages: %d to %d\n", min(single_ages), max(single_ages)))
cat(sprintf("  Basis: %s degree=%d knots at %s pct\n",
    varfun, vardegree, paste(varper, collapse = ", ")))
cat(sprintf("  Cold extreme pct: %.1f  Heat extreme pct: %.1f\n",
    cold_extreme_pct * 100, heat_extreme_pct * 100))
cat(sprintf("  Simulations: %d  Cores: %d\n", nsim, ncores))
cat(sprintf("  Adaptation levels: %s\n", paste(adapt_levels, collapse = ", ")))

#--- Check packages
required <- c("data.table", "arrow", "dplyr", "dlnm", "splines", "MASS",
              "ungroup", "ggplot2", "foreach", "doParallel", "eurostat")
missing <- setdiff(required, loadedNamespaces())
if (length(missing) > 0) {
  cat(sprintf("  WARNING: packages not loaded: %s\n", paste(missing, collapse = ", ")))
} else {
  cat("  All required packages loaded OK\n")
}

#--- Quick data file check
data_files <- c("coefs.csv", "coef_simu.csv", "city_results.csv", "tmeanproj.gz.parquet")
for (f in data_files) {
  exists <- file.exists(file.path(path_data, f))
  cat(sprintf("  Data file %s: %s\n", f, if (exists) "OK" else "MISSING"))
}

cat("\n--- Step 00 validation complete ---\n")