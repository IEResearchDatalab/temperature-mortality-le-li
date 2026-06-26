################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#   a multi-city, multi-GCM, multi-SSP health impact assessment
#
# Master pipeline script
#   Sources all pipeline scripts interleaved with validation diagnostics.
#   Each step produces diagnostic plots before proceeding to the next.
#
# Pipeline order:
#   00  Packages and parameters
#   01  Raw input preparation → validation 01
#   02  Climate-only AN computation → validation 02
#   03  Mortality backbone → validation 03
#   04  City mortality schedules → validation 04
#   05  AN disaggregation → validation 05
#   06  Cause-specific bookkeeping → validation 06
#   07  Life tables for ages 65+ → validation 07
#   08  Decomposition → validation 08
#   09  Aggregation to country and Europe
#   10  Absolute risk analysis
#
# Usage: Rscript 00_RunAll.R
#
################################################################################

# Load setup (packages, paths, parameters)
source("00_Packages_Parameters.R")

#------------------------------
# Step 01: Raw input preparation
#------------------------------
cat("\n========== Step 01: Prep ==========\n")
source("01_PrepData.R")
cat("\n--- Validating Step 01 ---\n")
source("../validation/01_validate_prep.R")

#------------------------------
# Step 02: AN computation
#------------------------------
cat("\n========== Step 02: AN computation ==========\n")
source("02_ComputeAN.R")
cat("\n--- Validating Step 02 ---\n")
source("../validation/02_validate_an.R")

#------------------------------
# Step 03: Mortality backbone
#------------------------------
cat("\n========== Step 03: Mortality backbone ==========\n")
source("03_MortalityBackbone.R")
cat("\n--- Validating Step 03 ---\n")
source("../validation/03_validate_backbone.R")

#------------------------------
# Step 04: City mortality schedules
#------------------------------
cat("\n========== Step 04: City schedules ==========\n")
source("04_CityMortalitySchedules.R")
cat("\n--- Validating Step 04 ---\n")
source("../validation/04_validate_city_schedules.R")

#------------------------------
# Step 05: AN disaggregation
#------------------------------
cat("\n========== Step 05: Disaggregation ==========\n")
source("05_Disaggregate.R")
cat("\n--- Validating Step 05 ---\n")
source("../validation/05_validate_disaggregate.R")

#------------------------------
# Step 06: Bookkeeping
#------------------------------
cat("\n========== Step 06: Bookkeeping ==========\n")
source("06_AnalysisDataset.R")
cat("\n--- Validating Step 06 ---\n")
source("../validation/06_validate_bookkeeping.R")

#------------------------------
# Step 07: Life tables
#------------------------------
cat("\n========== Step 07: Life tables ==========\n")
source("07_LifeTables.R")
cat("\n--- Validating Step 07 ---\n")
source("../validation/07_validate_lifetables.R")

#------------------------------
# Step 08: Decomposition
#------------------------------
cat("\n========== Step 08: Decomposition ==========\n")
source("08_Decomposition.R")
cat("\n--- Validating Step 08 ---\n")
source("../validation/08_validate_decomposition.R")

#------------------------------
# Step 09: Aggregation
#------------------------------
cat("\n========== Step 09: Aggregation ==========\n")
source("09_Aggregation.R")

#------------------------------
# Step 10: Absolute risk
#------------------------------
cat("\n========== Step 10: Absolute risk ==========\n")
source("10_AbsoluteRisk.R")

#------------------------------
# Done
#------------------------------
cat("\n========================================\n")
cat("=== Pipeline complete ===\n")
cat(sprintf("Output directory: %s\n", normalizePath(path_out)))
cat(sprintf("Number of cities: %d\n", n_cities))
cat(sprintf("Validation plots: %s\n",
  normalizePath("../results/validation_plots")))
cat("========================================\n")