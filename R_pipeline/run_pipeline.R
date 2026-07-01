################################################################################
# 
# Master script to run the full pipeline
#
################################################################################

message(">>> Starting Pipeline: Reproduction of Masselot et al. 2025 <<<")

# Step 2: Prep data (Step 1 is sourced within others)
message("\n[1/3] Preparing data...")
source("02_prep_data.R")

# Step 3: Attribution loop (Parallel)
message("\n[2/3] Running attribution simulations (this may take a while)...")
# Note: You can adjust ncores in 01_initialize.R
source("03_attribution.R")

# Step 4: Aggregate results
message("\n[3/3] Aggregating results...")
source("04_aggregate_results.R")

message("\n>>> Pipeline finished successfully! <<<")
message("Check 'data/europe_decadal_results.csv' for the final output.")
