################################################################################
# 
# Reproduction and extension of Masselot et al. 2025
# 
# Part 4: Aggregate results and generate summary tables
#
################################################################################

# Load parameters
if(file.exists("01_initialize.R")) source("01_initialize.R") else source("R_pipeline/01_initialize.R")

# Get list of city result files
res_files <- list.files(tdir, pattern = "\\.rds$", full.names = TRUE)

if(length(res_files) == 0) stop("No results found in ", tdir)

message("Aggregating results from ", length(res_files), " cities...")

# Function to process all cities for a specific scenario
process_all <- function() {
  all_res <- list()
  for (f in res_files) {
    all_res[[length(all_res) + 1]] <- readRDS(f)
  }
  rbindlist(all_res)
}

# Load all
full_data <- process_all()

#-------------------------------------------------------------------------------
# 1. TOTAL EUROPE RESULTS
#-------------------------------------------------------------------------------

# Sum across all cities
europe_sims <- full_data[, .(an = sum(an)), by = .(decade, range, sc, agegroup, gcm, sim)]

# Average across GCMs
europe_ens <- europe_sims[, .(an = mean(an)), by = .(decade, range, sc, agegroup, sim)]

# Compute quantiles across simulations
europe_final <- europe_ens[, as.list(quantile(an, probs = c(0.025, 0.5, 0.975))), 
                           by = .(decade, range, sc, agegroup)]

setnames(europe_final, c("2.5%", "50%", "97.5%"), c("low", "est", "high"))

# Save Europe total
fwrite(europe_final, "data/europe_decadal_results.csv")

#-------------------------------------------------------------------------------
# 2. DISAGGREGATED SUMMARY TABLE (Modern Period vs Late Century)
#-------------------------------------------------------------------------------

summary_table <- europe_final[decade %in% c(2020, 2090)] %>%
  arrange(sc, agegroup, range, decade)

fwrite(summary_table, "data/summary_impact_table.csv")

message("Aggregation complete. Final results in data/europe_decadal_results.csv")
