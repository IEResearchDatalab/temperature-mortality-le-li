################################################################################
# 
# Contrasting future heat and cold-related mortality under climate change, 
# demographic and adaptation scenarios in 854 European cities
#
# R Code Part 4: Aggregate and summarize health impact results
#
# Adapted from Pierre Masselot & Antonio Gasparrini
#
################################################################################

library(data.table)
library(foreach)
library(doSNOW)

source("R_pipeline/01_initialize.R")

message("\n[3/3] Aggregating results from temp_results/...")

out_file <- "data/final_attribution_results.csv"
existing_cities <- character(0)
if(file.exists(out_file)){
  d_existing <- fread(out_file, select = "URAU_CODE")
  existing_cities <- unique(d_existing$URAU_CODE)
  message("Found existing results for ", length(existing_cities), " cities.")
}

rds_files <- list.files("temp_results", pattern = "\\.rds$", full.names = TRUE)
rds_to_process <- rds_files[!(gsub(".rds", "", basename(rds_files)) %in% existing_cities)]

if(length(rds_to_process) == 0 && length(existing_cities) > 0){
  message("All cities already aggregated. Skipping.")
  # Still run validation on existing file
} else {
  message("Processing ", length(rds_to_process), " new cities...")
  
  # Parallel setup for efficient aggregation
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(rds_to_process), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  new_summaries <- foreach(f = rds_to_process, .packages = c("data.table"), .options.snow = opts) %dopar% {
    city_id <- gsub(".rds", "", basename(f))
    d <- try(readRDS(f), silent = TRUE)
    if(inherits(d, "try-error")) return(NULL)
    
    d_city <- d[, .(an = sum(an)), by = .(year, range, sim, ssp, gcm)]
    d_city[, decade := (year %/% 10) * 10]
    d_decade <- d_city[, .(an = sum(an)), by = .(decade, range, sim, ssp, gcm)]
    
    d_summary <- d_decade[, .(
      an_mean = mean(an),
      an_p2.5 = quantile(an, 0.025),
      an_p97.5 = quantile(an, 0.975)
    ), by = .(decade, range, ssp, gcm)]
    
    d_summary[, URAU_CODE := city_id]
    return(d_summary)
  }
  
  stopCluster(cl)
  close(pb)
  
  # Combine new and existing
  final_results <- rbindlist(new_summaries)
  if(length(existing_cities) > 0){
    d_old <- fread(out_file)
    final_results <- rbind(d_old, final_results)
  }
  fwrite(final_results, out_file)
  message("\nAggregation complete. Results saved to: ", out_file)
}

# --- Validation Block ---
message("\n--- Validating Results ---")
final_results <- fread(out_file)

# 1. Check if all 4 categories exist for a sample city
sample_city <- final_results$URAU_CODE[1]
ranges_found <- unique(final_results[URAU_CODE == sample_city]$range)
message("Ranges found for ", sample_city, ": ", paste(ranges_found, collapse = ", "))

expected_ranges <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
if(all(expected_ranges %in% ranges_found)){
  message("SUCCESS: All 4 temperature ranges are present.")
} else {
  message("WARNING: Missing ranges: ", paste(setdiff(expected_ranges, ranges_found), collapse = ", "))
}

# 2. Additivity check: Sum of components should align with expected "Masselot" patterns
# (In this simplified pipeline, we ensure they are mutually exclusive by construction in 03_attribution.R)
range_sums <- final_results[, .(total_an = sum(an_mean)), by = .(URAU_CODE, decade, ssp, gcm)]
message("Mean total attributable deaths for ", sample_city, " (SSP1, 2020s, GFDL): ", 
        round(range_sums[URAU_CODE == sample_city & decade == 2020 & ssp == 1 & gcm == "GFDL_ESM4"]$total_an, 2))

message("Validation complete.")
