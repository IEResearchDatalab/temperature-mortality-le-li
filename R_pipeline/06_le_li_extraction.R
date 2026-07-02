################################################################################
# 
# Contrasting future heat and cold-related mortality under climate change
#
# Extraction of ANs for LE/LI Decomposition (Aburto methodology)
# This script aggregates city-level results into country-level ANs for 
# age groups 65+ in the baseline (2010s) and future (2090s).
#
################################################################################

library(data.table)
library(doSNOW)
library(foreach)

#----- Parameters

n_cores <- 8
input_dir <- "temp_results"
output_file <- "data/le_li_input_ans.csv"

# Time periods
baseline_decade <- 2010
future_decade <- 2090

# Age groups of interest
age_groups_65plus <- c("65-74", "75-84", "85+")

#----- Load Metadata

city_meta <- fread("data/city_results.csv")
city_meta <- city_meta[, .(URAU_CODE, cntr_name, agegroup, pop = agepop, death_baseline = death)]

#----- Processing Function

process_city <- function(f, city_meta, baseline_decade, future_decade, groups) {
  city_id <- gsub(".rds", "", basename(f))
  d <- try(readRDS(f), silent = TRUE)
  if(inherits(d, "try-error")) return(NULL)
  
  # Filter for decades and age groups
  d[, decade := (year %/% 10) * 10]
  d_sub <- d[decade %in% c(baseline_decade, future_decade) & agegroup %in% groups]
  
  if(nrow(d_sub) == 0) return(NULL)
  
  # Group ANs by decade, range, ssp, gcm, agegroup (average over simulations)
  # Aggregating by sim first to get distribution, but Simon usually wants the mean.
  # Let's keep simulations for now to allow for CI calculation if needed, 
  # but maybe just the mean is enough for a prompt answer.
  # Actually, the lifetable code can handle simulations. 
  # But 500 simulations x 854 cities is too much data to keep in memory easily.
  # Let's take the mean AN over simulations per city/decade/range/ssp/gcm/agegroup.
  
  d_agg <- d_sub[, .(an = sum(an) / 10), by = .(decade, range, ssp, gcm, agegroup)]
  
  # Add city meta
  city_info <- city_meta[URAU_CODE == city_id]
  d_agg <- merge(d_agg, city_info, by = "agegroup")
  
  return(d_agg)
}

#----- Execution

rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
message("Extracting data for ", length(rds_files), " cities...")

cl <- makeCluster(n_cores)
registerDoSNOW(cl)

clusterExport(cl, c("process_city", "age_groups_65plus", "baseline_decade", "future_decade"))

all_results <- foreach(f = rds_files, .combine = function(a,b) rbindlist(list(a,b)), .packages = c("data.table")) %dopar% {
  process_city(f, city_meta, baseline_decade, future_decade, age_groups_65plus)
}

stopCluster(cl)

#----- Aggregate to Country Level

country_ans <- all_results[, .(
  an = sum(an),
  pop = sum(pop),
  death_baseline = sum(death_baseline)
), by = .(cntr_name, decade, ssp, gcm, range, agegroup)]

# Save results
fwrite(country_ans, output_file)
message("Saved country-level ANs to ", output_file)
