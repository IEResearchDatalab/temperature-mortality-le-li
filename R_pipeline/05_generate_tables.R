################################################################################
# 
# Contrasting future heat and cold-related mortality under climate change, 
# demographic and adaptation scenarios in 854 European cities
#
# R Code Part 5: Generate Tables (Masselot 2023 & Lloyd 2024 Reproductions)
#
################################################################################

library(data.table)
library(dplyr)
library(foreach)
library(doSNOW)

source("R_pipeline/01_initialize.R")

message("\n[5/5] Generating Summary Tables...")

# 1. Load City Metadata (Country, baseline population, deaths)
city_meta <- fread("data/city_results.csv")
# Keep unique city-age combinations
city_meta <- city_meta[, .(
  URAU_CODE, 
  cntr_name, 
  agegroup, 
  pop = agepop, 
  deaths_baseline = death
)]

# 2. ESP 2013 Weights for Age-Standardization (Relative to 20+ population)
# Weights per 100,000 for relevant groups:
# 20-44: 5*6500 + 4*6000? No, let's use the standard classes:
# 20-24: 6k, 25-29: 6k, 30-34: 6.5k, 35-39: 6.5k, 40-44: 6.5k -> 31.5k
# 45-49: 7k, 50-54: 7k, 55-59: 6.5k, 60-64: 6k -> 26.5k
# 65-69: 5.5k, 70-74: 5k -> 10.5k
# 75-79: 4k, 80-84: 2.5k -> 6.5k
# 85+: 2.5k
esp_weights <- data.table(
  agegroup = c("20-44", "45-64", "65-74", "75-84", "85+"),
  weight = c(31500, 26500, 10500, 6500, 2500)
)
esp_weights[, weight := weight / sum(weight)]

# 3. Load Granular Results (from temp_results instead of the aggregated csv)
# This is slow for all cities, so we aggregate per country iteratively.
rds_files <- list.files("temp_results", pattern = "\\.rds$", full.names = TRUE)

message("Aggregating country-level data from ", length(rds_files), " cities...")

# Setup parallel processing
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# Function to process one city and return country/decade/age/range summary
process_city_tables <- function(f, city_meta, esp_weights) {
  city_id <- gsub(".rds", "", basename(f))
  d <- try(readRDS(f), silent = TRUE)
  if(inherits(d, "try-error")) return(NULL)
  
  # Join with meta for country info
  meta_sub <- city_meta[URAU_CODE == city_id]
  setkey(meta_sub, agegroup)
  
  # Aggregate by decade/agegroup/range/ssp/gcm/sim
  d[, decade := (year %/% 10) * 10]
  
  d_agg <- d[, .(an = sum(an)), by = .(decade, range, ssp, gcm, agegroup, sim)]
  
  # Join Country
  d_agg[, country := meta_sub$cntr_name[1]]
  d_agg[, pop_baseline := meta_sub[agegroup == .BY$agegroup]$pop, by = agegroup]
  d_agg[, deaths_baseline := meta_sub[agegroup == .BY$agegroup]$deaths_baseline, by = agegroup]
  
  return(d_agg)
}

# Export functions and data to workers
clusterExport(cl, c("process_city_tables", "city_meta", "esp_weights"))

# We process in chunks to avoid memory overhead of return list
country_results <- foreach(f = rds_files, .combine = function(a,b) rbindlist(list(a,b)), .packages = c("data.table")) %dopar% {
  process_city_tables(f, city_meta, esp_weights)
}
stopCluster(cl)

# --- SAVE INTERMEDIATE COUNTRY DATA ---
# This is large but useful for the 3 tables
dir.create("tables", showWarnings = FALSE)

# --- TABLE 1: Country-level annual excess deaths (Masselot 2023 Table S6) ---
# Sum across ages and GCMS for a given SSP and Decade
# We use mean across GCMS and then summary stats across SIMS
# Or better: mean across SIMS per GCM, then mean of GCMS? 
# Usually, we treat (GCM x SIM) as the uncertainty pool.

t1_raw <- country_results[, .(
  an_total = sum(an),
  pop_total = sum(pop_baseline),
  deaths_total = sum(deaths_baseline)
), by = .(country, decade, ssp, range, gcm, sim)]

# Sum Cold and Heat separately
t1_raw[, type := ifelse(grepl("Cold", range), "Cold", "Heat")]
t1_combined <- t1_raw[, .(
  an = sum(an_total),
  pop = mean(pop_total),
  deaths_base = mean(deaths_total)
), by = .(country, decade, ssp, type, gcm, sim)]

# Summary stats across uncertainty pool (GCM x SIM)
table_s6 <- t1_combined[, .(
  AN = mean(an),
  AN_low = quantile(an, 0.025),
  AN_hi = quantile(an, 0.975),
  AF = mean(an / (deaths_base + 1e-6)) * 100,
  AF_low = quantile(an / (deaths_base + 1e-6), 0.025) * 100,
  AF_hi = quantile(an / (deaths_base + 1e-6), 0.975) * 100,
  Rate = mean(an / (pop + 1e-6)) * 100000,
  Rate_low = quantile(an / (pop + 1e-6), 0.025) * 100000,
  Rate_hi = quantile(an / (pop + 1e-6), 0.975) * 100000
), by = .(country, decade, ssp, type)]

fwrite(table_s6, "tables/Table_Masselot2023_S6_Country.csv")


# --- TABLE 2: Average annual number of deaths by Temperature Range (Lloyd 2024 Table S1) ---
# Group by Range and Decade/SSP, sum across countries
table_s1 <- t1_raw[, .(
  an_total = sum(an_total)
), by = .(decade, ssp, range, gcm, sim)]

table_s1_summary <- table_s1[, .(
  AN = mean(an_total),
  low = quantile(an_total, 0.025),
  hi = quantile(an_total, 0.975)
), by = .(decade, ssp, range)]

fwrite(table_s1_summary, "tables/Table_Lloyd2024_S1_Ranges.csv")


# --- TABLE 3: Percent change in Absolute Risk (Lloyd 2024 Table S2) ---
# Absolute Risk = AN / Pop
# We compare Decade 2020 vs 2090 for SSP585
table_s2_raw <- country_results[, .(
  an_sum = sum(an),
  pop_sum = sum(pop_baseline)
), by = .(agegroup, range, decade, ssp, gcm, sim)]

table_s2_ar <- table_s2_raw[, .(
  AR = mean(an_sum / (pop_sum + 1e-6)) * 100000
), by = .(agegroup, range, decade, ssp)]

# Compare 2020s vs 2090s
ar_2020 <- table_s2_ar[decade == 2020 & ssp == 5]
ar_2090 <- table_s2_ar[decade == 2090 & ssp == 5]

setnames(ar_2020, "AR", "AR_2020")
setnames(ar_2090, "AR", "AR_2090")

table_s2_change <- merge(ar_2020[, .(agegroup, range, AR_2020)], 
                         ar_2090[, .(agegroup, range, AR_2090)], 
                         by = c("agegroup", "range"))

table_s2_change[, pct_change := (AR_2090 - AR_2020) / (AR_2020 + 1e-6) * 100]

fwrite(table_s2_change, "tables/Table_Lloyd2024_S2_Change.csv")

message("\nTable generation complete. Files saved in tables/ folder.")
