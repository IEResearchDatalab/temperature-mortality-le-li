################################################################################
# 
# Contrasting future heat and cold-related mortality under climate change, 
# demographic and adaptation scenarios in 854 European cities
#
# R Code Part 5: Generate Tables (Masselot 2023 & Lloyd 2024 Reproductions)
#
# Adapted from Pierre Masselot & Antonio Gasparrini
#
################################################################################

#----- Libraries and environment

library(data.table)
library(dplyr)
library(foreach)
library(doSNOW)

# Source global parameters
source("R_pipeline/01_initialize.R")

message("\n[5/5] Generating Summary Tables...")

#----- Prepare metadata and weights

# 1. Load City Metadata (Country, baseline population, deaths)
city_meta <- fread("data/city_results.csv")
city_meta <- city_meta[, .(
  URAU_CODE, 
  cntr_name, 
  agegroup, 
  pop = agepop, 
  deaths_baseline = death
)]

# 2. ESP 2013 Weights for Age-Standardization (Relative to 20+ population)
esp_weights <- data.table(
  agegroup = c("20-44", "45-64", "65-74", "75-84", "85+"),
  weight = c(31500, 26500, 10500, 6500, 2500)
)
esp_weights[, weight := weight / sum(weight)]

#----- Perform granular country-level aggregation

# Identify all city-specific result files
rds_files <- list.files("temp_results", pattern = "\\.rds$", full.names = TRUE)

message("Aggregating country-level data from ", length(rds_files), " cities...")

# Setup parallel processing
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# Function to extract and aggregate results per city
process_city_tables <- function(f, city_meta, esp_weights) {
  city_id <- gsub(".rds", "", basename(f))
  d <- try(readRDS(f), silent = TRUE)
  if(inherits(d, "try-error")) return(NULL)
  
  # Filter metadata for this city
  meta_sub <- city_meta[URAU_CODE == city_id]
  setkey(meta_sub, agegroup)
  
  # Aggregate by decade/agegroup/range/ssp/gcm/sim
  d[, decade := (year %/% 10) * 10]
  d_agg <- d[, .(an = sum(an)), by = .(decade, range, ssp, gcm, agegroup, sim)]
  
  # Join country-level metadata
  d_agg[, country := meta_sub$cntr_name[1]]
  d_agg[, pop_baseline := meta_sub[agegroup == .BY$agegroup]$pop, by = agegroup]
  d_agg[, deaths_baseline := meta_sub[agegroup == .BY$agegroup]$deaths_baseline, by = agegroup]
  
  return(d_agg)
}

# Export environment to workers
clusterExport(cl, c("process_city_tables", "city_meta", "esp_weights"))

# Execute parallel aggregation
country_results <- foreach(f = rds_files, .combine = function(a,b) rbindlist(list(a,b)), .packages = c("data.table")) %dopar% {
  process_city_tables(f, city_meta, esp_weights)
}
stopCluster(cl)

#----- Process summary tables

dir.create("tables", showWarnings = FALSE)

# --- TABLE 1: Country-level annual excess deaths (Masselot 2023 Table S6) ---

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
table_s2_raw <- country_results[, .(
  an_sum = sum(an),
  pop_sum = sum(pop_baseline)
), by = .(agegroup, range, decade, ssp, gcm, sim)]

table_s2_ar <- table_s2_raw[, .(
  AR = mean(an_sum / (pop_sum + 1e-6)) * 100000
), by = .(agegroup, range, decade, ssp)]

# Compare 2020s vs 2090s for SSP5
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
