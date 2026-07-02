################################################################################
# 
# Contrasting future heat and cold-related mortality under climate change, 
# demographic and adaptation scenarios in 854 European cities
#
# R Code Part 2: Prepare city metadata and historical temperature data
#
# Adapted from Pierre Masselot & Antonio Gasparrini
#
################################################################################

#----- Libraries and environment

library(data.table)
library(dplyr)
library(arrow)

# Source global parameters
source("R_pipeline/01_initialize.R")

message("\n[1/3] Preparing data for 854 cities...")

#----- Load data

# Load city-specific thresholds and historical summaries
city_results <- fread(path_city)
coefs <- fread(path_coefs)

#----- Process metadata

# Filter for the relevant cities (854 cities in the main study)
cities <- unique(city_results$URAU_CODE)

# Calculate age-specific MMTs and temperature range thresholds (p2.5 and p97.5)
thresholds <- city_results[, .(
  URAU_CODE, 
  agegroup, 
  mmt, 
  p2_5, 
  p97_5, 
  death
)]

#----- Load historical temperatures

# Load historical daily temperature data for bias correction (ISIMIP3 mapping)
obs_ds <- open_dataset("data/tmean_obs") 
obs_data <- obs_ds %>% collect() %>% as.data.table()

#----- Save for the loop

# Save prepared objects for the Part 3 loop
save(cities, thresholds, obs_data, file = "data/prep_data.RData")

message("Preparation complete. ", length(cities), " cities ready for simulation.")
