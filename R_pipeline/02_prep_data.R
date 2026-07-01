################################################################################
# 
# Reproduction and extension of Masselot et al. 2025
# 
# Part 2: Prepare city metadata and baseline demographics
#
################################################################################

# Load parameters
if(file.exists("01_initialize.R")) source("01_initialize.R") else source("R_pipeline/01_initialize.R")

#-------------------------------------------------------------------------------
# 1. LOAD CITY METADATA
#-------------------------------------------------------------------------------

# Load city results (contains MMT, age groups, and baseline deaths)
city_meta <- fread(path_city_res)

# Select relevant columns and rename for convenience
city_data <- city_meta %>%
  select(
    URAU_CODE, 
    LABEL, 
    CNTR_CODE, 
    cntr_name, 
    region, 
    agegroup, 
    pop_baseline = agepop, 
    death_baseline = death, 
    mmt
  ) %>%
  mutate(URAU_CODE = as.character(URAU_CODE))

# Get list of unique cities
cities <- unique(city_data$URAU_CODE)
ncities <- length(cities)

#-------------------------------------------------------------------------------
# 2. LOAD OBSERVATIONS FOR BIAS CORRECTION
#-------------------------------------------------------------------------------

# Load ERA5-Land series for training bias correction
obs_data <- fread(path_obs)
obs_data[, date := as.Date(date)]

# Filter for relevant years (hist_years)
obs_data <- obs_data[year(date) %in% hist_years]

#-------------------------------------------------------------------------------
# 3. PREPARE COEFFICIENTS
#-------------------------------------------------------------------------------

# We have coef_simu.csv with 500 simulations per (city, agegroup)
# This will be loaded during the loop to save memory, or pre-loaded if feasible.
# For now, let's just check the number of rows.
# n_coef_simu <- fread(path_coef_simu, select = "URAU_CODE") %>% nrow()
# message("Found ", n_coef_simu, " rows in coef_simu.csv")

#-------------------------------------------------------------------------------
# 4. SAVE PREPARED DATA
#-------------------------------------------------------------------------------

save(city_data, cities, ncities, obs_data, file = "data/prep_data.RData")

message("Data preparation complete. ", ncities, " cities identified.")
