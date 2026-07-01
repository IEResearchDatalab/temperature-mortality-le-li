################################################################################
# 
# Reproduction and extension of Masselot et al. 2025
# 
# Part 1: Packages and analysis parameters
#
################################################################################

#------------------------
# LOAD THE PACKAGES
#------------------------

#----- Data management
library(dplyr)
library(tidyr)
library(data.table)
library(dtplyr)
library(arrow)
library(stringr)
library(sf)
library(openxlsx)

#----- Statistical analysis
library(dlnm)
library(splines)
library(doParallel)
library(doSNOW)
library(mixmeta)
library(collapse)

#----- Custom functions
if(file.exists("functions/isimip3.R")){
  source("functions/isimip3.R")
  source("functions/impact.R")
} else {
  source("R_pipeline/functions/isimip3.R")
  source("R_pipeline/functions/impact.R")
}

#------------------------
# PARAMETERS
#------------------------

#----- Computation
nsim <- 500
ncores <- pmax(detectCores() - 1, 1) |> pmin(32) # Allow more cores if available
grpsize <- 10
tdir <- "temp_results"
if(!dir.exists(tdir)) dir.create(tdir)

#----- Analysis
agebreaks <- c(20, 45, 65, 75, 85, Inf)
agelabs <- gsub("-Inf", "+", 
  paste(agebreaks[-length(agebreaks)], agebreaks[-1] - 1, sep = "-"))

# Specification of the exposure-response function
varfun <- "bs"
vardegree <- 2
varper <- c(10, 75, 90)
vardf <- vardegree + length(varper)

# Temperature percentiles for disaggregation
# Original: Cold/Heat split at MMT
# This project: 
# 1. Extreme Cold (T < P2.5)
# 2. Moderate Cold (P2.5 <= T < MMT)
# 3. Moderate Heat (MMT <= T <= P97.5)
# 4. Extreme Heat (T > P97.5)
range_thresholds <- c(2.5, 97.5) 
range_labels <- c("Extreme Cold", "Moderate Cold", "Moderate Heat", "Extreme Heat")

# Scenarios
scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")
gcms <- c("GFDL_ESM4", "IPSL_CM6A_LR", "MPI_ESM1_2_HR", "MRI_ESM2_0")

# Year ranges
hist_years <- 1990:2014
fut_years <- 2015:2099
pred_years <- 2015:2099

# Paths (assuming run from project root)
dir_data <- "data/"
dir_ref <- "references/2025-masselot-zenodo/"
if(!dir.exists(dir_data)) {
  # Try relative if run from R_pipeline/
  dir_data <- "../data/"
  dir_ref <- "../references/2025-masselot-zenodo/"
}
path_tmean <- paste0(dir_data, "tmeanproj.gz.parquet")
path_obs <- paste0(dir_ref, "additional_data/era5series.csv")
path_coefs <- paste0(dir_data, "coefs.csv")
path_vcov <- paste0(dir_data, "vcov.csv")
path_coef_simu <- paste0(dir_data, "coef_simu.csv")
path_city_res <- paste0(dir_data, "city_results.csv")
