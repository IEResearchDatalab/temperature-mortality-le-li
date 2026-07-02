################################################################################
# 
# Contrasting future heat and cold-related mortality under climate change, 
# demographic and adaptation scenarios in 854 European cities
#
# R Code Part 1: Global parameters and environment setup
#
# Adapted from Pierre Masselot & Antonio Gasparrini
#
################################################################################

# Libraries
library(data.table)
library(dplyr)
library(arrow)

# Simulation parameters
nsim <- 500
n_cores <- 8

# Study periods and scenarios
hist_years <- 1990:2014
scenarios <- c(1, 2, 3) # Mapping to SSP 126, 245, 585
gcms <- c("GFDL_ESM4", "IPSL_CM6A_LR", "MPI_ESM1_2_HR", "MRI_ESM2_0")

# Exposure-response parameters (Masselot 2025)
knots_percentiles <- c(10, 75, 90)

# Paths
path_tmean <- "data/tmeanproj.gz.parquet"  # Daily mean temperature parquet
path_coef_simu <- "data/coef_simu.csv"
path_coefs <- "data/coefs.csv"
path_vcov <- "data/vcov.csv"
path_city <- "data/city_results.csv"

# Ensure reproducibility
set.seed(1308)

message("Environment initialized with ", n_cores, " cores and nsim=", nsim)
