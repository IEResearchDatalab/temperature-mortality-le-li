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

#----- Libraries

library(data.table)
library(dplyr)
library(arrow)

#----- Simulation parameters

# Number of Monte Carlo simulations and cores for parallelization
nsim <- 500
n_cores <- 8

#----- Study periods and scenarios

# Reference period for historical bias correction
hist_years <- 1990:2014

# Scenarios mapping (SSP 126, 245, 585)
scenarios <- c(1, 2, 3) 

# General Circulation Models (GCMs) from CMIP6
gcms <- c("GFDL_ESM4", "IPSL_CM6A_LR", "MPI_ESM1_2_HR", "MRI_ESM2_0")

#----- Exposure-response parameters (Masselot 2025)

# Internal knots for the natural cubic spline cross-basis
knots_percentiles <- c(10, 75, 90)

#----- Paths

# Data and coefficient directory paths
path_tmean <- "data/tmeanproj.gz.parquet"  
path_coef_simu <- "data/coef_simu.csv"
path_coefs <- "data/coefs.csv"
path_vcov <- "data/vcov.csv"
path_city <- "data/city_results.csv"

#----- Final initialization

# Ensure reproducibility of the random draws
set.seed(1308)

message("Environment initialized with ", n_cores, " cores and nsim=", nsim)
