################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#   a multi-city, multi-GCM, multi-SSP health impact assessment
#
# Pipeline: 00_Packages_Parameters.R
#   Loading packages and analysis parameters
#
# Based on the framework of:
#   Masselot et al. (2023) Lancet Planetary Health
#   Gasparrini & Leone (2014) BMC Medical Research Methodology
#
################################################################################

#----------------------
# Necessary packages
#----------------------

# Data management
library(data.table); library(arrow); library(dplyr)
library(doParallel); library(doSNOW); library(doRNG)

# Analysis
library(dlnm); library(splines); library(MASS)
library(ungroup)

# Plotting
library(ggplot2); library(scales)

#----------------------
# Parameters
#----------------------

# Paths
path_data <- "../data"
path_out <- "../results/production"

dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
seed_main <- 12345
set.seed(seed_main)

# City selection (all 854 European cities)
# If NULL, all cities in coefs.csv are processed
# If a single city code (e.g. "ES001C" for Madrid), only that city is processed
city_subset <- "ES001C"  # Madrid demo; set to NULL for all 854 cities

# Scenario selection
ssp_keep <- c("hist", "1", "2", "3")
ssp_labels <- c(
  "hist" = "Historical",
  "1"  = "SSP1-2.6",
  "2"  = "SSP3-7.0",
  "3"  = "SSP5-8.5"
)

# Time windows
hist_year_min <- 1990
hist_year_max <- 2019
proj_year_min <- 2020
proj_year_max <- 2100

# Adaptation levels (0 = no adaptation, 0.5 = 50%, 0.9 = 90%)
adapt_levels <- c(0)

# Temperature range thresholds
cold_extreme_pct <- 0.025
heat_extreme_pct <- 0.975

# Age groups (Masselot 2023)
agebreaks <- c(20, 45, 65, 75, 85)
agelabs <- c("20-44", "45-64", "65-74", "75-84", "85+")
age_midpoints <- c(32, 54.5, 69.5, 79.5, 92.5)
single_ages <- 20:100
nage <- length(single_ages)
pclm_nlast <- 16  # ages 85-100

# Basis parameters (matching Masselot 2023)
varfun <- "bs"
vardegree <- 2
varper <- c(10, 75, 90)

# Number of simulations for CIs
nsim <- 1000

# Number of cores for parallel computation
ncores <- max(1, parallel::detectCores() - 2)