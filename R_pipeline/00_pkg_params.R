################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# R Pipeline Part 0: Packages and analysis parameters
#   Sourced at the top of every other script. Structure and parameter names
#   follow Masselot & Gasparrini (2025) 01_pkg_params.R wherever the same
#   quantity is used.
#
################################################################################

#------------------------
# LOAD THE PACKAGES
#------------------------

suppressPackageStartupMessages({

  #----- Data management
  library(data.table) # For very large databases
  library(dplyr) # Data.frame management (arrow queries)
  library(arrow) # To deal with datasets that cannot be loaded all at once

  #----- Statistical analysis
  library(dlnm); library(splines) # Create bases for RR computation
  library(ungroup) # PCLM disaggregation to single ages (Rizzi et al. 2015)
  library(DemoDecomp) # Horiuchi decomposition (Aburto et al. 2022)

  #----- Plotting
  library(ggplot2) # For plots
  library(patchwork) # Put ggplots together
})

#------------------------
# PARAMETERS
#------------------------

#----- Study unit
# Defaults run the single-city validation (Madrid, SSP3-7.0, GFDL-ESM4). Batch
# runs (R_pipeline/run_batch.sh) override them through environment variables.

city_id <- Sys.getenv("CITY_ID", "ES001C") # URAU code
ssp_name <- Sys.getenv("SSP", "3") # as in the `ssp` column of tmeanproj / Wittgenstein
gcm_name <- Sys.getenv("GCM", "GFDL_ESM4") # or "ENSEMBLE" (mean ANs over the 19 GCMs, Part 01b)
city_name <- local({
  lab <- unique(fread("data/city_results.csv", select = c("URAU_CODE", "LABEL"))[URAU_CODE == city_id]$LABEL)
  if (length(lab) != 1L) stop(sprintf("Unknown city %s.", city_id), call. = FALSE)
  lab
})

# GCMs to exclude (Masselot 2025): 19 of the 21 GCMs in tmeanproj remain
gcmexcl <- c("CMCC_CM2_SR5", "TaiESM1")
gcmlist <- c("ACCESS_CM2", "ACCESS_ESM1_5", "BCC_CSM2_MR", "CMCC_ESM2", "CanESM5",
  "EC_Earth3", "EC_Earth3_Veg_LR", "GFDL_ESM4", "IITM_ESM", "INM_CM4_8", "INM_CM5_0",
  "IPSL_CM6A_LR", "KACE_1_0_G", "MIROC6", "MPI_ESM1_2_HR", "MPI_ESM1_2_LR",
  "MRI_ESM2_0", "NorESM2_LM", "NorESM2_MM")
if (gcm_name %in% gcmexcl) stop(sprintf("GCM %s is excluded in Masselot (2025).", gcm_name), call. = FALSE)

#----- Ages

# Age groups of the ERFs used here (Masselot's agelabs restricted to 65+)
agelabs <- c("65-74", "75-84", "85+")
age_slices <- list("65-74" = 65:74, "75-84" = 75:84, "85+" = 85:100)

# Single ages of the life tables; 100 is the open interval 100+
age_levels <- 65:100
nx <- c(rep(1, 100 - 65), Inf)

#----- Periods

# Historical (calibration) period and projection calibration periods (Masselot)
histrange <- c(2000, 2014)
projrange <- c(2015, seq(2030, 2100, by = 10))

# Length of the periods of the demographic projections (in years)
perlen <- 5

# Years reported
future_years <- 2020:2099

#----- Exposure-response function (follows Masselot et al 2023 Lancet Plan. Health)

varfun <- "bs"
vardegree <- 2
varper <- c(10, 75, 90)

# Temperature percentiles (used for the ERF basis and the MMT)
predper <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))

#----- Temperature ranges (Lloyd et al. 2024)

range_levels <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
cause_levels <- c(range_levels, "rest")

# Percentiles delimiting extreme cold / heat, and the ERA5 period they come from
# (the ERF estimation window). Open point: Simon's methods draft says 2000-2014.
extreme_probs <- c(p2_5 = 0.025, p97_5 = 0.975)
threshold_years <- c(1990L, 2019L)

#----- Scenarios

branch_levels <- c("with_cc", "without_cc")

# Without-climate-change series: "masselot_demo" (Masselot 2025: each 5-year
# block recalibrated to the 2010-2014 distribution) or "era5_cycle" (observed
# 2000-2019 repeated forward, sensitivity analysis)
counterfactual <- "masselot_demo"
counterfactual_ref_year5 <- 2010L
hist_years_counterfactual <- 2000:2019

#----- Single-age disaggregation (PCLM)

pclm_input_scale <- 1000 # Wittgenstein counts are in thousands; fit on persons
pclm_open_nlast <- 11L # spread 100+ over 100-110, then collapse back to 100+

#----- Decomposition

N_HORIUCHI <- as.integer(Sys.getenv("N_HORIUCHI", "400"))
closure_tol <- 1e-6

# Part 04 options: year-on-year decomposition (158 Horiuchi steps, the costly
# part) and levels-only runs (per-GCM uncertainty in batch runs)
decomp_annual <- Sys.getenv("DECOMP_ANNUAL", "1") == "1"
levels_only <- Sys.getenv("LEVELS_ONLY", "0") == "1"

#------------------------
# OUTPUT DIRECTORIES
#------------------------

out_dir <- Sys.getenv("OUT_DIR", "results/phase1_madrid")
check_dir <- Sys.getenv("CHECK_DIR", "results/checks")
fig_dir <- Sys.getenv("FIG_DIR", "results/figures")
# Part 00 outputs (shared by all GCMs of a city in batch runs)
dem_dir <- Sys.getenv("DEMOG_DIR", out_dir)
for (d in c(out_dir, check_dir, fig_dir, dem_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

#------------------------
# AESTHETICS
#------------------------

range_labels <- c(ExtrCold = "Extreme cold", ModCold = "Moderate cold",
  ModHeat = "Moderate heat", ExtrHeat = "Extreme heat")
range_colors <- c(ExtrCold = "#2166ac", ModCold = "#67a9cf",
  ModHeat = "#ef8a62", ExtrHeat = "#b2182b")
branch_labels <- c(with_cc = "With climate change", without_cc = "Without climate change")
branch_colors <- c(with_cc = "#b2182b", without_cc = "#2166ac")
ssplabs <- c("1" = "SSP1-2.6", "2" = "SSP2-4.5", "3" = "SSP3-7.0")
scenario_label <- sprintf("%s, %s, central ERF estimates", ssplabs[ssp_name], gsub("_", "-", gcm_name))
