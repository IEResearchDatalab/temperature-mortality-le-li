################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#   a multi-city, multi-GCM, multi-SSP health impact assessment
#
# Master pipeline script
#   Sources all scripts in cascading chain (Masselot-style)
#
# Usage: Rscript 00_RunAll.R
#
################################################################################

#----------------------
# Run the full pipeline
#----------------------

source("01_PrepData.R")
source("02_ComputeAN.R")
source("03_Disaggregate.R")
source("04_AnalysisDataset.R")
source("05_LifeTables.R")
source("06_Decomposition.R")

cat("\n=== Pipeline complete ===\n")
cat(sprintf("Output directory: %s\n",
  normalizePath("../results/production")))