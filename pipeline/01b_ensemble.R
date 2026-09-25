#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# Pipeline Part 01b: Ensemble-mean attributable numbers across GCMs
#   Averages the Part 01 grouped ANs of the 19 GCMs (gcmlist) for one city and
#   SSP, as the GCM ensemble mean of point estimates in Masselot & Gasparrini
#   (2025) (functions/impact.R), and writes them with gcm = "ENSEMBLE" so that
#   Parts 02-04 can run on the ensemble (Simon's methods draft 2.4.5).
#   Expected layout (pipeline/run_batch.sh): <city dir>/<GCM>/ for each GCM
#   and <city dir>/ENSEMBLE/ = OUT_DIR, with GCM = "ENSEMBLE".
#
################################################################################

source("pipeline/00_pkg_params.R")

if (gcm_name != "ENSEMBLE") stop("Run Part 01b with GCM=ENSEMBLE.", call. = FALSE)
message(sprintf("\n[01b] Ensemble-mean ANs for %s %s...", city_name, ssplabs[ssp_name]))

#----- Load the per-GCM ANs

gcm_files <- file.path(dirname(out_dir), gcmlist, "01_attribution_grouped.csv")
missing <- gcmlist[!file.exists(gcm_files)]
if (length(missing)) stop(sprintf("Missing Part 01 outputs for: %s", paste(missing, collapse = ", ")), call. = FALSE)
an_all <- rbindlist(lapply(gcm_files, fread))

keys <- c("branch", "year", "agegroup", "range")
domain <- an_all[, .N, by = keys]
if (any(domain$N != length(gcmlist)) || uniqueN(an_all$gcm) != length(gcmlist)) {
  stop("Per-GCM AN domains differ; cannot average.", call. = FALSE)
}

#----- Ensemble mean of point estimates

ens <- an_all[, .(an = mean(an)), by = c(keys, "geo_id", "label", "ssp", "days_in_year", "annualization_rule")]
ens[, gcm := "ENSEMBLE"]
setcolorder(ens, names(an_all))
setorderv(ens, c("branch", "year", "agegroup", "range"))

# Spread across GCMs, kept for reporting uncertainty (Simon, 10 Sep: SE across GCMs)
spread <- an_all[, .(an_mean = mean(an), an_sd = sd(an), an_min = min(an), an_max = max(an)), by = keys]

fwrite(ens, file.path(out_dir, "01_attribution_grouped.csv"))
fwrite(spread, file.path(out_dir, "01b_gcm_spread.csv"))
message("Saved ensemble-mean ANs (", length(gcmlist), " GCMs) to ", file.path(out_dir, "01_attribution_grouped.csv"))
