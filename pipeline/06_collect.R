#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# Pipeline Part 06: Collect batch results into the three data objects
#   Reads results/europe/ssp<k>/<city>/ written by run_batch.sh and stores,
#   at the most disaggregated level (meeting 10 Sep 2026, §6), as parquet:
#     object1_dataset.parquet  city x ssp x scenario x year x age x cause:
#                              deaths and population (ensemble-mean ANs;
#                              Lloyd et al. 2024 Fig S1 "dataset for analysis")
#     object2_contributions.parquet  Horiuchi contributions by city x ssp x
#                              age x cause: with - without CC per 5-year period
#                              ("between") and change between consecutive
#                              5-year periods within each scenario ("within")
#     object3_levels.parquet   LE65 and LI65+ by city x ssp x scenario x year,
#                              for the ensemble and each of the 19 GCMs
#   Country and region are attached from city_results.csv (Masselot regions);
#   higher levels are derived by aggregating cities.
#
################################################################################

source("pipeline/00_pkg_params.R")

root <- Sys.getenv("BATCH_ROOT", "results/europe")
out <- file.path(root, "collected")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta <- unique(fread("data/city_results.csv")[, .(city = URAU_CODE, city_name = LABEL, country = CNTR_CODE, region)])

#----- Completed city x SSP runs

ens_dirs <- Sys.glob(file.path(root, "ssp*", "*", "ENSEMBLE"))
ens_dirs <- ens_dirs[file.exists(file.path(ens_dirs, ".done"))]
if (!length(ens_dirs)) stop("No completed ENSEMBLE runs under ", root, call. = FALSE)
runs <- data.table(dir = ens_dirs, city = basename(dirname(ens_dirs)),
  ssp = as.integer(sub("ssp", "", basename(dirname(dirname(ens_dirs))))))
message(sprintf("\n[06] Collecting %d city x SSP runs (%d cities)...", nrow(runs), uniqueN(runs$city)))

#----- Object 1: dataset for analysis

obj1 <- rbindlist(lapply(seq_len(nrow(runs)), function(i) {
  m <- fread(file.path(runs$dir[i], "03_master_table.csv"))
  m[, .(city = geo_id, ssp = runs$ssp[i], scenario = branch, year, age, cause = range, deaths = deaths_component, pop)]
}))
obj1 <- merge(meta, obj1, by = "city")
setorder(obj1, city, ssp, scenario, year, age, cause)
write_parquet(obj1, file.path(out, "object1_dataset.parquet"))

#----- Object 2: decomposition contributions

obj2 <- rbindlist(lapply(seq_len(nrow(runs)), function(i) {
  b <- fread(file.path(runs$dir[i], "04_between_branch_decomposition.csv"))
  w <- fread(file.path(runs$dir[i], "04_within_branch_period_decomposition.csv"))
  rbind(
    b[, .(city = geo_id, ssp = runs$ssp[i], type = "between", scenario = "with_cc - without_cc",
      period_from = period, period_to = period, age, cause, le_contribution, li_contribution)],
    w[, .(city = geo_id, ssp = runs$ssp[i], type = "within", scenario = branch,
      period_from, period_to, age, cause, le_contribution, li_contribution)]
  )
}))
obj2 <- merge(meta, obj2, by = "city")
write_parquet(obj2, file.path(out, "object2_contributions.parquet"))

#----- Object 3: LE65 / LI65+ levels (ensemble and each GCM)

lev_files <- unlist(lapply(dirname(runs$dir), function(d) file.path(d, c("ENSEMBLE", gcmlist), "04_le_li_levels.csv")))
lev_files <- lev_files[file.exists(lev_files)]
obj3 <- rbindlist(lapply(lev_files, function(f) {
  l <- fread(f)
  l[, .(city = geo_id, ssp = as.integer(sub("ssp", "", basename(dirname(dirname(dirname(f)))))),
    gcm, scenario = branch, year, LE65, LI65)]
}))
obj3 <- merge(meta, obj3, by = "city")
write_parquet(obj3, file.path(out, "object3_levels.parquet"))

#----- Completeness report

report <- obj3[, .(n_gcm = uniqueN(gcm[gcm != "ENSEMBLE"]), has_ensemble = any(gcm == "ENSEMBLE")), by = .(city, ssp)]
fwrite(report, file.path(out, "completeness.csv"))
message(sprintf("Object 1: %d rows | Object 2: %d rows | Object 3: %d rows | runs with 19 GCMs: %d/%d",
  nrow(obj1), nrow(obj2), nrow(obj3), sum(report$n_gcm == length(gcmlist)), nrow(report)))
message("Saved to ", out)
