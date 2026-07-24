################################################################################
# Baseline attributable numbers using ERA5 observations (Masselot 2023 method)
# 4-range disaggregation: ExtrCold, ModCold, ModHeat, ExtrHeat
################################################################################

library(data.table)
library(dlnm)
library(splines)
library(doSNOW)
library(foreach)
library(dplyr)

source("R_pipeline/01_initialize.R")
load("data/prep_data.RData")

dir.create("temp_results_baseline", showWarnings = FALSE)

city_filter <- trimws(Sys.getenv("CITY_FILTER", unset = ""))
if (nzchar(city_filter)) {
  cities <- intersect(cities, trimws(strsplit(city_filter, ",", fixed = TRUE)[[1]]))
}

message("\n=== ERA5 Baseline Attribution (2023 formula, 4 ranges) ===\n")

cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterExport(cl, c("nsim", "knots_percentiles", "thresholds", "obs_data", 
                   "varfun", "vardegree", "path_coef_simu", "path_coefs"))

pb <- txtProgressBar(max = length(cities), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results <- foreach(city_id = cities, .packages = c("data.table", "dlnm", "splines", "dplyr"), .options.snow = opts) %dopar% {
  out_path <- paste0("temp_results_baseline/", city_id, ".rds")
  if (file.exists(out_path)) return(paste0("Skipped: ", city_id))

  city_thresholds <- thresholds[URAU_CODE == city_id]
  if (nrow(city_thresholds) == 0) return(paste0("No thresholds: ", city_id))

  dt <- obs_data[URAU_CODE == city_id]
  dt[, yr := as.integer(format(date, "%Y"))]
  if (nrow(dt) == 0) return(paste0("No ERA5: ", city_id))

  p2_5 <- city_thresholds$p2_5[1]
  p97_5 <- city_thresholds$p97_5[1]
  temp_vals <- dt$tmean_obs
  knots <- quantile(temp_vals, knots_percentiles / 100, na.rm = TRUE)
  bound <- range(temp_vals, na.rm = TRUE)

  coef_pt <- fread(cmd = paste0("grep '", city_id, "' ", path_coefs))
  coef_sims <- fread(cmd = paste0("grep '", city_id, "' ", path_coef_simu))
  if (nrow(coef_pt) == 0 || nrow(coef_sims) == 0) return(paste0("No coefficients: ", city_id))
  setnames(coef_pt, c("URAU_CODE", "agegroup", paste0("b", 1:5)))
  setnames(coef_sims, c("URAU_CODE", "agegroup", "sim", paste0("b", 1:5)))
  coef_pt[, sim := 0L]
  coef_sims <- coef_sims[sim <= nsim]
  city_coefs <- rbindlist(list(
    coef_pt[, .(URAU_CODE, agegroup, sim, b1, b2, b3, b4, b5)],
    coef_sims[, .(URAU_CODE, agegroup, sim, b1, b2, b3, b4, b5)]
  ), use.names = TRUE)

  city_results_list <- list()

  for (agegrp in unique(city_thresholds$agegroup)) {
    age_thresh <- city_thresholds[agegroup == agegrp]
    mmt <- age_thresh$mmt
    death_annual <- age_thresh$death

    b_temp <- onebasis(dt$tmean_obs, fun = "bs", degree = 2, knots = knots)
    b_mmt <- onebasis(mmt, fun = "bs", degree = 2, knots = knots,
              Boundary.knots = bound)
    b_centered <- scale(b_temp, center = b_mmt, scale = FALSE)

    age_coefs <- as.matrix(city_coefs[agegroup == agegrp][order(sim), .(b1, b2, b3, b4, b5)])

    log_rr <- b_centered %*% t(age_coefs)
    af <- 1 - exp(-log_rr)
    an <- af * death_annual

    range_idx <- case_when(
      dt$tmean_obs < p2_5 ~ "ExtrCold",
      dt$tmean_obs < mmt ~ "ModCold",
      dt$tmean_obs < p97_5 ~ "ModHeat",
      TRUE ~ "ExtrHeat"
    )

    groups <- paste(dt$yr, range_idx, sep = "::")
    an_agg <- rowsum(an, groups)

    yr_days <- as.numeric(table(dt$yr))
    yr_labels <- names(table(dt$yr))
    grp_yr <- sub("::.*", "", rownames(an_agg))
    an_agg <- an_agg / yr_days[match(grp_yr, yr_labels)]

    an_dt <- as.data.table(an_agg, keep.rownames = "group")
    an_long <- melt(an_dt, id.vars = "group", variable.name = "sim", value.name = "an")
    an_long[, sim := as.integer(sub("^V", "", sim)) - 1L]
    an_long[, c("year", "range") := tstrsplit(group, "::")]
    an_long[, `:=`(
      year = as.integer(year),
      group = NULL,
      ssp = 0,
      gcm = "ERA5",
      agegroup = agegrp
    )]

    city_results_list[[length(city_results_list) + 1]] <- an_long
  }

  if (length(city_results_list) > 0) {
    city_total <- rbindlist(city_results_list)
    saveRDS(city_total, out_path)
    return(paste0("Success: ", city_id))
  }
  return(paste0("Empty: ", city_id))
}

stopCluster(cl)
close(pb)

errors <- results[!grepl("^(Success|Skipped)", results)]
message("\nBaseline attribution complete. ", length(cities) - length(errors), "/", length(cities), " cities processed.")
if (length(errors) > 0) message("Errors: ", length(errors))