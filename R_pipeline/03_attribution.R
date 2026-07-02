################################################################################
# 
# Contrasting future heat and cold-related mortality under climate change, 
# demographic and adaptation scenarios in 854 European cities
#
# R Code Part 3: Perform health impact projections across cities
#
# Adapted from Pierre Masselot & Antonio Gasparrini
#
################################################################################

library(data.table)
library(arrow)
library(doSNOW)
library(dlnm)
library(splines)
library(dplyr)

source("R_pipeline/01_initialize.R")

# Inlining isimip3 function to avoid external dependencies
isimip3 <- function(obshist, simhist, simfut, 
  yearobshist, yearsimhist, yearsimfut, detrend = T)
{
  if (detrend){
    obstrend <- lm(obshist ~ yearobshist, na.action = na.exclude) |> 
      predict() |> scale(scale = F)
    simhisttrend <- lm(simhist ~ yearsimhist, na.action = na.exclude) |> 
        predict() |> scale(scale = F)
    simfuttrend <- lm(simfut ~ yearsimfut, na.action = na.exclude) |> 
      predict() |> scale(scale = F)
    
    obshist <- obshist - obstrend
    simhist <- simhist - simhisttrend
    simfut <- simfut - simfuttrend
  }
  
  ecdfobs <- ecdf(obshist)(obshist)
  deltaadd <- quantile(simfut, ecdfobs, na.rm = T) - quantile(simhist, ecdfobs, na.rm = T)
  obsfut <- deltaadd + obshist
  simfutcdf <- pnorm(simfut, mean(simfut, na.rm = T), sd(simfut, na.rm = T))
  calsimfut <- qnorm(p = simfutcdf, mean = mean(obsfut, na.rm = T), sd = sd(obsfut, na.rm = T))
  
  if (detrend){
    calsimfut <- calsimfut + simfuttrend
  }
  calsimfut
}

message("\n[2/3] Starting attribution simulations...")

load("data/prep_data.RData")
dir.create("temp_results", showWarnings = FALSE)

cl <- makeCluster(n_cores)
registerDoSNOW(cl)

clusterExport(cl, c("isimip3", "nsim", "scenarios", "gcms", "hist_years", 
                   "knots_percentiles", "path_tmean", "path_coef_simu",
                   "thresholds", "obs_data"))

pb <- txtProgressBar(max = length(cities), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results <- foreach(city_id = cities, .packages = c("data.table", "arrow", "dlnm", "splines", "dplyr"), .options.snow = opts) %dopar% {
  
  # Check if result already exists to allow resuming
  out_path <- paste0("temp_results/", city_id, ".rds")
  if(file.exists(out_path)) return(paste0("Skipped: ", city_id))

  city_thresholds <- thresholds[URAU_CODE == city_id]
  if(nrow(city_thresholds) == 0) return(paste0("No thresholds found for ", city_id))
  
  tmean_obs_city <- obs_data[URAU_CODE == city_id, .(date, tmean_obs, year = year(date), month = month(date))]
  p2_5 <- city_thresholds$p2_5[1]
  p97_5 <- city_thresholds$p97_5[1]
  
  proj_ds <- open_dataset(path_tmean)
  tmean_proj_all <- proj_ds %>% filter(URAU_CODE == city_id) %>% collect() %>% as.data.table()
  
  city_coefs_raw <- fread(cmd = paste0("grep '", city_id, "' ", path_coef_simu))
  if(nrow(city_coefs_raw) == 0) return(paste0("No coefficients found for ", city_id))
  setnames(city_coefs_raw, c("URAU_CODE", "agegroup", "sim", paste0("b", 1:5)))
  city_coefs_raw <- city_coefs_raw[sim <= nsim]
  
  city_results_list <- list()
  
  for (sc in scenarios) {
    for (gcm in gcms) {
      col_name <- paste0("tas_", gcm)
      if(!(col_name %in% colnames(tmean_proj_all))) next
      
      t_proj_sc_gcm <- tmean_proj_all[ssp %in% c("hist", sc), c("date", "ssp", col_name), with = FALSE]
      setnames(t_proj_sc_gcm, col_name, "tmean")
      t_proj_sc_gcm[, `:=`(year = year(date), month = month(date))]
      
      t_proj_sc_gcm[, tmean_bc := {
        m <- .BY$month
        obs_m <- tmean_obs_city[month == m & year %in% hist_years]
        sim_h_m <- t_proj_sc_gcm[ssp == "hist" & month == m & year %in% hist_years]
        
        if(nrow(obs_m) < 10 || nrow(sim_h_m) < 10) return(as.numeric(NA))
        
        isimip3(
          obshist = obs_m$tmean_obs,
          simhist = sim_h_m$tmean,
          simfut = tmean,
          yearobshist = obs_m$year,
          yearsimhist = sim_h_m$year,
          yearsimfut = year
        )
      }, by = month]
      
      for (agegrp in unique(city_thresholds$agegroup)) {
        age_thresh <- city_thresholds[agegroup == agegrp]
        mmt <- age_thresh$mmt
        daily_deaths <- age_thresh$death / 365.25
        
        obs_temp_vals <- tmean_obs_city$tmean_obs
        knots <- quantile(obs_temp_vals, knots_percentiles/100, na.rm=TRUE)
        bound <- range(obs_temp_vals, na.rm=TRUE)
        
        b_fut <- onebasis(t_proj_sc_gcm$tmean_bc, fun="ns", knots=knots, Bound=bound, intercept = TRUE)
        b_mmt <- onebasis(mmt, fun="ns", knots=knots, Bound=bound, intercept = TRUE)
        b_fut_centered <- scale(b_fut, center = b_mmt, scale = FALSE)
        
        age_coefs <- as.matrix(city_coefs_raw[agegroup == agegrp, .(b1, b2, b3, b4, b5)])
        rr <- exp(b_fut_centered %*% t(age_coefs))
        an <- (1 - 1/rr) * daily_deaths
        
        range_idx <- case_when(
          t_proj_sc_gcm$tmean_bc < p2_5 ~ "ExtrCold",
          t_proj_sc_gcm$tmean_bc < mmt ~ "ModCold",
          t_proj_sc_gcm$tmean_bc < p97_5 ~ "ModHeat",
          TRUE ~ "ExtrHeat"
        )
        
        groups <- paste(t_proj_sc_gcm$year, range_idx, sep = "::")
        an_agg <- rowsum(an, groups)
        
        an_agg_dt <- as.data.table(an_agg, keep.rownames = "group")
        an_agg_long <- melt(an_agg_dt, id.vars = "group", variable.name = "sim", value.name = "an")
        an_agg_long[, c("year", "range") := tstrsplit(group, "::")]
        an_agg_long[, `:=`(year = as.integer(year), group = NULL, ssp = sc, gcm = gcm, agegroup = agegrp)]
        
        city_results_list[[length(city_results_list) + 1]] <- an_agg_long
      }
    }
  }
  
  if(length(city_results_list) > 0) {
    city_total <- rbindlist(city_results_list)
    saveRDS(city_total, paste0("temp_results/", city_id, ".rds"))
    return(paste0("Success: ", city_id))
  } else {
    return(paste0("Empty: ", city_id))
  }
}

stopCluster(cl)
close(pb)

errors <- results[!grepl("^(Success|Skipped)", results)]
if(length(errors) > 0) {
  message("\nNotice: ", length(errors), " cities were not processed successfully.")
} else {
  message("\nSimulations completed for all cities.")
}
