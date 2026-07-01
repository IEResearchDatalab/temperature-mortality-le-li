################################################################################
# 
# Reproduction and extension of Masselot et al. 2025
# 
# Part 3: Main Attribution Loop
#
################################################################################

# Load parameters and prepared data
if(file.exists("01_initialize.R")) source("01_initialize.R") else source("R_pipeline/01_initialize.R")
load(ifelse(file.exists("data/prep_data.RData"), "data/prep_data.RData", "../data/prep_data.RData"))

# Load knot information
tmean_dist <- fread("references/2025-masselot-zenodo/tmean_distribution.csv")

# Set up parallel cluster
cl <- makeCluster(ncores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = ncities, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Results storage
message("\nStarting attribution loop for ", ncities, " cities...")

# Loop over cities
foreach(i = 1:ncities, .packages = c("data.table", "dplyr", "arrow", "dlnm", "splines"), 
        .options.snow = opts) %dopar% {
  
  city_id <- cities[i]
  
  # 1. Load city-specific data
  city_meta_sub <- city_data[URAU_CODE == city_id]
  if(nrow(city_meta_sub) == 0) return(NULL)
  
  # Get knots for this city from tmean_dist
  # Percentiles 10%, 75%, 90% are at columns 21, 86, 101
  city_dist <- tmean_dist[URAU_CODE == city_id]
  knots <- as.numeric(city_dist[, c(21, 86, 101), with=F])
  bound_knots <- as.numeric(city_dist[, c(1, 118), with=F]) # 0% and 100%
  
  # 2. Get observations for bias correction
  obs_city <- obs_data[URAU_CODE == city_id, .(date, tmean = era5landtmean)]
  
  # 3. Load temperature projections
  # We use the parquet file and filter for this city
  proj_city <- open_dataset(path_tmean) %>%
    filter(URAU_CODE == city_id) %>%
    collect() %>%
    as.data.table()
  
  # Load coefficients simulations for this city
  # We read it from the large CSV using a filter (this might be slow, but safe)
  # Better: read the whole file once outside the loop? No, too large.
  # Use 'grep' to extract city-specific rows quickly
  cmd <- paste0("grep '", city_id, "' ", path_coef_simu)
  city_coef_simu <- fread(cmd = cmd)
  # Set column names back if they were lost by grep
  setnames(city_coef_simu, c("URAU_CODE", "agegroup", "sim", "b1", "b2", "b3", "b4", "b5"))
  
  # Results for this city
  city_res_list <- list()
  
  # Loop over Scenarios
  for(sc in scenarios){
    
    # Loop over GCMs
    for(gcm in gcms){
      
      # Prepare temperature series
      t_col <- paste0("tas_", gcm)
      if(!(t_col %in% names(proj_city))) next
      
      t_raw <- proj_city[ssp == sc, .(date, tmean = get(t_col))]
      
      # Bias correction
      # Train on historical period (overlapping obs and GCM historical)
      t_hist_gcm <- t_raw[year(date) %in% hist_years]
      t_obs_sub <- obs_city[year(date) %in% hist_years]
      
      # Merge to ensure alignment
      hist_merge <- merge(t_hist_gcm, t_obs_sub, by = "date", suffixes = c("_gcm", "_obs"))
      
      # Apply ISIMIP3 bias correction
      t_bc <- isimip3(t_raw$tmean, hist_merge$tmean_gcm, hist_merge$tmean_obs)
      
      # Calculate thresholds for disaggregation from BC historical distribution
      t_bc_hist <- t_bc[year(t_raw$date) %in% hist_years]
      p25 <- quantile(t_bc_hist, range_thresholds[1]/100)
      p975 <- quantile(t_bc_hist, range_thresholds[2]/100)
      
      # Loop over Age Groups
      for(age in unique(city_meta_sub$agegroup)){
        
        meta_age <- city_meta_sub[agegroup == age]
        mmt_age <- meta_age$mmt
        death_base_age <- meta_age$death_baseline / 365.25 # Daily baseline
        
        # Get 500 simulations for this age group
        sims <- city_coef_simu[agegroup == age]
        
        # Basis matrix for the temperatures
        # Use bs() with knots from tmean_dist
        B <- as.matrix(bs(t_bc, knots = knots, Boundary.knots = bound_knots, degree = 2))
        
        # Center at MMT
        B_mmt <- as.matrix(bs(mmt_age, knots = knots, Boundary.knots = bound_knots, degree = 2))
        B_centered <- sweep(B, 2, B_mmt)
        
        # Ranges masks
        mask_exc <- t_bc < p25
        mask_moc <- t_bc >= p25 & t_bc < mmt_age
        mask_moh <- t_bc >= mmt_age & t_bc <= p975
        mask_exh <- t_bc > p975
        
        # loop over 500 simulations (matrix multiplication for speed)
        coef_mat <- as.matrix(sims[, .(b1, b2, b3, b4, b5)]) # 500 x 5
        
        # RR for all days and all sims: B_centered (Days x 5) %*% t(coef_mat) (5 x 500)
        # logRR matrix: Days x 500
        logRR <- B_centered %*% t(coef_mat)
        AF <- 1 - exp(-logRR)
        
        # AN = AF * deaths
        # We assume constant daily deaths for now
        AN_all <- AF * death_base_age
        
        # Aggregate by range, decade, and sim
        years <- year(t_raw$date)
        decades <- floor(years / 10) * 10
        
        calc_range_an <- function(mask, label) {
          an_range <- AN_all
          an_range[!mask, ] <- 0
          # Sum by decade and sim
          an_dt <- as.data.table(cbind(decade = decades, an_range))
          an_dec <- an_dt[, lapply(.SD, sum), by = decade]
          
          dt_long <- melt(an_dec, id.vars = "decade", variable.name = "sim", value.name = "an")
          dt_long[, range := label]
          return(dt_long)
        }
        
        res_exc <- calc_range_an(mask_exc, "Extreme Cold")
        res_moc <- calc_range_an(mask_moc, "Moderate Cold")
        res_moh <- calc_range_an(mask_moh, "Moderate Heat")
        res_exh <- calc_range_an(mask_exh, "Extreme Heat")
        
        res_city_age_sc_gcm <- rbind(res_exc, res_moc, res_moh, res_exh)
        res_city_age_sc_gcm[, `:=`(URAU_CODE = city_id, agegroup = age, sc = sc, gcm = gcm)]
        
        city_res_list[[length(city_res_list) + 1]] <- res_city_age_sc_gcm
      }
    }
  }
  
  # Combine results for this city and save to temp
  city_res_dt <- rbindlist(city_res_list)
  saveRDS(city_res_dt, file = paste0(tdir, "/", city_id, ".rds"))
  
  return(NULL)
}

stopCluster(cl)
message("\nAttribution loop finished. Results saved in ", tdir)
