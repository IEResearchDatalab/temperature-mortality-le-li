################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 02_ComputeAN.R
#   Compute attributable numbers for all cities
#   Sequential over cities (parquet I/O bound), parallel within city
#
################################################################################

if (length(ls()) == 0) source("01_PrepData.R")

#---------------------------
# Log file
#---------------------------

dir.create("temp", showWarnings = FALSE)
writeLines(c(""), "temp/log_an.txt")
cat(as.character(as.POSIXct(Sys.time())), "\n", file = "temp/log_an.txt", append = TRUE)

cat("Computing AN for", n_cities, "cities...\n")

#---------------------------
# Sequential loop over cities (parquet read is I/O bound)
#---------------------------

all_results <- list()

for (ic in seq_len(n_cities)) {

  cc <- city_codes[ic]
  cat("\n", "city ", ic, "/", n_cities, ": ", cc, " ", 
    as.character(Sys.time()), "\n", sep = "",
    file = "temp/log_an.txt", append = TRUE)

  #----- Load temperature data (once per city)
  temp_all <- as.data.table(read_parquet(
    file.path(path_data, "tmeanproj.gz.parquet")))
  temp_city <- temp_all[URAU_CODE == cc]
  if (nrow(temp_city) == 0) next
  rm(temp_all)

  temp_city[, year := year(date)]
  if (!"doy" %in% names(temp_city))
    temp_city[, doy := as.integer(format(date, "%j"))]

  id_cols <- c("URAU_CODE", "date", "year", "doy", "ssp")
  gcm_cols <- setdiff(names(temp_city), id_cols)
  if (length(gcm_cols) == 0) next

  # Historical temperatures for basis
  hist_data <- temp_city[ssp == "hist" & year >= hist_year_min &
    year <= hist_year_max]
  hist_temps <- unlist(hist_data[, ..gcm_cols], use.names = FALSE)
  hist_temps <- hist_temps[!is.na(hist_temps)]
  if (length(hist_temps) < 100) next

  varknots <- quantile(hist_temps, varper / 100, na.rm = TRUE)
  varbound <- range(hist_temps, na.rm = TRUE)
  argvar <- list(fun = varfun, degree = vardegree,
    knots = varknots, Bound = varbound)
  p025_hist <- as.numeric(quantile(hist_temps, cold_extreme_pct, na.rm = TRUE))
  p975_hist <- as.numeric(quantile(hist_temps, heat_extreme_pct, na.rm = TRUE))

  coefs_city <- coefs_all[URAU_CODE == cc]
  coef_simu_city <- coef_simu[URAU_CODE == cc]
  city_base <- city_res[URAU_CODE == cc & agegroup %in% agelabs]
  base_death <- setNames(city_base$death, city_base$agegroup)

  city_results <- list()

  for (ssp_v in intersect(ssp_keep[-1], unique(temp_city$ssp))) {

    proj_data <- temp_city[ssp == ssp_v & year >= proj_year_min &
      year <= proj_year_max]

    for (gcm_col in gcm_cols) {

      proj <- proj_data[, .(date, year, doy, temp = get(gcm_col))]
      proj <- proj[!is.na(temp)]
      if (nrow(proj) == 0) next
      proj[, n_days_year := .N, by = year]

      for (ag in agelabs) {
        coef_row <- coefs_city[agegroup == ag]
        if (nrow(coef_row) == 0) next
        coef_vec <- as.numeric(coef_row[, .(b1, b2, b3, b4, b5)])

        # MMT
        temp_seq <- seq(varbound[1], varbound[2], by = 0.5)
        basis_seq <- do.call(onebasis, c(list(x = temp_seq), argvar))
        log_rr_seq <- as.vector(basis_seq %*% coef_vec)
        mmt_idx <- which.min(log_rr_seq[
          temp_seq >= quantile(temp_seq, 0.25) &
          temp_seq <= quantile(temp_seq, 0.99)])
        mmt <- temp_seq[mmt_idx]

        # Daily AN (Masselot-style)
        bvar <- do.call(onebasis, c(list(x = proj$temp), argvar))
        cenvec <- do.call(onebasis, c(list(x = mmt), argvar))
        bvarcen <- scale(bvar, center = cenvec, scale = FALSE)

        af_day <- 1 - exp(-bvarcen %*% coef_vec)
        annual_deaths <- base_death[ag]
        daily_deaths <- annual_deaths / 365
        an_day <- as.vector(af_day * daily_deaths)

        # Temperature range classification
        temp_range <- fifelse(proj$temp < p025_hist, "extreme_cold",
          fifelse(proj$temp < mmt, "moderate_cold",
            fifelse(proj$temp <= p975_hist, "moderate_heat", "extreme_heat")))

        # Annual aggregation
        an_annual <- tapply(an_day, list(proj$year, temp_range), sum, default = 0)
        an_annual_dt <- as.data.table(an_annual)
        an_annual_dt[, year := as.integer(rownames(an_annual))]

        # CIs from pre-simulated coefficients
        sim_coefs <- coef_simu_city[agegroup == ag]
        n_sims <- uniqueN(sim_coefs$sim)

        an_sim <- t(sapply(seq_len(n_sims), function(s) {
          cf <- as.numeric(sim_coefs[sim == s, .(b1, b2, b3, b4, b5)])
          af_s <- 1 - exp(-bvarcen %*% cf)
          an_s <- as.vector(af_s * daily_deaths)
          c(total = sum(an_s),
            extreme_cold = sum(an_s[temp_range == "extreme_cold"]),
            moderate_cold = sum(an_s[temp_range == "moderate_cold"]),
            moderate_heat = sum(an_s[temp_range == "moderate_heat"]),
            extreme_heat = sum(an_s[temp_range == "extreme_heat"]))
        }))

        ci <- apply(an_sim, 2, quantile, c(.025, .975), na.rm = TRUE)

        for (tr in c("extreme_cold", "moderate_cold", "moderate_heat",
          "extreme_heat")) {
          an_val <- if (tr %in% colnames(an_annual)) an_annual_dt[[tr]] else 0
          an_low_val <- if (tr %in% colnames(ci)) ci[1, tr] else 0
          an_hi_val <- if (tr %in% colnames(ci)) ci[2, tr] else 0
          city_results[[length(city_results) + 1]] <- data.table(
            city_code = cc, gcm = gcm_col, ssp = ssp_v,
            age_group = ag, temp_range = tr,
            year = an_annual_dt$year,
            an_est = an_val,
            an_low = an_low_val,
            an_hi = an_hi_val)
        }
      }
    }
  }

  if (length(city_results) > 0) {
    all_results[[length(all_results) + 1]] <- rbindlist(city_results, fill = TRUE)
    cat(sprintf("  City %s: %d rows\n", cc, nrow(all_results[[length(all_results)]])))
  }
}

all_results <- rbindlist(all_results, fill = TRUE)
cat("AN complete:", format(nrow(all_results), big.mark = ","), "rows\n")
fwrite(all_results, file.path(path_out, "ans_annual_all_cities.csv"))
cat("Written:", file.path(path_out, "ans_annual_all_cities.csv"), "\n")