################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 02_ComputeAN.R
#   Compute climate-only attributable numbers (ANs) under fixed baseline deaths.
#
# This is the epidemiological attribution layer only.
# It uses fixed baseline deaths for all future years, producing ANs driven
# solely by changing temperatures (climate effect).
#
# Output: annual ANs by city × year × GCM × SSP × age group × temperature range
#   with annual uncertainty intervals.
#
# Critical implementation rules:
#   - basis definition matches the ERF estimation source (Masselot 2023)
#   - MMT is found using the correct empirical percentile search
#   - uncertainty is computed per year, not as one CI pasted onto all years
#   - AF handling is consistent between point estimate and simulations
#   - output is clearly labeled as climate-only attribution
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")

#---------------------------
# Load data prepared by 01_PrepData.R
#---------------------------

coefs_all <- readRDS(file.path(path_out, "coefs_all.rds"))
if (!is.null(city_subset)) {
  coefs_all <- coefs_all[URAU_CODE %in% city_subset]
}
coef_simu <- readRDS(file.path(path_out, "coef_simu.rds"))
if (!is.null(city_subset)) {
  coef_simu <- coef_simu[URAU_CODE %in% city_subset]
}
city_res <- readRDS(file.path(path_out, "city_res.rds"))
if (!is.null(city_subset)) {
  city_res <- city_res[URAU_CODE %in% city_subset]
}

city_codes <- sort(unique(coefs_all$URAU_CODE))
n_cities <- length(city_codes)

cat(sprintf("Loaded data for %d cities\n", n_cities))

#---------------------------
# Log file
#---------------------------

dir.create("temp", showWarnings = FALSE)
writeLines(c(""), "temp/log_an.txt")
cat(as.character(as.POSIXct(Sys.time())), "\n", file = "temp/log_an.txt", append = TRUE)

cat("Computing climate-only ANs for", n_cities, "cities...\n")

#----- Checkpoint directory (per-city CSVs for resumability)
an_checkpoint <- file.path(path_out, "an_checkpoint")
dir.create(an_checkpoint, showWarnings = FALSE)

#---------------------------
# Sequential loop over cities
#---------------------------

for (ic in seq_len(n_cities)) {

  cc <- city_codes[ic]
  checkpoint_file <- file.path(an_checkpoint, paste0(cc, ".csv"))

  if (file.exists(checkpoint_file)) {
    cat(sprintf("  Skipping %s (checkpoint found)\n", cc))
    next
  }

  cat("\n", "city ", ic, "/", n_cities, ": ", cc, " ",
    as.character(Sys.time()), "\n", sep = "",
    file = "temp/log_an.txt", append = TRUE)

  #----- Load temperature data via arrow (predicate pushdown — reads only this city)
  temp_city <- as.data.table(
    arrow::open_dataset(file.path(path_data, "tmeanproj.gz.parquet")) |>
      dplyr::filter(URAU_CODE == cc) |>
      dplyr::collect()
  )

  temp_city[, year := year(date)]
  if (!"doy" %in% names(temp_city))
    temp_city[, doy := as.integer(format(date, "%j"))]

  id_cols <- c("URAU_CODE", "date", "year", "doy", "ssp")
  gcm_cols <- setdiff(names(temp_city), id_cols)
  if (length(gcm_cols) == 0) next

  # Historical temperatures for basis construction
  hist_data <- temp_city[ssp == "hist" & year >= hist_year_min &
    year <= hist_year_max]
  hist_temps <- unlist(hist_data[, ..gcm_cols], use.names = FALSE)
  hist_temps <- hist_temps[!is.na(hist_temps)]
  if (length(hist_temps) < 100) next

  # Basis parameters (matching Masselot 2023)
  varknots <- quantile(hist_temps, varper / 100, na.rm = TRUE)
  varbound <- range(hist_temps, na.rm = TRUE)
  argvar <- list(fun = varfun, degree = vardegree,
    knots = varknots, Bound = varbound)

  # Historical percentiles for temperature range classification
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

        #----- MMT: search over empirical temperature distribution
        temp_seq <- seq(varbound[1], varbound[2], by = 0.5)
        basis_seq <- do.call(onebasis, c(list(x = temp_seq), argvar))
        log_rr_seq <- as.vector(basis_seq %*% coef_vec)

        # Restrict search to the central 75% of the empirical temp distribution
        # (quantiles of hist_temps, not of temp_seq)
        mmt_lower <- quantile(hist_temps, 0.25, na.rm = TRUE)
        mmt_upper <- quantile(hist_temps, 0.99, na.rm = TRUE)
        mmt_search_idx <- which(temp_seq >= mmt_lower & temp_seq <= mmt_upper)
        mmt <- temp_seq[mmt_search_idx][which.min(log_rr_seq[mmt_search_idx])]

        #----- Daily attributable fraction and AN
        bvar <- do.call(onebasis, c(list(x = proj$temp), argvar))
        cenvec <- do.call(onebasis, c(list(x = mmt), argvar))
        bvarcen <- scale(bvar, center = cenvec, scale = FALSE)

        # Point estimate: AF, truncated at zero (Gasparrini & Leone 2014)
        af_day <- 1 - exp(-bvarcen %*% coef_vec)
        af_day[af_day < 0] <- 0

        annual_deaths <- base_death[ag]
        daily_deaths <- annual_deaths / 365
        an_day <- as.vector(af_day * daily_deaths)

        # Temperature range classification
        temp_range <- fifelse(proj$temp < p025_hist, "extreme_cold",
          fifelse(proj$temp < mmt, "moderate_cold",
            fifelse(proj$temp <= p975_hist, "moderate_heat", "extreme_heat")))

        #----- Annual aggregation (point estimate)
        years_in_data <- sort(unique(proj$year))
        an_annual <- tapply(an_day, list(proj$year, temp_range), sum, default = 0)
        an_annual_dt <- as.data.table(an_annual)
        an_annual_dt[, year := as.integer(rownames(an_annual))]

        #----- Annual uncertainty: vectorised over all 1000 sims at once
        sim_coefs <- coef_simu_city[agegroup == ag][order(sim)]
        coef_sims_mat <- as.matrix(sim_coefs[, .(b1, b2, b3, b4, b5)])  # n_sims × 5

        # All AFs in one matrix multiply: n_days × n_sims
        af_all_mat <- 1 - exp(-bvarcen %*% t(coef_sims_mat))
        af_all_mat[af_all_mat < 0] <- 0
        an_all_mat <- af_all_mat * daily_deaths  # broadcast scalar

        range_names <- c("extreme_cold", "moderate_cold", "moderate_heat", "extreme_heat")
        range_idx_vec <- match(temp_range, range_names)

        # Compute annual CIs for each range using rowsum
        for (r_idx in seq_along(range_names)) {
          rname <- range_names[r_idx]
          row_sel <- which(range_idx_vec == r_idx)

          if (length(row_sel) == 0) {
            # No days in this range: zero CI for all years
            for (yr_idx in seq_along(years_in_data)) {
              yr <- years_in_data[yr_idx]
              an_val <- 0
              city_results[[length(city_results) + 1]] <- data.table(
                city_code = cc, gcm = gcm_col, ssp = ssp_v,
                age_group = ag, temp_range = rname, year = yr,
                an_est = an_val, an_low = 0, an_hi = 0)
            }
            next
          }

          # Annual sims for this range: rowsum over year × sims submatrix
          an_range_mat <- an_all_mat[row_sel, , drop = FALSE]  # n_range_days × n_sims
          yr_range_vec <- proj$year[row_sel]
          ann_sims <- rowsum(an_range_mat, yr_range_vec)  # n_years_with_data × n_sims

          # Point-estimate annual values for this range
          an_pe_vec <- tapply(an_day[row_sel], proj$year[row_sel], sum)

          for (yr_char in rownames(ann_sims)) {
            yr <- as.integer(yr_char)
            sims_yr <- ann_sims[yr_char, ]
            ci <- quantile(sims_yr, c(0.025, 0.975), na.rm = TRUE)
            an_val <- if (yr_char %in% names(an_pe_vec)) an_pe_vec[yr_char] else 0
            city_results[[length(city_results) + 1]] <- data.table(
              city_code = cc, gcm = gcm_col, ssp = ssp_v,
              age_group = ag, temp_range = rname, year = yr,
              an_est = as.numeric(an_val),
              an_low = as.numeric(ci[1]),
              an_hi  = as.numeric(ci[2]))
          }
        }
      }
    }
  }

  if (length(city_results) > 0) {
    city_dt <- rbindlist(city_results, fill = TRUE)
    fwrite(city_dt, checkpoint_file)
    cat(sprintf("  City %s: %d rows -> %s\n", cc, nrow(city_dt),
      basename(checkpoint_file)))
  }
}

#----- Aggregate all checkpoints
cat("Aggregating checkpoints...\n")
checkpoint_files <- list.files(an_checkpoint, pattern = "\\.csv$", full.names = TRUE)
all_results <- rbindlist(lapply(checkpoint_files, fread), fill = TRUE)
cat("AN complete:", format(nrow(all_results), big.mark = ","), "rows\n")

# Label clearly as climate-only under fixed baseline deaths
all_results[, attribution_type := "climate_only_fixed_baseline"]

fwrite(all_results, file.path(path_out, "ans_annual_all_cities.csv"))
cat("Written:", file.path(path_out, "ans_annual_all_cities.csv"), "\n")