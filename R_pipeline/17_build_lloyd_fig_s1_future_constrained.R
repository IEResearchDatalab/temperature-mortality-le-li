#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ungroup)
  library(parallel)
})

source("R_pipeline/functions/pclm_utils.R")

message("\n[17] Building constrained Lloyd-style single-age future AN table...")

input_dir <- trimws(Sys.getenv("INPUT_DIR", unset = "temp_results_future"))
output_file <- trimws(Sys.getenv("OUTPUT_FILE", unset = "results/le_li_input/lloyd_fig_s1_future_constrained_city_year_single_age.csv"))
check_file <- trimws(Sys.getenv("CHECK_FILE", unset = "results/le_li_input/lloyd_fig_s1_future_constrained_checks.csv"))
city_filter <- trimws(Sys.getenv("CITY_FILTER", unset = ""))
country_filter <- trimws(Sys.getenv("COUNTRY_FILTER", unset = ""))
year_filter <- trimws(Sys.getenv("YEAR_FILTER", unset = ""))
ssp_filter <- trimws(Sys.getenv("SSP_FILTER", unset = ""))
gcm_filter <- trimws(Sys.getenv("GCM_FILTER", unset = ""))
sim_filter <- trimws(Sys.getenv("SIM_FILTER", unset = ""))
n_cores_env <- suppressWarnings(as.integer(trimws(Sys.getenv("N_CORES", unset = ""))))

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(check_file), recursive = TRUE, showWarnings = FALSE)

range_levels <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
age_groups_65plus <- c("65-74", "75-84", "85+")
overshoot_tol <- 1e-8

parse_int_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

parse_chr_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

normalize_sim_id <- function(x) {
  x_chr <- as.character(x)
  suppressWarnings(as.integer(sub("^V", "", x_chr)))
}

redistribute_with_caps <- function(raw, cap, target_total = NULL, tol = 1e-10, max_iter = 1000L) {
  raw <- as.numeric(raw)
  cap <- as.numeric(cap)
  if (is.null(target_total)) {
    target_total <- sum(raw, na.rm = TRUE)
  }

  if (target_total <= tol) return(rep(0, length(raw)))
  if (target_total - sum(cap, na.rm = TRUE) > tol) {
    stop("Target total exceeds available capacity.")
  }

  raw_sum <- sum(raw, na.rm = TRUE)
  if (raw_sum > tol) {
    raw <- raw * (target_total / raw_sum)
  } else {
    raw <- rep(target_total / length(raw), length(raw))
  }

  adjusted <- pmin(raw, cap)
  deficit <- target_total - sum(adjusted, na.rm = TRUE)
  iter <- 0L

  while (deficit > tol && iter < max_iter) {
    room <- pmax(cap - adjusted, 0)
    room_sum <- sum(room, na.rm = TRUE)
    if (room_sum <= tol) break
    adjusted <- pmin(adjusted + deficit * room / room_sum, cap)
    deficit <- target_total - sum(adjusted, na.rm = TRUE)
    iter <- iter + 1L
  }

  if (abs(target_total - sum(adjusted, na.rm = TRUE)) > 1e-8) {
    stop("Failed to redistribute constrained totals within tolerance.")
  }

  adjusted
}

ipf_with_margins <- function(seed, target_rows, target_cols, tol = 1e-8, max_iter = 5000L) {
  seed <- as.matrix(seed)
  target_rows <- as.numeric(target_rows)
  target_cols <- as.numeric(target_cols)

  if (abs(sum(target_rows) - sum(target_cols)) > tol) {
    stop("Row and column targets do not sum to the same total.")
  }

  n_rows <- length(target_rows)
  n_cols <- length(target_cols)
  out <- matrix(0, nrow = n_rows, ncol = n_cols)

  active_rows <- target_rows > tol
  active_cols <- target_cols > tol
  if (!any(active_rows) || !any(active_cols)) return(out)

  mat <- pmax(seed[active_rows, active_cols, drop = FALSE], 1e-8)
  row_targets <- target_rows[active_rows]
  col_targets <- target_cols[active_cols]

  for (iter in seq_len(max_iter)) {
    row_sums <- rowSums(mat)
    mat <- mat * (row_targets / row_sums)

    col_sums <- colSums(mat)
    mat <- sweep(mat, 2L, col_targets / col_sums, `*`)

    row_err <- max(abs(rowSums(mat) - row_targets))
    col_err <- max(abs(colSums(mat) - col_targets))
    if (max(row_err, col_err) <= tol) break

    if (iter == max_iter) {
      stop("IPF did not converge within the iteration limit.")
    }
  }

  out[active_rows, active_cols] <- mat
  out
}

wanted_years <- parse_int_filter(year_filter)
wanted_ssps <- parse_int_filter(ssp_filter)
wanted_sims <- parse_int_filter(sim_filter)
wanted_gcms <- parse_chr_filter(gcm_filter)
wanted_countries <- parse_chr_filter(country_filter)

city_meta <- fread("data/city_results.csv")
city_meta <- unique(city_meta[agegroup %in% age_groups_65plus, .(
  URAU_CODE,
  LABEL,
  CNTR_CODE,
  cntr_name,
  region,
  lon,
  lat,
  agegroup,
  age_start = fifelse(agegroup == "65-74", 65L, fifelse(agegroup == "75-84", 75L, 85L)),
  age_width = fifelse(agegroup == "85+", 16L, 10L),
  pop = agepop,
  death_baseline = death
)])

future_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
if (!is.null(wanted_countries)) {
  wanted_country_cities <- unique(city_meta[CNTR_CODE %in% wanted_countries, URAU_CODE])
  future_files <- future_files[sub("\\.rds$", "", basename(future_files)) %in% wanted_country_cities]
}
if (nzchar(city_filter)) {
  wanted_cities <- parse_chr_filter(city_filter)
  future_files <- future_files[sub("\\.rds$", "", basename(future_files)) %in% wanted_cities]
}
if (!length(future_files)) {
  stop(sprintf("No matching future files found in %s", input_dir), call. = FALSE)
}

process_city <- function(file_path) {
  city_id <- sub("\\.rds$", "", basename(file_path))
  d <- as.data.table(readRDS(file_path))[agegroup %in% age_groups_65plus]
  d[, sim_id := normalize_sim_id(sim)]

  if (!is.null(wanted_years)) d <- d[year %in% wanted_years]
  if (!is.null(wanted_ssps)) d <- d[ssp %in% wanted_ssps]
  if (!is.null(wanted_gcms)) d <- d[gcm %in% wanted_gcms]
  if (!is.null(wanted_sims)) d <- d[sim_id %in% wanted_sims]
  if (!nrow(d)) return(NULL)

  meta_city <- city_meta[URAU_CODE == city_id][order(age_start)]
  if (nrow(meta_city) != 3) {
    stop(sprintf("Metadata for %s does not contain the three expected 65+ age groups.", city_id))
  }

  x <- meta_city$age_start
  widths <- meta_city$age_width
  nlast <- tail(widths, 1)
  ages_single <- 65:(65 + pclm_expected_length(x, nlast) - 1)

  pop_single <- pclm_disaggregate_nonnegative(x = x, y = meta_city$pop, nlast = nlast)
  death_single <- pclm_disaggregate_nonnegative(x = x, y = meta_city$death_baseline, nlast = nlast)

  grouped <- d[, .(an = sum(an)), by = .(year, ssp, gcm, sim = sim_id, agegroup, range)]
  grouped <- dcast(grouped, year + ssp + gcm + sim + agegroup ~ range, value.var = "an", fill = 0)
  grouped <- merge(grouped, meta_city[, .(agegroup, age_start, age_width)], by = "agegroup")
  setorder(grouped, year, ssp, gcm, sim, age_start)

  band_starts <- cumsum(c(1L, head(widths, -1L)))
  band_indices <- Map(seq.int, band_starts, band_starts + widths - 1L)

  slices <- unique(grouped[, .(year, ssp, gcm, sim)])
  rows <- vector("list", nrow(slices))
  checks <- vector("list", nrow(slices))

  for (i in seq_len(nrow(slices))) {
    slice <- slices[i]
    grp <- grouped[
      year == slice$year &
      ssp == slice$ssp &
      gcm == slice$gcm &
      sim == slice$sim
    ][order(age_start)]

    if (any(as.matrix(grp[, ..range_levels]) < -1e-10, na.rm = TRUE)) {
      stop(sprintf("Negative grouped AN encountered for %s in constrained prototype.", city_id))
    }

    unconstrained_mat <- sapply(range_levels, function(range_name) {
      pclm_disaggregate_nonnegative(x = x, y = grp[[range_name]], nlast = nlast)
    })
    unconstrained_mat <- as.matrix(unconstrained_mat)
    colnames(unconstrained_mat) <- range_levels
    unconstrained_total <- rowSums(unconstrained_mat)

    constrained_mat <- matrix(0, nrow = length(ages_single), ncol = length(range_levels))
    colnames(constrained_mat) <- range_levels

    for (band_idx in seq_along(band_indices)) {
      idx <- band_indices[[band_idx]]
      seed_band <- unconstrained_mat[idx, , drop = FALSE]
      cap_band <- death_single[idx]
      target_cols <- as.numeric(grp[band_idx, ..range_levels])
      target_rows <- redistribute_with_caps(rowSums(seed_band), cap_band, target_total = sum(target_cols))
      constrained_mat[idx, ] <- ipf_with_margins(seed_band, target_rows, target_cols)
    }

    constrained_total <- rowSums(constrained_mat)

    single <- data.table(
      URAU_CODE = city_id,
      LABEL = meta_city$LABEL[1],
      CNTR_CODE = meta_city$CNTR_CODE[1],
      cntr_name = meta_city$cntr_name[1],
      region = meta_city$region[1],
      lon = meta_city$lon[1],
      lat = meta_city$lat[1],
      year = slice$year,
      ssp = slice$ssp,
      gcm = slice$gcm,
      sim = slice$sim,
      age = ages_single,
      age_label = fifelse(ages_single == 100L, "100+", as.character(ages_single)),
      pop = pop_single,
      death_baseline = death_single,
      AN_ExtrCold = constrained_mat[, "ExtrCold"],
      AN_ModCold = constrained_mat[, "ModCold"],
      AN_ModHeat = constrained_mat[, "ModHeat"],
      AN_ExtrHeat = constrained_mat[, "ExtrHeat"],
      AN_total = constrained_total
    )
    rows[[i]] <- single

    checks[[i]] <- data.table(
      URAU_CODE = city_id,
      year = slice$year,
      ssp = slice$ssp,
      gcm = slice$gcm,
      sim = slice$sim,
      diff_ExtrCold = sum(single$AN_ExtrCold) - sum(grp$ExtrCold),
      diff_ModCold = sum(single$AN_ModCold) - sum(grp$ModCold),
      diff_ModHeat = sum(single$AN_ModHeat) - sum(grp$ModHeat),
      diff_ExtrHeat = sum(single$AN_ExtrHeat) - sum(grp$ExtrHeat),
      max_pre_overshoot = max(unconstrained_total - death_single),
      max_post_overshoot = max(constrained_total - death_single),
      n_pre_overshoot = sum(unconstrained_total > death_single + overshoot_tol),
      n_post_overshoot = sum(constrained_total > death_single + overshoot_tol),
      any_na = anyNA(single[, .(AN_ExtrCold, AN_ModCold, AN_ModHeat, AN_ExtrHeat, pop, death_baseline)])
    )
  }

  list(data = rbindlist(rows), checks = rbindlist(checks))
}

num_cores <- if (is.na(n_cores_env) || n_cores_env < 1L) {
  min(8L, max(1L, detectCores() - 1L))
} else {
  n_cores_env
}

results <- mclapply(future_files, process_city, mc.cores = num_cores)
results <- Filter(Negate(is.null), results)

if (!length(results)) {
  stop("All selected future files were filtered out; nothing to write.", call. = FALSE)
}

lloyd_s1_future <- rbindlist(lapply(results, `[[`, "data"), use.names = TRUE)
checks <- rbindlist(lapply(results, `[[`, "checks"), use.names = TRUE)

setcolorder(lloyd_s1_future, c(
  "URAU_CODE", "LABEL", "CNTR_CODE", "cntr_name", "region", "lon", "lat",
  "year", "ssp", "gcm", "sim", "age", "age_label", "pop", "death_baseline",
  "AN_ExtrCold", "AN_ModCold", "AN_ModHeat", "AN_ExtrHeat", "AN_total"
))

fwrite(lloyd_s1_future, output_file)
fwrite(checks, check_file)

message("Saved constrained Lloyd-style future table to ", output_file)
message("Saved constrained additivity/cap checks to ", check_file)
message(sprintf(
  "Rows: %d | Cities: %d | Years: %d | SSPs: %d | GCMs: %d | Sims: %d | Any NA rows in checks: %s | Max post-overshoot: %g",
  nrow(lloyd_s1_future),
  uniqueN(lloyd_s1_future$URAU_CODE),
  uniqueN(lloyd_s1_future$year),
  uniqueN(lloyd_s1_future$ssp),
  uniqueN(lloyd_s1_future$gcm),
  uniqueN(lloyd_s1_future$sim),
  any(checks$any_na),
  max(checks$max_post_overshoot)
))