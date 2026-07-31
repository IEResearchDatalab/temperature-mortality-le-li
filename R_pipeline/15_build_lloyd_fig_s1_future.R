#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ungroup)
  library(parallel)
})

source("R_pipeline/functions/pclm_utils.R")

message("\n[15] Building Lloyd-style single-age future AN table...")

input_dir <- trimws(Sys.getenv("INPUT_DIR", unset = "temp_results_future"))
output_file <- trimws(Sys.getenv("OUTPUT_FILE", unset = "results/le_li_input/lloyd_fig_s1_future_city_year_single_age.csv"))
check_file <- trimws(Sys.getenv("CHECK_FILE", unset = "results/le_li_input/lloyd_fig_s1_future_checks.csv"))
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
  nlast <- tail(meta_city$age_width, 1)
  ages_single <- 65:(65 + pclm_expected_length(x, nlast) - 1)

  pop_single <- pclm_disaggregate_nonnegative(x = x, y = meta_city$pop, nlast = nlast)
  death_single <- pclm_disaggregate_nonnegative(x = x, y = meta_city$death_baseline, nlast = nlast)

  grouped <- d[, .(an = sum(an)), by = .(year, ssp, gcm, sim = sim_id, agegroup, range)]
  grouped <- dcast(grouped, year + ssp + gcm + sim + agegroup ~ range, value.var = "an", fill = 0)
  grouped <- merge(grouped, meta_city[, .(agegroup, age_start)], by = "agegroup")
  setorder(grouped, year, ssp, gcm, sim, age_start)

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
      AN_ExtrCold = pclm_disaggregate_signed(x = x, y = grp$ExtrCold, nlast = nlast),
      AN_ModCold = pclm_disaggregate_signed(x = x, y = grp$ModCold, nlast = nlast),
      AN_ModHeat = pclm_disaggregate_signed(x = x, y = grp$ModHeat, nlast = nlast),
      AN_ExtrHeat = pclm_disaggregate_signed(x = x, y = grp$ExtrHeat, nlast = nlast)
    )
    single[, AN_total := AN_ExtrCold + AN_ModCold + AN_ModHeat + AN_ExtrHeat]
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

message("Saved Lloyd-style future table to ", output_file)
message("Saved additivity checks to ", check_file)
message(sprintf(
  "Rows: %d | Cities: %d | Years: %d | SSPs: %d | GCMs: %d | Sims: %d | Any NA rows in checks: %s",
  nrow(lloyd_s1_future),
  uniqueN(lloyd_s1_future$URAU_CODE),
  uniqueN(lloyd_s1_future$year),
  uniqueN(lloyd_s1_future$ssp),
  uniqueN(lloyd_s1_future$gcm),
  uniqueN(lloyd_s1_future$sim),
  any(checks$any_na)
))