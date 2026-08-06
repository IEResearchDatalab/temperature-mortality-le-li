#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

message("\n[21] Collapsing future simulation draws to point-estimate RDS files...")

input_dir <- trimws(Sys.getenv("INPUT_DIR", unset = "temp_results_future"))
output_dir <- trimws(Sys.getenv("OUTPUT_DIR", unset = "temp_results_future_collapsed"))
city_filter <- trimws(Sys.getenv("CITY_FILTER", unset = ""))
country_filter <- trimws(Sys.getenv("COUNTRY_FILTER", unset = ""))
year_filter <- trimws(Sys.getenv("YEAR_FILTER", unset = ""))
ssp_filter <- trimws(Sys.getenv("SSP_FILTER", unset = ""))
gcm_filter <- trimws(Sys.getenv("GCM_FILTER", unset = ""))

parse_int_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

parse_chr_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

wanted_years <- parse_int_filter(year_filter)
wanted_ssps <- parse_int_filter(ssp_filter)
wanted_gcms <- parse_chr_filter(gcm_filter)
wanted_cities <- parse_chr_filter(city_filter)
wanted_countries <- parse_chr_filter(country_filter)

if (!dir.exists(input_dir)) {
  stop(sprintf("INPUT_DIR not found: %s", input_dir), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

city_meta <- fread("data/city_results.csv", select = c("URAU_CODE", "CNTR_CODE"))
city_meta <- unique(city_meta)

future_files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
if (!length(future_files)) {
  stop(sprintf("No .rds files found in %s", input_dir), call. = FALSE)
}

if (!is.null(wanted_countries)) {
  keep_city <- unique(city_meta[CNTR_CODE %in% wanted_countries, URAU_CODE])
  future_files <- future_files[sub("\\.rds$", "", basename(future_files)) %in% keep_city]
}
if (!is.null(wanted_cities)) {
  future_files <- future_files[sub("\\.rds$", "", basename(future_files)) %in% wanted_cities]
}

if (!length(future_files)) {
  stop("No matching city files after filters.", call. = FALSE)
}

summary_rows <- vector("list", length(future_files))

for (i in seq_along(future_files)) {
  file_path <- future_files[[i]]
  city_id <- sub("\\.rds$", "", basename(file_path))
  t0 <- Sys.time()

  d <- as.data.table(readRDS(file_path))

  if (!is.null(wanted_years)) d <- d[year %in% wanted_years]
  if (!is.null(wanted_ssps)) d <- d[ssp %in% wanted_ssps]
  if (!is.null(wanted_gcms)) d <- d[gcm %in% wanted_gcms]

  if (!nrow(d)) {
    message(sprintf("[21] %s: empty after filtering, skipping write", city_id))
    summary_rows[[i]] <- data.table(
      URAU_CODE = city_id,
      rows_in = 0L,
      rows_out = 0L,
      unique_sim_in = 0L,
      unique_sim_out = 0L,
      elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    next
  }

  rows_in <- nrow(d)
  unique_sim_in <- uniqueN(d$sim)

  collapsed <- d[, .(an = mean(an, na.rm = TRUE)), by = .(year, range, ssp, gcm, agegroup)]
  collapsed[, sim := 0L]
  setcolorder(collapsed, c("sim", "an", "year", "range", "ssp", "gcm", "agegroup"))

  out_file <- file.path(output_dir, basename(file_path))
  saveRDS(collapsed, out_file)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf(
    "[21] %s: rows_in=%d sims_in=%d rows_out=%d -> %s (%.2fs)",
    city_id,
    rows_in,
    unique_sim_in,
    nrow(collapsed),
    out_file,
    elapsed
  ))

  summary_rows[[i]] <- data.table(
    URAU_CODE = city_id,
    rows_in = rows_in,
    rows_out = nrow(collapsed),
    unique_sim_in = unique_sim_in,
    unique_sim_out = uniqueN(collapsed$sim),
    elapsed_sec = elapsed
  )
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
summary_file <- file.path(output_dir, "collapse_summary.csv")
fwrite(summary_dt, summary_file)

message(sprintf(
  "[21] Completed. Cities processed=%d | Cities with output=%d | Total rows out=%d",
  nrow(summary_dt),
  sum(summary_dt$rows_out > 0),
  sum(summary_dt$rows_out)
))
message("[21] Summary written to ", summary_file)
