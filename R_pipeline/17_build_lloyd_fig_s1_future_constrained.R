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
merge_input_dir <- trimws(Sys.getenv("MERGE_INPUT_DIR", unset = ""))
merge_output_file <- trimws(Sys.getenv("MERGE_OUTPUT_FILE", unset = ""))
merge_check_file <- trimws(Sys.getenv("MERGE_CHECK_FILE", unset = ""))
city_filter <- trimws(Sys.getenv("CITY_FILTER", unset = ""))
country_filter <- trimws(Sys.getenv("COUNTRY_FILTER", unset = ""))
year_filter <- trimws(Sys.getenv("YEAR_FILTER", unset = ""))
ssp_filter <- trimws(Sys.getenv("SSP_FILTER", unset = ""))
gcm_filter <- trimws(Sys.getenv("GCM_FILTER", unset = ""))
sim_filter <- trimws(Sys.getenv("SIM_FILTER", unset = ""))
n_cores_env <- suppressWarnings(as.integer(trimws(Sys.getenv("N_CORES", unset = ""))))
chunk_tag <- trimws(Sys.getenv("STEP17_CHUNK_TAG", unset = ""))
checkpoint_root <- trimws(Sys.getenv("STEP17_CHECKPOINTS_DIR", unset = ""))
timing_file <- trimws(Sys.getenv("STEP17_TIMING_FILE", unset = ""))
use_checkpoints <- nzchar(chunk_tag) && nzchar(checkpoint_root)

timing_state <- new.env(parent = emptyenv())
timing_state$rows <- list()

append_timing <- function(stage, city_id = NA_character_, slice_id = NA_character_, seconds) {
  timing_state$rows[[length(timing_state$rows) + 1L]] <- data.table(
    stage = stage,
    city_id = as.character(city_id),
    slice_id = as.character(slice_id),
    seconds = as.numeric(seconds)
  )
  invisible(NULL)
}

if (nzchar(merge_input_dir)) {
  if (!nzchar(merge_output_file)) merge_output_file <- output_file
  if (!nzchar(merge_check_file)) merge_check_file <- check_file
  dir.create(dirname(merge_output_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(merge_check_file), recursive = TRUE, showWarnings = FALSE)

  data_files <- sort(list.files(merge_input_dir, pattern = "^lloyd_fig_s1_future_constrained_chunk_[0-9]+\\.csv$", full.names = TRUE))
  check_files <- sort(list.files(merge_input_dir, pattern = "^lloyd_fig_s1_future_constrained_chunk_[0-9]+_checks\\.csv$", full.names = TRUE))

  if (!length(data_files) || !length(check_files)) {
    stop(sprintf("No chunk outputs found in %s", merge_input_dir), call. = FALSE)
  }
  if (length(data_files) != length(check_files)) {
    stop(sprintf("Mismatched chunk data/check files found in %s", merge_input_dir), call. = FALSE)
  }

  merged_data <- rbindlist(lapply(data_files, fread), use.names = TRUE, fill = TRUE)
  merged_checks <- rbindlist(lapply(check_files, fread), use.names = TRUE, fill = TRUE)

  setcolorder(merged_data, c(
    "URAU_CODE", "LABEL", "CNTR_CODE", "cntr_name", "region", "lon", "lat",
    "year", "ssp", "gcm", "sim", "age", "age_label", "pop", "death_baseline",
    "AN_ExtrCold", "AN_ModCold", "AN_ModHeat", "AN_ExtrHeat", "AN_total"
  ))

  fwrite(merged_data, merge_output_file)
  fwrite(merged_checks, merge_check_file)
  message("Merged chunk output into ", merge_output_file)
  message("Merged chunk checks into ", merge_check_file)
  quit(status = 0)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(check_file), recursive = TRUE, showWarnings = FALSE)

if (use_checkpoints) {
  dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(checkpoint_root, chunk_tag, "cities"), recursive = TRUE, showWarnings = FALSE)
}

range_levels <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
age_groups_65plus <- c("65-74", "75-84", "85+")
overshoot_tol <- 1e-8

parse_int_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  tokens <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  as.integer(sub("^V", "", tokens))
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

append_manifest_entry <- function(manifest_file, row) {
  if (!file.exists(manifest_file)) {
    fwrite(row, manifest_file)
  } else {
    existing <- fread(manifest_file, showProgress = FALSE)
    fwrite(rbind(existing, row), manifest_file)
  }
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

future_files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
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
  city_t0 <- Sys.time()
  city_timing_rows <- list()

  append_city_timing <- function(stage, city_id = NA_character_, slice_id = NA_character_, seconds) {
    city_timing_rows[[length(city_timing_rows) + 1L]] <<- data.table(
      stage = stage,
      city_id = as.character(city_id),
      slice_id = as.character(slice_id),
      seconds = as.numeric(seconds)
    )
    invisible(NULL)
  }

  read_t0 <- Sys.time()
  d <- as.data.table(readRDS(file_path))[agegroup %in% age_groups_65plus]
  d[, sim_id := normalize_sim_id(sim)]

  if (!is.null(wanted_years)) d <- d[year %in% wanted_years]
  if (!is.null(wanted_ssps)) d <- d[ssp %in% wanted_ssps]
  if (!is.null(wanted_gcms)) d <- d[gcm %in% wanted_gcms]
  if (!is.null(wanted_sims)) d <- d[sim_id %in% wanted_sims]
  if (!nrow(d)) {
    append_timing("city_total", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), city_t0, units = "secs")))
    append_city_timing("city_total", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), city_t0, units = "secs")))
    return(list(empty = TRUE, city_id = city_id, timing_rows = data.table(
      stage = character(),
      city_id = character(),
      slice_id = character(),
      seconds = numeric()
    )))
  }
  append_timing("readRDS", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), read_t0, units = "secs")))
  append_city_timing("readRDS", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), read_t0, units = "secs")))

  meta_t0 <- Sys.time()
  meta_city <- city_meta[URAU_CODE == city_id][order(age_start)]
  if (nrow(meta_city) != 3) {
    stop(sprintf("Metadata for %s does not contain the three expected 65+ age groups.", city_id))
  }
  append_timing("city_meta", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), meta_t0, units = "secs")))
  append_city_timing("city_meta", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), meta_t0, units = "secs")))

  disagg_t0 <- Sys.time()
  x <- meta_city$age_start
  widths <- meta_city$age_width
  nlast <- tail(widths, 1)
  ages_single <- 65:(65 + pclm_expected_length(x, nlast) - 1)

  pop_single <- pclm_disaggregate_nonnegative(x = x, y = meta_city$pop, nlast = nlast)
  death_single <- pclm_disaggregate_nonnegative(x = x, y = meta_city$death_baseline, nlast = nlast)
  append_timing("city_disagg", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), disagg_t0, units = "secs")))
  append_city_timing("city_disagg", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), disagg_t0, units = "secs")))

  reshape_t0 <- Sys.time()
  grouped <- d[, .(an = sum(an)), by = .(year, ssp, gcm, sim = sim_id, agegroup, range)]
  grouped <- dcast(grouped, year + ssp + gcm + sim + agegroup ~ range, value.var = "an", fill = 0)
  for (range_name in range_levels) {
    if (!range_name %in% names(grouped)) {
      grouped[[range_name]] <- 0
    }
    grouped[[range_name]] <- as.numeric(grouped[[range_name]])
  }
  grouped <- merge(grouped, meta_city[, .(agegroup, age_start, age_width)], by = "agegroup")
  setorder(grouped, year, ssp, gcm, sim, age_start)
  append_timing("city_reshape", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), reshape_t0, units = "secs")))
  append_city_timing("city_reshape", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), reshape_t0, units = "secs")))

  band_starts <- cumsum(c(1L, head(widths, -1L)))
  band_indices <- Map(seq.int, band_starts, band_starts + widths - 1L)

  slices <- unique(grouped[, .(year, ssp, gcm, sim)])
  rows <- vector("list", nrow(slices))
  checks <- vector("list", nrow(slices))

  for (i in seq_len(nrow(slices))) {
    slice <- slices[i]
    slice_t0 <- Sys.time()
    slice_id <- sprintf("%s|%s|%s|%s", slice$year, slice$ssp, slice$gcm, slice$sim)

    grp <- grouped[
      year == slice$year &
      ssp == slice$ssp &
      gcm == slice$gcm &
      sim == slice$sim
    ][order(age_start)]

    precompute_t0 <- Sys.time()
    grp_range_values <- as.data.frame(grp[, range_levels, with = FALSE])
    if (any(as.matrix(grp_range_values) < -1e-10, na.rm = TRUE)) {
      stop(sprintf("Negative grouped AN encountered for %s in constrained prototype.", city_id))
    }

    unconstrained_mat <- sapply(range_levels, function(range_name) {
      pclm_disaggregate_nonnegative(x = x, y = grp[[range_name]], nlast = nlast)
    })
    unconstrained_mat <- as.matrix(unconstrained_mat)
    colnames(unconstrained_mat) <- range_levels
    unconstrained_total <- rowSums(unconstrained_mat)
    append_timing("slice_precompute", city_id = city_id, slice_id = slice_id, seconds = as.numeric(difftime(Sys.time(), precompute_t0, units = "secs")))
    append_city_timing("slice_precompute", city_id = city_id, slice_id = slice_id, seconds = as.numeric(difftime(Sys.time(), precompute_t0, units = "secs")))

    ipf_t0 <- Sys.time()
    constrained_mat <- matrix(0, nrow = length(ages_single), ncol = length(range_levels))
    colnames(constrained_mat) <- range_levels

    for (band_idx in seq_along(band_indices)) {
      idx <- band_indices[[band_idx]]
      seed_band <- unconstrained_mat[idx, , drop = FALSE]
      cap_band <- death_single[idx]
      target_cols <- as.numeric(grp_range_values[band_idx, range_levels])
      target_rows <- redistribute_with_caps(rowSums(seed_band), cap_band, target_total = sum(target_cols))
      constrained_mat[idx, ] <- ipf_with_margins(seed_band, target_rows, target_cols)
    }

    constrained_total <- rowSums(constrained_mat)
    append_timing("slice_ipf", city_id = city_id, slice_id = slice_id, seconds = as.numeric(difftime(Sys.time(), ipf_t0, units = "secs")))
    append_city_timing("slice_ipf", city_id = city_id, slice_id = slice_id, seconds = as.numeric(difftime(Sys.time(), ipf_t0, units = "secs")))

    output_t0 <- Sys.time()
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
    append_timing("slice_output", city_id = city_id, slice_id = slice_id, seconds = as.numeric(difftime(Sys.time(), output_t0, units = "secs")))
    append_city_timing("slice_output", city_id = city_id, slice_id = slice_id, seconds = as.numeric(difftime(Sys.time(), output_t0, units = "secs")))
    append_timing("slice_total", city_id = city_id, slice_id = slice_id, seconds = as.numeric(difftime(Sys.time(), slice_t0, units = "secs")))
    append_city_timing("slice_total", city_id = city_id, slice_id = slice_id, seconds = as.numeric(difftime(Sys.time(), slice_t0, units = "secs")))
  }

  append_timing("city_total", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), city_t0, units = "secs")))
  append_city_timing("city_total", city_id = city_id, seconds = as.numeric(difftime(Sys.time(), city_t0, units = "secs")))

  timing_rows <- if (length(city_timing_rows)) rbindlist(city_timing_rows, fill = TRUE) else data.table(
    stage = character(),
    city_id = character(),
    slice_id = character(),
    seconds = numeric()
  )

  list(empty = FALSE, city_id = city_id, data = rbindlist(rows), checks = rbindlist(checks), timing_rows = timing_rows)
}

if (use_checkpoints) {
  manifest_file <- file.path(checkpoint_root, "manifest.csv")
  if (file.exists(manifest_file)) {
    manifest <- tryCatch(fread(manifest_file, showProgress = FALSE), error = function(e) data.table(
      chunk_tag = character(),
      output_file = character(),
      check_file = character(),
      completed_at = character(),
      completed = logical(),
      cities = integer(),
      rows = integer()
    ))
  } else {
    manifest <- data.table(
      chunk_tag = character(),
      output_file = character(),
      check_file = character(),
      completed_at = character(),
      completed = logical(),
      cities = integer(),
      rows = integer()
    )
  }

  chunk_tag_value <- chunk_tag
  if (nrow(manifest[chunk_tag == chunk_tag_value & completed]) > 0L && file.exists(output_file) && file.exists(check_file)) {
    message(sprintf("[17] Chunk %s already completed; skipping recomputation.", chunk_tag))
    quit(status = 0)
  }
}

city_results <- vector("list", length(future_files))
city_checks <- vector("list", length(future_files))
processed_count <- 0L

parallel_cores <- if (!is.na(n_cores_env) && n_cores_env > 1L) n_cores_env else 1L
if (parallel_cores > 1L) {
  message(sprintf("[17] chunk %s: using city-level parallelism with %d cores", chunk_tag, parallel_cores))
}

city_plan <- vector("list", length(future_files))
for (idx in seq_along(future_files)) {
  file_path <- future_files[[idx]]
  city_id <- sub("\\.rds$", "", basename(file_path))
  checkpoint_file <- if (use_checkpoints) file.path(checkpoint_root, chunk_tag, "cities", paste0(city_id, ".rds")) else ""
  city_plan[[idx]] <- list(file_path = file_path, city_id = city_id, checkpoint_file = checkpoint_file)
}

run_city_job <- function(job) {
  file_path <- job$file_path
  city_id <- job$city_id
  checkpoint_file <- job$checkpoint_file

  if (use_checkpoints && file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    message(sprintf("[17] chunk %s city %s: resumed from checkpoint", chunk_tag, city_id))
    return(list(
      empty = FALSE,
      city_id = city_id,
      data = checkpoint$data,
      checks = checkpoint$checks,
      resumed = TRUE,
      timing_rows = data.table(
        stage = character(),
        city_id = character(),
        slice_id = character(),
        seconds = numeric()
      )
    ))
  }

  message(sprintf("[17] chunk %s city %s: processing", chunk_tag, city_id))
  city_result <- process_city(file_path)
  if (isTRUE(city_result$empty)) {
    message(sprintf("[17] chunk %s city %s: empty after filtering", chunk_tag, city_id))
    return(list(empty = TRUE, city_id = city_id, resumed = FALSE))
  }

  if (use_checkpoints) {
    saveRDS(list(data = city_result$data, checks = city_result$checks), checkpoint_file)
  }

  message(sprintf("[17] chunk %s city %s: completed", chunk_tag, city_id))
  list(empty = FALSE, city_id = city_id, data = city_result$data, checks = city_result$checks, resumed = FALSE, timing_rows = city_result$timing_rows)
}

if (parallel_cores > 1L) {
  worker_results <- mclapply(city_plan, run_city_job, mc.cores = parallel_cores)
} else {
  worker_results <- lapply(city_plan, run_city_job)
}

for (idx in seq_along(worker_results)) {
  result <- worker_results[[idx]]
  if (isTRUE(result$empty)) next
  city_results[[idx]] <- result$data
  city_checks[[idx]] <- result$checks
  if (nrow(result$timing_rows)) {
    timing_state$rows[[length(timing_state$rows) + 1L]] <- result$timing_rows
  }
  processed_count <- processed_count + 1L
}

city_results <- city_results[!vapply(city_results, is.null, logical(1))]
city_checks <- city_checks[!vapply(city_checks, is.null, logical(1))]

if (!length(city_results)) {
  filter_msg <- c(
    if (nzchar(country_filter)) sprintf("COUNTRY_FILTER=%s", country_filter) else NULL,
    if (nzchar(city_filter)) sprintf("CITY_FILTER=%s", city_filter) else NULL,
    if (nzchar(year_filter)) sprintf("YEAR_FILTER=%s", year_filter) else NULL,
    if (nzchar(ssp_filter)) sprintf("SSP_FILTER=%s", ssp_filter) else NULL,
    if (nzchar(gcm_filter)) sprintf("GCM_FILTER=%s", gcm_filter) else NULL,
    if (nzchar(sim_filter)) sprintf("SIM_FILTER=%s", sim_filter) else NULL
  )
  filter_text <- if (length(filter_msg)) paste(filter_msg, collapse = ", ") else "(none)"

  stop(
    sprintf(
      paste0(
        "All selected future files were filtered out; nothing to write.\n",
        "  files scheduled: %d\n",
        "  active filters: %s"
      ),
      length(future_files), filter_text
    ),
    call. = FALSE
  )
}

lloyd_s1_future <- rbindlist(city_results, use.names = TRUE)
checks <- rbindlist(city_checks, use.names = TRUE)

setcolorder(lloyd_s1_future, c(
  "URAU_CODE", "LABEL", "CNTR_CODE", "cntr_name", "region", "lon", "lat",
  "year", "ssp", "gcm", "sim", "age", "age_label", "pop", "death_baseline",
  "AN_ExtrCold", "AN_ModCold", "AN_ModHeat", "AN_ExtrHeat", "AN_total"
))

write_t0 <- Sys.time()
fwrite(lloyd_s1_future, output_file)
fwrite(checks, check_file)
append_timing("writing", seconds = as.numeric(difftime(Sys.time(), write_t0, units = "secs")))

if (use_checkpoints) {
  manifest_entry <- data.table(
    chunk_tag = chunk_tag,
    output_file = output_file,
    check_file = check_file,
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"),
    completed = TRUE,
    cities = length(city_results),
    rows = nrow(lloyd_s1_future)
  )
  append_manifest_entry(manifest_file, manifest_entry)
}

message(sprintf("[17] chunk %s completed: wrote %d city outputs and %d rows", chunk_tag, length(city_results), nrow(lloyd_s1_future)))
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

if (length(timing_state$rows)) {
  timing_dt <- rbindlist(timing_state$rows, fill = TRUE)
  phase_summary <- timing_dt[stage != "city_total" & stage != "slice_total", .(seconds = sum(seconds, na.rm = TRUE)), by = stage]
  setorder(phase_summary, -seconds)
  city_summary <- timing_dt[stage == "city_total", .(seconds = sum(seconds, na.rm = TRUE)), by = city_id]
  setorder(city_summary, -seconds)
  slice_summary <- timing_dt[stage == "slice_total", .(seconds = sum(seconds, na.rm = TRUE)), by = .(city_id, slice_id)]
  setorder(slice_summary, -seconds)

  message("[17] Timing summary by phase:")
  print(phase_summary)
  message("[17] Timing summary by city:")
  print(city_summary)
  message("[17] Timing summary by slice:")
  print(slice_summary)

  if (nzchar(timing_file)) {
    dir.create(dirname(timing_file), recursive = TRUE, showWarnings = FALSE)
    summary_dt <- rbind(
      city_summary[, .(level = "city", city_id, slice_id = NA_character_, seconds)],
      slice_summary[, .(level = "slice", city_id, slice_id, seconds)]
    )
    fwrite(summary_dt, timing_file)
    message("[17] Timing breakdown written to ", timing_file)
  }
}
