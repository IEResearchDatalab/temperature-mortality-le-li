#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

default_test_cities <- c("AT001C", "ES001C", "UK001C", "BG001C", "SE001C")
expected_ranges <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")

args <- commandArgs(trailingOnly = TRUE)
component <- if (length(args) >= 1) tolower(args[[1]]) else "all"
test_cities <- if (length(args) >= 2) trimws(strsplit(args[[2]], ",", fixed = TRUE)[[1]]) else default_test_cities
input_dir <- if (length(args) >= 3) args[[3]] else if (dir.exists("temp_results_baseline")) "temp_results_baseline" else "temp_results"

fail <- function(message) {
  stop(message, call. = FALSE)
}

pass <- function(message) {
  cat(sprintf("PASS: %s\n", message))
}

skip <- function(message) {
  cat(sprintf("SKIP: %s\n", message))
}

leap_days <- function(year) {
  ifelse((year %% 400L == 0L) | (year %% 4L == 0L & year %% 100L != 0L), 366L, 365L)
}

sample_files <- function(dir_path, cities) {
  if (!dir.exists(dir_path)) {
    fail(sprintf("Input directory '%s' does not exist.", dir_path))
  }

  files <- file.path(dir_path, paste0(cities, ".rds"))
  files <- files[file.exists(files)]
  if (!length(files)) {
    fail(sprintf("No sample RDS files found in '%s' for the requested cities.", dir_path))
  }
  files
}

validate_structure <- function(dt, city_id) {
  required_cols <- c("sim", "an", "year", "range", "ssp", "gcm", "agegroup")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols)) {
    fail(sprintf("%s is missing required columns: %s", city_id, paste(missing_cols, collapse = ", ")))
  }

  ranges_found <- sort(unique(dt$range))
  if (!all(expected_ranges %in% ranges_found)) {
    fail(sprintf("%s is missing expected ranges. Found: %s", city_id, paste(ranges_found, collapse = ", ")))
  }
}

baseline_point_estimate <- function(dt_city) {
  point <- copy(dt_city[gcm == "ERA5"])
  if (!nrow(point)) {
    fail("Baseline exact check requires ERA5 rows.")
  }
  if (!any(point$sim == 0L, na.rm = TRUE)) {
    fail("Baseline exact check requires point-estimate rows with sim == 0.")
  }

  point <- point[sim == 0L, .(an = sum(an)), by = .(agegroup, year, range)]
  full_index <- CJ(
    agegroup = unique(point$agegroup),
    year = sort(unique(point$year)),
    range = expected_ranges,
    unique = TRUE
  )
  point <- point[full_index, on = .(agegroup, year, range)]
  point[is.na(an), an := 0]
  point[, year_days := leap_days(year)]
  point[, .(an_est = sum(an * year_days) / sum(year_days)), by = .(agegroup, range)]
}

test_script03_output <- function(cities, dir_path) {
  files <- sample_files(dir_path, cities)
  checked <- character(0)

  for (f in files) {
    city_id <- sub("\\.rds$", "", basename(f))
    dt <- as.data.table(readRDS(f))
    validate_structure(dt, city_id)
    checked <- c(checked, city_id)

    if (basename(dir_path) == "temp_results_baseline") {
      ref <- fread("references/2025-masselot-zenodo/results/cityage.csv")
      est <- baseline_point_estimate(dt)
      est_wide <- dcast(est, agegroup ~ range, value.var = "an_est", fill = 0)
      est_wide[, `:=`(
        our_cold = ExtrCold + ModCold,
        our_heat = ModHeat + ExtrHeat
      )]

      ref_city <- ref[URAU_CODE == city_id, .(
        agegroup,
        masselot_cold = excess_cold_est,
        masselot_heat = excess_heat_est
      )]
      cmp <- merge(est_wide, ref_city, by = "agegroup")
      cmp[, `:=`(
        diff_cold = abs(our_cold - masselot_cold),
        diff_heat = abs(our_heat - masselot_heat)
      )]

      if (max(cmp$diff_cold, cmp$diff_heat) > 1e-8) {
        fail(sprintf("Baseline exact check failed for %s (max diff %.6g).", city_id, max(cmp$diff_cold, cmp$diff_heat)))
      }
    }
  }

  pass(sprintf("script03_output checked %d city files in %s", length(checked), dir_path))
}

test_script04_aggregation <- function(cities, dir_path) {
  files <- sample_files(dir_path, cities)
  parts <- lapply(files, function(f) {
    city_id <- sub("\\.rds$", "", basename(f))
    dt <- as.data.table(readRDS(f))
    validate_structure(dt, city_id)
    dt[, .(an = sum(an)), by = .(year, range, sim, ssp, gcm)]
  })

  city_year <- rbindlist(parts, idcol = "city_index")
  city_year[, decade := (year %/% 10) * 10]
  decade_summary <- city_year[, .(
    an_mean = mean(an),
    an_p2_5 = quantile(an, 0.025),
    an_p97_5 = quantile(an, 0.975)
  ), by = .(decade, range, ssp, gcm)]

  if (!nrow(decade_summary)) {
    fail("Script 04 sample aggregation produced no rows.")
  }

  if (!all(expected_ranges %in% unique(decade_summary$range))) {
    fail("Script 04 sample aggregation is missing one or more temperature ranges.")
  }

  pass(sprintf("script04_aggregation produced %d summary rows from %d sample cities", nrow(decade_summary), length(files)))
}

test_script05_tables <- function(cities, dir_path) {
  files <- sample_files(dir_path, cities)
  city_meta <- fread("data/city_results.csv")[, .(
    URAU_CODE,
    cntr_name,
    agegroup,
    pop = agepop,
    deaths_baseline = death
  )]

  parts <- lapply(files, function(f) {
    city_id <- sub("\\.rds$", "", basename(f))
    dt <- as.data.table(readRDS(f))
    validate_structure(dt, city_id)
    meta_v <- city_meta[URAU_CODE == city_id]
    if (!nrow(meta_v)) {
      fail(sprintf("Missing city metadata for %s", city_id))
    }

    dt[, decade := (year %/% 10) * 10]
    out <- dt[, .(an = sum(an)), by = .(decade, range, ssp, gcm, agegroup, sim)]
    out[, pop_baseline := meta_v$pop[match(agegroup, meta_v$agegroup)]]
    out[, deaths_baseline := meta_v$deaths_baseline[match(agegroup, meta_v$agegroup)]]
    out[, country := meta_v$cntr_name[1]]
    out
  })

  sample_country <- rbindlist(parts)
  if (!nrow(sample_country)) {
    fail("Script 05 sample preparation produced no rows.")
  }

  table_s6 <- sample_country[, .(
    an = sum(an),
    pop = sum(pop_baseline),
    deaths = sum(deaths_baseline)
  ), by = .(country, decade, ssp, gcm, sim, range)]
  table_s6[, type := ifelse(grepl("Cold", range), "Cold", "Heat")]
  table_s6 <- table_s6[, .(
    AN = mean(an),
    AN_low = quantile(an, 0.025),
    AN_hi = quantile(an, 0.975)
  ), by = .(country, decade, ssp, type)]

  if (!nrow(table_s6)) {
    fail("Script 05 sample table logic produced no summary rows.")
  }

  pass(sprintf("script05_tables produced %d summary rows from %d sample cities", nrow(table_s6), length(files)))
}

test_pclm <- function(cities) {
  if (!requireNamespace("ungroup", quietly = TRUE)) {
    skip("pclm test skipped because package 'ungroup' is not available")
    return(invisible(NULL))
  }

  city_id <- cities[[1]]
  city_meta <- fread("data/city_results.csv")[
    URAU_CODE == city_id & agegroup %in% c("65-74", "75-84", "85+"),
    .(agegroup, pop = agepop, death = death)
  ]

  if (nrow(city_meta) != 3) {
    fail(sprintf("pclm test needs three 65+ age groups for %s", city_id))
  }

  city_meta[agegroup == "65-74", age_start := 65]
  city_meta[agegroup == "75-84", age_start := 75]
  city_meta[agegroup == "85+", age_start := 85]
  setorder(city_meta, age_start)

  pop_fit <- ungroup::pclm(x = city_meta$age_start, y = city_meta$pop, nlast = 15)$fitted
  death_fit <- ungroup::pclm(x = city_meta$age_start, y = city_meta$death, nlast = 15)$fitted

  if (length(pop_fit) != 35 || length(death_fit) != 35) {
    fail("pclm test returned an unexpected fitted-length output.")
  }
  if (anyNA(pop_fit) || anyNA(death_fit)) {
    fail("pclm test returned NA values for the sample city.")
  }

  pass(sprintf("pclm disaggregation succeeded for %s", city_id))
}

run_test <- function(name, fn) {
  cat(sprintf("\n== Running %s ==\n", name))
  fn()
}

available_tests <- list(
  script03 = function() test_script03_output(test_cities, input_dir),
  script04 = function() test_script04_aggregation(test_cities, input_dir),
  script05 = function() test_script05_tables(test_cities, input_dir),
  pclm = function() test_pclm(test_cities)
)

if (component == "all") {
  for (name in names(available_tests)) {
    run_test(name, available_tests[[name]])
  }
} else if (component %in% names(available_tests)) {
  run_test(component, available_tests[[component]])
} else {
  fail(sprintf("Unknown component '%s'. Use one of: %s", component, paste(c("all", names(available_tests)), collapse = ", ")))
}

cat("\nAll requested checks completed.\n")