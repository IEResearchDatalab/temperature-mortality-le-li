#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: Rscript analysis/pooled_gcm_le_li_spread.R <source_repo> <output_dir> [ssp]", call. = FALSE)
}

source_repo <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
ssp <- if (length(args) == 3L) as.integer(args[[3]]) else 3L
if (!ssp %in% 1:3) stop("ssp must be 1, 2, or 3", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

years <- 2020:2099
ages <- 65:100
agegroups <- c("65-74", "75-84", "85+")
age_slices <- list("65-74" = 65:74, "75-84" = 75:84, "85+" = 85:100)
branches <- c("with_cc", "without_cc")
ranges <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
gcms <- c(
  "ACCESS_CM2", "ACCESS_ESM1_5", "BCC_CSM2_MR", "CMCC_ESM2", "CanESM5",
  "EC_Earth3", "EC_Earth3_Veg_LR", "GFDL_ESM4", "IITM_ESM", "INM_CM4_8",
  "INM_CM5_0", "IPSL_CM6A_LR", "KACE_1_0_G", "MIROC6", "MPI_ESM1_2_HR",
  "MPI_ESM1_2_LR", "MRI_ESM2_0", "NorESM2_LM", "NorESM2_MM"
)
regions <- c("Eastern", "Northern", "Southern", "Western")
geographies <- c("Europe", paste(regions, "Europe"))

ssp_root <- file.path(source_repo, "results", "europe", sprintf("ssp%d", ssp))
meta <- unique(fread(
  file.path(source_repo, "data", "city_results.csv"),
  select = c("URAU_CODE", "LABEL", "region")
)[, .(city = URAU_CODE, city_name = LABEL, region)])
if (anyDuplicated(meta$city)) stop("City metadata is not unique", call. = FALSE)
if (!setequal(unique(meta$region), regions)) stop("Unexpected region domain", call. = FALSE)
setorder(meta, city)
cities <- meta$city

expected <- CJ(city = cities, gcm = gcms)
expected[, file := file.path(ssp_root, city, gcm, "01_attribution_grouped.csv")]
expected[, `:=`(size = file.info(file)$size, mtime = as.character(file.info(file)$mtime))]
dem_inventory <- rbind(
  data.table(city = cities, gcm = "DEMOGRAPHY_SINGLE", file = file.path(ssp_root, cities, "demography", "00_demography_single_age.csv")),
  data.table(city = cities, gcm = "DEMOGRAPHY_GROUPED", file = file.path(ssp_root, cities, "demography", "00_demography_grouped.csv"))
)
dem_inventory[, `:=`(size = file.info(file)$size, mtime = as.character(file.info(file)$mtime))]
input_inventory <- rbind(expected, dem_inventory, fill = TRUE)
input_inventory[, exists := file.exists(file)]
fwrite(input_inventory, file.path(output_dir, "input_inventory.csv"))
if (any(!input_inventory$exists)) {
  fwrite(input_inventory[exists == FALSE], file.path(output_dir, "missing_inputs.csv"))
  stop("Required inputs are missing; see missing_inputs.csv", call. = FALSE)
}

life_expectancy_from_mx_65plus <- function(mx, x = ages, nx = c(rep(1, 35), Inf), age = 0) {
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  ex[age + 1]
}

sd_from_mx_fun_65_plus <- function(mx, x = ages, nx = c(rep(1, 35), Inf), age = 0) {
  x_conditional <- x - 65
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  sqrt(sum(dx * (x_conditional + 0.5 - ex[age + 1])^2))
}

an_dims <- c(length(gcms), length(geographies), length(branches), length(years), length(ages), length(ranges))
an_dimnames <- list(gcm = gcms, geography = geographies, branch = branches,
  year = as.character(years), age = as.character(ages), range = ranges)
dem_dims <- c(length(geographies), length(years), length(ages))
dem_dimnames <- list(geography = geographies, year = as.character(years), age = as.character(ages))

worker <- function(city_subset) {
  an_sum <- array(0, dim = an_dims, dimnames = an_dimnames)
  pop_sum <- array(0, dim = dem_dims, dimnames = dem_dimnames)
  death_sum <- array(0, dim = dem_dims, dimnames = dem_dimnames)
  stats <- list(files = 0L, grouped_rows = 0L, max_weight_error = 0,
    max_provided_weight_diff = 0, max_grouped_dem_diff = 0, max_an_reconstruction_diff = 0)

  for (city_code in city_subset) {
    city_region <- meta[city == city_code]$region
    geo_idx <- c(1L, match(paste(city_region, "Europe"), geographies))

    dem <- fread(
      file.path(ssp_root, city_code, "demography", "00_demography_single_age.csv"),
      select = c("year", "agegroup", "age", "pop", "death", "death_weight")
    )[year %in% years & age %in% ages]
    if (nrow(dem) != length(years) * length(ages) ||
        anyDuplicated(dem, by = c("year", "age")) ||
        !setequal(dem$year, years) || !setequal(dem$age, ages)) {
      stop(sprintf("Incomplete single-age demography for %s", city_code), call. = FALSE)
    }
    setorder(dem, year, age)
    pop_mat <- matrix(dem$pop, nrow = length(years), ncol = length(ages), byrow = TRUE)
    death_mat <- matrix(dem$death, nrow = length(years), ncol = length(ages), byrow = TRUE)
    for (gi in geo_idx) {
      pop_sum[gi, , ] <- pop_sum[gi, , ] + pop_mat
      death_sum[gi, , ] <- death_sum[gi, , ] + death_mat
    }

    dem[, group_death := sum(death), by = .(year, agegroup)]
    dem[, computed_weight := fifelse(group_death > 0, death / group_death, 0)]
    weight_check <- dem[, .(group_death = first(group_death),
      weight_sum = sum(computed_weight)), by = .(year, agegroup)]
    stats$max_weight_error <- max(stats$max_weight_error,
      max(abs(weight_check[group_death > 0]$weight_sum - 1)))
    stats$max_provided_weight_diff <- max(stats$max_provided_weight_diff,
      max(abs(dem$computed_weight - dem$death_weight)))

    grouped_dem <- fread(
      file.path(ssp_root, city_code, "demography", "00_demography_grouped.csv"),
      select = c("year", "agegroup", "death")
    )[year %in% years & agegroup %in% agegroups]
    calc_grouped <- dem[, .(calculated_death = sum(death)), by = .(year, agegroup)]
    grouped_compare <- merge(grouped_dem, calc_grouped, by = c("year", "agegroup"))
    if (nrow(grouped_compare) != length(years) * length(agegroups)) {
      stop(sprintf("Incomplete grouped demography for %s", city_code), call. = FALSE)
    }
    stats$max_grouped_dem_diff <- max(stats$max_grouped_dem_diff,
      max(abs(grouped_compare$death - grouped_compare$calculated_death)))

    weights <- lapply(agegroups, function(grp) {
      x <- dem[agegroup == grp][order(year, age)]
      expected_ages <- age_slices[[grp]]
      if (nrow(x) != length(years) * length(expected_ages) || !setequal(x$age, expected_ages)) {
        stop(sprintf("Bad age-band mapping for %s/%s", city_code, grp), call. = FALSE)
      }
      matrix(x$computed_weight, nrow = length(years), ncol = length(expected_ages), byrow = TRUE)
    })
    names(weights) <- agegroups

    for (gcm_i in seq_along(gcms)) {
      gcm <- gcms[[gcm_i]]
      ga <- fread(
        file.path(ssp_root, city_code, gcm, "01_attribution_grouped.csv"),
        select = c("year", "branch", "agegroup", "range", "an")
      )[year %in% years & branch %in% branches & agegroup %in% agegroups & range %in% ranges]
      expected_rows <- length(years) * length(branches) * length(agegroups) * length(ranges)
      if (nrow(ga) != expected_rows || anyDuplicated(ga, by = c("year", "branch", "agegroup", "range")) ||
          !setequal(ga$year, years) || !setequal(ga$branch, branches) ||
          !setequal(ga$agegroup, agegroups) || !setequal(ga$range, ranges) || any(!is.finite(ga$an))) {
        stop(sprintf("Incomplete grouped AN domain for %s/%s", city_code, gcm), call. = FALSE)
      }

      for (branch_i in seq_along(branches)) for (range_i in seq_along(ranges)) for (grp in agegroups) {
        vals_dt <- ga[branch == branches[[branch_i]] & range == ranges[[range_i]] & agegroup == grp][order(year)]
        if (!identical(vals_dt$year, years)) stop("Year ordering failed", call. = FALSE)
        vals <- vals_dt$an
        weighted <- weights[[grp]] * vals
        recon_error <- max(abs(rowSums(weighted) - vals))
        stats$max_an_reconstruction_diff <- max(stats$max_an_reconstruction_diff, recon_error)
        age_idx <- match(age_slices[[grp]], ages)
        for (gi in geo_idx) {
          an_sum[gcm_i, gi, branch_i, , age_idx, range_i] <-
            an_sum[gcm_i, gi, branch_i, , age_idx, range_i] + weighted
        }
      }
      stats$files <- stats$files + 1L
      stats$grouped_rows <- stats$grouped_rows + nrow(ga)
    }
  }
  list(an = an_sum, pop = pop_sum, death = death_sum, stats = stats)
}

ncores <- as.integer(Sys.getenv("NCORES", "8"))
ncores <- max(1L, min(ncores, length(cities)))
chunks <- split(cities, cut(seq_along(cities), breaks = ncores, labels = FALSE))
message(sprintf("Reading %d city-GCM files using %d workers", nrow(expected), ncores))
parts <- parallel::mclapply(chunks, worker, mc.cores = ncores)
if (any(vapply(parts, inherits, logical(1), "try-error"))) stop("At least one worker failed", call. = FALSE)

an_sum <- Reduce(`+`, lapply(parts, `[[`, "an"))
pop_sum <- Reduce(`+`, lapply(parts, `[[`, "pop"))
death_sum <- Reduce(`+`, lapply(parts, `[[`, "death"))
stat_names <- names(parts[[1]]$stats)
stats <- lapply(stat_names, function(nm) {
  values <- vapply(parts, function(x) x$stats[[nm]], numeric(1))
  if (grepl("^max_", nm)) max(values) else sum(values)
})
names(stats) <- stat_names

make_levels <- function(an_values) {
  out <- vector("list", length(gcms) * length(geographies) * length(branches))
  k <- 0L
  min_rest <- Inf
  for (gcm_i in seq_along(gcms)) for (geo_i in seq_along(geographies)) {
    without_temp <- apply(an_values[gcm_i, geo_i, "without_cc", , , , drop = FALSE], c(4, 5), sum)
    rest <- death_sum[geo_i, , ] - without_temp
    min_rest <- min(min_rest, rest)
    for (branch_i in seq_along(branches)) {
      temp <- apply(an_values[gcm_i, geo_i, branch_i, , , , drop = FALSE], c(4, 5), sum)
      total_deaths <- rest + temp
      mx <- total_deaths / pop_sum[geo_i, , ]
      k <- k + 1L
      out[[k]] <- data.table(
        gcm = gcms[[gcm_i]], geography = geographies[[geo_i]], branch = branches[[branch_i]], year = years,
        LE65 = apply(mx, 1, life_expectancy_from_mx_65plus),
        LI65 = apply(mx, 1, sd_from_mx_fun_65_plus)
      )
    }
  }
  list(levels = rbindlist(out), min_rest = min_rest)
}

annual_result <- make_levels(an_sum)
annual_levels <- annual_result$levels

make_period_levels <- function(annual_levels_source, an_values) {
  rows <- list(); k <- 0L
  period_starts <- seq(2020, 2095, by = 5)
  min_period_rate <- Inf
  for (gcm_i in seq_along(gcms)) for (geo_i in seq_along(geographies)) {
    without_temp <- apply(an_values[gcm_i, geo_i, "without_cc", , , , drop = FALSE], c(4, 5), sum)
    rest <- death_sum[geo_i, , ] - without_temp
    for (branch_i in seq_along(branches)) {
      temp <- apply(an_values[gcm_i, geo_i, branch_i, , , , drop = FALSE], c(4, 5), sum)
      annual_mx <- (rest + temp) / pop_sum[geo_i, , ]
      for (p in period_starts) {
        idx <- match(p:(p + 4L), years)
        period_mx <- colMeans(annual_mx[idx, , drop = FALSE])
        min_period_rate <- min(min_period_rate, period_mx)
        k <- k + 1L
        rows[[k]] <- data.table(
          gcm = gcms[[gcm_i]], geography = geographies[[geo_i]], branch = branches[[branch_i]],
          period_start = p, period = sprintf("%d-%d", p, p + 4L),
          LE65 = life_expectancy_from_mx_65plus(period_mx),
          LI65 = sd_from_mx_fun_65_plus(period_mx)
        )
      }
    }
  }
  list(levels = rbindlist(rows), min_rate = min_period_rate)
}

period_result <- make_period_levels(annual_levels, an_sum)
period_levels <- period_result$levels

climate_effect <- function(levels, time_cols) {
  wide <- dcast(levels, as.formula(paste(paste(c("gcm", "geography", time_cols), collapse = " + "), "~ branch")),
    value.var = c("LE65", "LI65"))
  wide[, `:=`(
    delta_LE65 = LE65_with_cc - LE65_without_cc,
    delta_LI65 = LI65_with_cc - LI65_without_cc
  )]
  wide
}

annual_effect <- climate_effect(annual_levels, "year")
period_effect <- climate_effect(period_levels, c("period_start", "period"))
period_change <- copy(period_effect)
period_change[, `:=`(
  delta_LE65_change_from_2020_2024 = delta_LE65 - delta_LE65[period_start == 2020L],
  delta_LI65_change_from_2020_2024 = delta_LI65 - delta_LI65[period_start == 2020L]
), by = .(gcm, geography)]

spread_one <- function(dt, group_cols, metric) {
  dt[, .(
    n_gcms = .N,
    mean = mean(get(metric)), median = median(get(metric)),
    p02_5 = as.numeric(quantile(get(metric), 0.025, names = FALSE)),
    p05 = as.numeric(quantile(get(metric), 0.05, names = FALSE)),
    p95 = as.numeric(quantile(get(metric), 0.95, names = FALSE)),
    p97_5 = as.numeric(quantile(get(metric), 0.975, names = FALSE)),
    minimum = min(get(metric)), maximum = max(get(metric)),
    n_positive = sum(get(metric) > 0), n_negative = sum(get(metric) < 0), n_zero = sum(get(metric) == 0)
  ), by = group_cols][, metric := metric][]
}

annual_spread <- rbindlist(lapply(c("delta_LE65", "delta_LI65"),
  function(m) spread_one(annual_effect, c("geography", "year"), m)))
period_spread <- rbindlist(lapply(c("delta_LE65", "delta_LI65"),
  function(m) spread_one(period_effect, c("geography", "period_start", "period"), m)))
period_change_spread <- rbindlist(lapply(c("delta_LE65_change_from_2020_2024", "delta_LI65_change_from_2020_2024"),
  function(m) spread_one(period_change, c("geography", "period_start", "period"), m)))

quantile_tables <- list(annual_spread, period_spread, period_change_spread)
quantiles_sensible <- all(vapply(quantile_tables, function(x) all(
  x$minimum <= x$p02_5 & x$p02_5 <= x$p05 & x$p05 <= x$median &
    x$median <= x$p95 & x$p95 <= x$p97_5 & x$p97_5 <= x$maximum & x$n_gcms == length(gcms)
), logical(1)))

# Reconstruct the production ENSEMBLE by averaging pooled cause counts before rates/life tables.
ensemble_an <- apply(an_sum, c(2, 3, 4, 5, 6), mean)
ensemble_levels_rows <- list(); ensemble_comparisons <- list(); count_comparisons <- list()
for (geo_i in seq_along(geographies)) {
  geography <- geographies[[geo_i]]
  ref_dir <- if (geography == "Europe") {
    file.path(source_repo, "results", "europe", "geographies", sprintf("ssp%d", ssp), "europe", "EUROPE")
  } else {
    file.path(source_repo, "results", "europe", "geographies", sprintf("ssp%d", ssp), "region", sub(" Europe$", "", geography))
  }
  without_temp <- apply(ensemble_an[geo_i, "without_cc", , , , drop = FALSE], c(3, 4), sum)
  rest <- death_sum[geo_i, , ] - without_temp
  calc_levels <- list()
  for (branch_i in seq_along(branches)) {
    temp <- apply(ensemble_an[geo_i, branch_i, , , , drop = FALSE], c(3, 4), sum)
    mx <- (rest + temp) / pop_sum[geo_i, , ]
    calc_levels[[branch_i]] <- data.table(
      geography = geography, branch = branches[[branch_i]], year = years,
      LE65_calculated = apply(mx, 1, life_expectancy_from_mx_65plus),
      LI65_calculated = apply(mx, 1, sd_from_mx_fun_65_plus)
    )
  }
  calc_levels <- rbindlist(calc_levels)
  ref_levels <- fread(file.path(ref_dir, "04_le_li_levels.csv"))[, .(
    branch, year, LE65_reference = LE65, LI65_reference = LI65
  )]
  cmp <- merge(calc_levels, ref_levels, by = c("branch", "year"))
  cmp[, `:=`(
    LE65_abs_diff = abs(LE65_calculated - LE65_reference),
    LI65_abs_diff = abs(LI65_calculated - LI65_reference)
  )]
  ensemble_comparisons[[geo_i]] <- cmp

  ref_master <- fread(file.path(ref_dir, "03_master_table.csv"))
  ref_temp <- ref_master[range %in% ranges, .(branch, year, age, range, an_reference = an)]
  calc_temp <- as.data.table(as.table(ensemble_an[geo_i, , , , , drop = FALSE]))
  setnames(calc_temp, c("geography", "branch", "year", "age", "range", "N"))
  calc_temp[, `:=`(geography = NULL, year = as.integer(as.character(year)), age = as.integer(as.character(age)))]
  setnames(calc_temp, "N", "an_calculated")
  count_cmp <- merge(calc_temp, ref_temp, by = c("branch", "year", "age", "range"))
  count_cmp[, abs_diff := abs(an_calculated - an_reference)]
  count_cmp[, geography := geography]
  count_comparisons[[geo_i]] <- count_cmp
}
ensemble_comparison <- rbindlist(ensemble_comparisons)
count_comparison <- rbindlist(count_comparisons)

all_numeric_finite <- all(vapply(list(annual_levels, annual_effect, period_levels, period_effect, period_change,
  annual_spread, period_spread, period_change_spread), function(x) {
    cols <- names(x)[vapply(x, is.numeric, logical(1))]
    all(vapply(x[, ..cols], function(v) all(is.finite(v)), logical(1)))
  }, logical(1)))

checks <- data.table(
  check = c(
    "input_file_coverage", "processed_city_gcm_files", "grouped_an_domain",
    "death_weights_sum_to_one", "provided_death_weights_match", "grouped_demography_reconstruction",
    "grouped_an_reconstruction", "positive_pooled_population", "nonnegative_rest_mortality",
    "annual_level_coverage", "period_level_coverage", "finite_outputs",
    "ensemble_count_reconstruction", "ensemble_LE65_reconstruction", "ensemble_LI65_reconstruction",
    "sensible_quantiles"
  ),
  value = c(
    sprintf("%d/%d", sum(input_inventory$exists), nrow(input_inventory)),
    sprintf("%d/%d", stats$files, nrow(expected)),
    sprintf("%d rows", stats$grouped_rows),
    sprintf("max_abs_error=%.3e", stats$max_weight_error),
    sprintf("max_abs_diff=%.3e", stats$max_provided_weight_diff),
    sprintf("max_abs_diff=%.3e", stats$max_grouped_dem_diff),
    sprintf("max_abs_diff=%.3e", stats$max_an_reconstruction_diff),
    sprintf("minimum=%.6f", min(pop_sum)),
    sprintf("minimum=%.6f", annual_result$min_rest),
    sprintf("%d rows", nrow(annual_levels)),
    sprintf("%d rows", nrow(period_levels)),
    as.character(all_numeric_finite),
    sprintf("max_abs_diff=%.3e", max(count_comparison$abs_diff)),
    sprintf("max_abs_diff=%.3e", max(ensemble_comparison$LE65_abs_diff)),
    sprintf("max_abs_diff=%.3e", max(ensemble_comparison$LI65_abs_diff)),
    as.character(quantiles_sensible)
  ),
  threshold = c(
    sprintf("%d files", nrow(input_inventory)), sprintf("%d files", nrow(expected)),
    sprintf("%d rows", nrow(expected) * length(years) * length(branches) * length(agegroups) * length(ranges)),
    "<= 1e-12", "<= 1e-12", "<= 1e-8", "<= 1e-9", "> 0", ">= -1e-9",
    sprintf("%d rows", length(gcms) * length(geographies) * length(branches) * length(years)),
    sprintf("%d rows", length(gcms) * length(geographies) * length(branches) * 16L),
    "TRUE", "<= 1e-6", "<= 1e-10", "<= 1e-10", "TRUE"
  ),
  status = c(
    if (all(input_inventory$exists)) "PASS" else "FAIL",
    if (stats$files == nrow(expected)) "PASS" else "FAIL",
    if (stats$grouped_rows == nrow(expected) * 1920L) "PASS" else "FAIL",
    if (stats$max_weight_error <= 1e-12) "PASS" else "FAIL",
    if (stats$max_provided_weight_diff <= 1e-12) "PASS" else "FAIL",
    if (stats$max_grouped_dem_diff <= 1e-8) "PASS" else "FAIL",
    if (stats$max_an_reconstruction_diff <= 1e-9) "PASS" else "FAIL",
    if (min(pop_sum) > 0) "PASS" else "FAIL",
    if (annual_result$min_rest >= -1e-9) "PASS" else "FAIL",
    if (nrow(annual_levels) == length(gcms) * length(geographies) * length(branches) * length(years)) "PASS" else "FAIL",
    if (nrow(period_levels) == length(gcms) * length(geographies) * length(branches) * 16L) "PASS" else "FAIL",
    if (all_numeric_finite) "PASS" else "FAIL",
    if (max(count_comparison$abs_diff) <= 1e-6) "PASS" else "FAIL",
    if (max(ensemble_comparison$LE65_abs_diff) <= 1e-10) "PASS" else "FAIL",
    if (max(ensemble_comparison$LI65_abs_diff) <= 1e-10) "PASS" else "FAIL",
    if (quantiles_sensible) "PASS" else "FAIL"
  )
)

fwrite(annual_levels, file.path(output_dir, "per_gcm_annual_levels.csv"))
fwrite(annual_effect, file.path(output_dir, "per_gcm_annual_climate_effect.csv"))
fwrite(annual_spread, file.path(output_dir, "annual_gcm_spread.csv"))
fwrite(period_levels, file.path(output_dir, "per_gcm_period_levels.csv"))
fwrite(period_effect, file.path(output_dir, "per_gcm_period_climate_effect.csv"))
fwrite(period_spread, file.path(output_dir, "period_gcm_spread.csv"))
fwrite(period_change, file.path(output_dir, "per_gcm_period_change_from_baseline.csv"))
fwrite(period_change_spread, file.path(output_dir, "period_change_gcm_spread.csv"))
fwrite(ensemble_comparison, file.path(output_dir, "ensemble_level_reconstruction.csv"))
fwrite(count_comparison[, .(geography, branch, year, age, range, an_calculated, an_reference, abs_diff)],
  file.path(output_dir, "ensemble_count_reconstruction.csv"))
fwrite(checks, file.path(output_dir, "validation_checks.csv"))

headline <- rbind(
  period_spread[period_start == 2095L],
  period_change_spread[period_start == 2095L],
  fill = TRUE
)
fwrite(headline, file.path(output_dir, "headline_findings.csv"))

geo_order <- c("Europe", "Southern Europe", "Eastern Europe", "Western Europe", "Northern Europe")
end_change <- period_change[period_start == 2095L]
end_change[, geography_plot := factor(geography, levels = rev(geo_order))]

p1 <- ggplot(end_change, aes(x = delta_LE65_change_from_2020_2024, y = geography_plot)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.5) +
  geom_boxplot(width = 0.52, outlier.shape = NA, fill = "#DCE6F1", colour = "#17365D") +
  geom_point(position = position_jitter(height = 0.10, width = 0), alpha = 0.65,
    size = 1.8, colour = "#17365D") +
  labs(
    title = "Climate-model spread in the erosion of LE65 gains",
    subtitle = "SSP3-7.0: 2095-2099 climate effect minus the 2020-2024 climate effect",
    x = "Change in LE65 climate effect (years)", y = NULL,
    caption = "Each point is one of 19 GCMs. Counts are pooled before life tables. Central ERF; not a confidence interval."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15), panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0, colour = "grey35"))
ggsave(file.path(output_dir, "01_gcm_spread_le65_change.png"), p1, width = 10, height = 5.8, dpi = 180)

eu_annual <- annual_spread[geography == "Europe" & metric == "delta_LE65"]
p2 <- ggplot(eu_annual, aes(x = year)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.5) +
  geom_ribbon(aes(ymin = p02_5, ymax = p97_5), fill = "#9ECAE1", alpha = 0.45) +
  geom_line(aes(y = median), colour = "#17365D", linewidth = 0.9) +
  geom_line(aes(y = mean), colour = "#B2182B", linewidth = 0.75, linetype = "dashed") +
  labs(
    title = "Europe's LE65 climate effect becomes increasingly negative",
    subtitle = "SSP3-7.0 across 19 GCMs; ribbon is empirical 2.5th-97.5th percentile spread",
    x = NULL, y = "LE65 with climate change minus without (years)",
    caption = "Solid: median GCM. Dashed: mean GCM. Counts pooled across 854 cities before life tables; central ERF."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15),
    plot.caption = element_text(hjust = 0, colour = "grey35"))
ggsave(file.path(output_dir, "02_europe_le65_effect_trajectory.png"), p2, width = 10, height = 5.8, dpi = 180)

p3 <- ggplot(end_change, aes(x = delta_LI65_change_from_2020_2024, y = geography_plot)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.5) +
  geom_boxplot(width = 0.52, outlier.shape = NA, fill = "#FCE4D6", colour = "#A13D2D") +
  geom_point(position = position_jitter(height = 0.10, width = 0), alpha = 0.65,
    size = 1.8, colour = "#A13D2D") +
  labs(
    title = "Climate-model spread in the LI65 change",
    subtitle = "SSP3-7.0: 2095-2099 climate effect minus the 2020-2024 climate effect",
    x = "Change in LI65 climate effect (SD years)", y = NULL,
    caption = "Each point is one of 19 GCMs. Central ERF; empirical climate-model spread, not a confidence interval."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15), panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0, colour = "grey35"))
ggsave(file.path(output_dir, "03_gcm_spread_li65_change.png"), p3, width = 10, height = 5.8, dpi = 180)

git_value <- function(...) {
  value <- tryCatch(system2("git", c("-C", normalizePath(".", mustWork = TRUE), ...), stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
  paste(value, collapse = " ")
}
tracked_status <- git_value("status", "--porcelain", "--untracked-files=no")
full_status <- git_value("status", "--porcelain")
invocation <- sprintf(
  "NCORES=%d Rscript analysis/pooled_gcm_le_li_spread.R %s %s %d",
  ncores, shQuote(source_repo), shQuote(output_dir), ssp
)
manifest <- c(
  sprintf("generated_at=%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  sprintf("server=%s", Sys.info()[["nodename"]]),
  sprintf("repository=%s", normalizePath(".", mustWork = TRUE)),
  sprintf("branch=%s", git_value("branch", "--show-current")),
  sprintf("commit=%s", git_value("rev-parse", "HEAD")),
  sprintf("tracked_git_status=%s", ifelse(nzchar(tracked_status), "dirty", "clean")),
  sprintf("full_git_status=%s", ifelse(nzchar(full_status), "dirty", "clean")),
  "script=analysis/pooled_gcm_le_li_spread.R",
  sprintf("source_repository=%s", source_repo),
  sprintf("input_pattern=%s", file.path(ssp_root, "<CITY>", "<GCM>", "01_attribution_grouped.csv")),
  sprintf("demography_pattern=%s", file.path(ssp_root, "<CITY>", "demography", "00_demography_{single_age,grouped}.csv")),
  sprintf("output_directory=%s", output_dir),
  sprintf("command=%s", invocation),
  sprintf("ssp=%d", ssp), sprintf("cities=%d", length(cities)), sprintf("gcms=%d", length(gcms)), sprintf("ncores=%d", ncores),
  sprintf("validation=%s", if (all(checks$status == "PASS")) "PASS" else "FAIL")
)
writeLines(manifest, file.path(output_dir, "run_manifest.txt"))

if (any(checks$status == "FAIL")) stop("Validation failed; see validation_checks.csv", call. = FALSE)
message("Completed pooled per-GCM LE65/LI65 spread analysis: ", output_dir)
