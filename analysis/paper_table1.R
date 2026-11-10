#!/usr/bin/env Rscript

################################################################################
# Paper Table 1 and supplement tables S1-S4.
# Presentation only: reads completed Part 00-04 outputs. It never writes into
# the batch root. Stages (TABLE1_STAGE): gcm_time, gcm, build, test, all.
#   gcm_time  time one per-GCM case (all 35 pooled geographies, SSP3)
#   gcm       pooled per-GCM decompositions for the 35 geographies, SSP1-3
#   build     Table 1, S1-S4 (CSV, HTML, PDF) in TABLE1_DIR
#   test      checks on the built table
# Run from the repository root: Rscript analysis/paper_table1.R
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(DemoDecomp)
  library(parallel)
})
source("analysis/decade_decomposition.R")

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
table_dir <- Sys.getenv("TABLE1_DIR", "results/figures_simon_20oct/table1")
scratch_dir <- Sys.getenv("TABLE1_SCRATCH", "results/figures_simon_20oct/_scratch_gcm")
n_horiuchi <- as.integer(Sys.getenv("N_HORIUCHI", "400"))
n_cores <- as.integer(Sys.getenv("TABLE1_CORES", "16"))
stage <- Sys.getenv("TABLE1_STAGE", "all")

temp_causes <- decade_cause_levels[1:5]
heat_causes <- c("ModHeat", "NormExtrHeat", "ExtrHeat")
cold_causes <- c("ExtrCold", "ModCold")
start_years <- 2020:2024
end_years <- 2095:2099
gcm_names <- c(
  "ACCESS_CM2", "ACCESS_ESM1_5", "BCC_CSM2_MR", "CanESM5", "CMCC_ESM2",
  "EC_Earth3", "EC_Earth3_Veg_LR", "GFDL_ESM4", "IITM_ESM", "INM_CM4_8",
  "INM_CM5_0", "IPSL_CM6A_LR", "KACE_1_0_G", "MIROC6", "MPI_ESM1_2_HR",
  "MPI_ESM1_2_LR", "MRI_ESM2_0", "NorESM2_LM", "NorESM2_MM"
)
scenario_labels <- c(`1` = "SSP1\u20132.6", `2` = "SSP2\u20134.5", `3` = "SSP3\u20137.0")
minus <- function(x) gsub("-", "\u2212", x, fixed = TRUE)

city_meta <- function() {
  unique(fread("data/city_results.csv", select = c("URAU_CODE", "CNTR_CODE", "cntr_name", "region")))
}

geo_members <- function(meta) {
  rbindlist(list(
    meta[, .(geo_id = URAU_CODE, level = "europe", id = "EUROPE")],
    meta[, .(geo_id = URAU_CODE, level = "region", id = region)],
    meta[, .(geo_id = URAU_CODE, level = "country", id = CNTR_CODE)]
  ))
}

#------------------------
# POOLED PER-GCM INPUTS (single-age ANs rebuilt as in Parts 02 and 03)
#------------------------

read_demography <- function(ssp_value, cities) {
  dem <- rbindlist(lapply(cities, function(city) {
    f <- file.path(batch_root, paste0("ssp", ssp_value), city, "demography", "00_demography_single_age.csv")
    x <- fread(f, select = c("geo_id", "ssp", "year", "agegroup", "age", "pop", "death"))
    x[ssp == ssp_value & year %in% c(start_years, end_years) & age %in% decade_age_levels]
  }))
  dem[, group_death := sum(death), by = .(geo_id, year, agegroup)]
  dem[, weight := fifelse(group_death > 0, death / group_death, 0)]
  dem
}

# Pooled deaths and population by geography, branch, year, age and cause.
# Part 02: single-age AN = grouped AN x within-group death weight.
# Part 03: rest = all-cause deaths minus without-CC temperature deaths.
pooled_case <- function(ssp_value, gcm_name, dem, members) {
  an <- rbindlist(lapply(unique(dem$geo_id), function(city) {
    f <- file.path(batch_root, paste0("ssp", ssp_value), city, gcm_name, "01_attribution_grouped.csv")
    x <- fread(f, select = c("year", "branch", "agegroup", "range", "an"))
    x[year %in% c(start_years, end_years)][, geo_id := city]
  }))
  single <- merge(an, dem[, .(geo_id, year, agegroup, age, weight)],
                  by = c("geo_id", "year", "agegroup"), allow.cartesian = TRUE)
  single[, deaths := an * weight]
  without_temp <- single[branch == "without_cc", .(without_temp = sum(deaths)), by = .(geo_id, year, age)]
  rest <- merge(dem[, .(geo_id, year, age, death)], without_temp, by = c("geo_id", "year", "age"))
  rest[, deaths := death - without_temp]
  rest <- rbindlist(lapply(c("with_cc", "without_cc"), function(b) {
    rest[, .(geo_id, year, age, branch = b, range = "rest", deaths)]
  }))
  city_deaths <- rbind(single[, .(geo_id, year, age, branch, range, deaths)], rest)
  deaths <- merge(city_deaths, members, by = "geo_id", allow.cartesian = TRUE)[
    , .(deaths = sum(deaths)), by = .(level, id, branch, year, age, cause = range)]
  pop <- merge(dem[, .(geo_id, year, age, pop)], members, by = "geo_id", allow.cartesian = TRUE)[
    , .(pop = sum(pop)), by = .(level, id, year, age)]
  merge(deaths, pop, by = c("level", "id", "year", "age"))
}

period_names <- c(start = "2020-2024", end = "2095-2099")
horiuchi_steps <- c(400L, 1600L, 6400L)

# Gradient follows the cause-major parameter order used by decade_cause_vector().
life_expectancy_gradient <- function(mx) {
  mx <- matrix(mx, nrow = length(decade_age_levels), ncol = length(decade_cause_levels))
  total_mx <- rowSums(mx)
  px <- exp(-total_mx * decade_nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(total_mx == 0, lx * decade_nx, dx / total_mx)
  future_L <- rev(cumsum(rev(Lx))) - Lx

  derivative_interval <- numeric(length(total_mx))
  finite <- is.finite(decade_nx)
  finite_index <- which(finite)
  rates <- total_mx[finite_index]
  intervals <- decade_nx[finite_index]
  small <- abs(rates * intervals) < 1e-4
  small_index <- finite_index[small]
  small_rates <- rates[small]
  small_intervals <- intervals[small]
  derivative_interval[small_index] <-
    -small_intervals^2 / 2 +
    small_intervals^3 * small_rates / 3 -
    small_intervals^4 * small_rates^2 / 8 +
    small_intervals^5 * small_rates^3 / 30 -
    small_intervals^6 * small_rates^4 / 144
  regular <- !small
  regular_index <- finite_index[regular]
  regular_rates <- rates[regular]
  regular_intervals <- intervals[regular]
  derivative_interval[regular_index] <-
    (regular_intervals * regular_rates *
       exp(-regular_intervals * regular_rates) +
       expm1(-regular_intervals * regular_rates)) /
    regular_rates^2
  derivative_interval[!finite] <- -1 / total_mx[!finite]^2

  gradient_by_age <- lx * derivative_interval
  gradient_by_age[finite] <- gradient_by_age[finite] -
    decade_nx[finite] * future_L[finite]
  rep(gradient_by_age, times = length(decade_cause_levels))
}

gauss_legendre_01 <- function(n) {
  off_diagonal <- seq_len(n - 1L) / sqrt(4 * seq_len(n - 1L)^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- off_diagonal
  jacobi[cbind(2:n, seq_len(n - 1L))] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  order_nodes <- order(decomposition$values)
  list(
    nodes = (decomposition$values[order_nodes] + 1) / 2,
    weights = decomposition$vectors[1, order_nodes]^2
  )
}

horiuchi_analytic_gradient <- function(pars1, pars2, n_nodes) {
  # Integrate the straight-path cause derivatives without differencing LE values.
  quadrature <- gauss_legendre_01(n_nodes)
  change <- pars2 - pars1
  mean_gradient <- numeric(length(change))
  for (i in seq_along(quadrature$nodes)) {
    mean_gradient <- mean_gradient +
      quadrature$weights[i] *
      life_expectancy_gradient(pars1 + quadrature$nodes[i] * change)
  }
  change * mean_gradient
}

period_mean_schedule <- function(x) {
  x <- copy(x)
  x[, period := fifelse(
    year %in% start_years, period_names[["start"]],
    fifelse(year %in% end_years, period_names[["end"]], NA_character_)
  )]
  x <- x[!is.na(period)]
  keys <- c("level", "id", "branch", "year", "age", "cause")
  if (anyDuplicated(x, by = keys)) {
    stop("Pooled annual counts contain duplicate cause-specific rows.", call. = FALSE)
  }

  annual_rates <- x[, .(
    mx_cause = sum(deaths) / unique(pop),
    n_population_values = uniqueN(pop)
  ), by = c("level", "id", "branch", "period", "year", "age", "cause")]
  if (any(annual_rates$n_population_values != 1L) ||
      any(!is.finite(annual_rates$mx_cause)) ||
      any(annual_rates$mx_cause < 0)) {
    stop("Pooled annual rates have invalid populations or mortality rates.", call. = FALSE)
  }

  schedule <- annual_rates[, .(
    mx_cause = mean(mx_cause),
    years = uniqueN(year)
  ), by = c("level", "id", "branch", "period", "age", "cause")]
  expected_per_geo <- 2L * 2L * length(decade_age_levels) * length(decade_cause_levels)
  geo_counts <- schedule[, .N, by = .(level, id)]
  if (any(geo_counts$N != expected_per_geo) ||
      any(schedule$years != 5L) ||
      anyDuplicated(schedule, by = c("level", "id", "branch", "period", "age", "cause"))) {
    stop("Period-mean schedules do not form complete five-year grids.", call. = FALSE)
  }
  schedule[, years := NULL]
  schedule
}

decompose_between_branches <- function(schedule, period, context) {
  period_name <- period
  by_branch <- dcast(
    schedule[period == period_name & cause == "rest"],
    age ~ branch, value.var = "mx_cause"
  )
  rest_difference <- max(abs(by_branch$with_cc - by_branch$without_cc))
  if (rest_difference > 1e-12) {
    stop(context, " ", period, ": rest mortality differs between branches by ",
         format(rest_difference, scientific = TRUE), ".", call. = FALSE)
  }

  mx_without <- decade_cause_vector(schedule[period == period_name & branch == "without_cc"])
  mx_with <- decade_cause_vector(schedule[period == period_name & branch == "with_cc"])
  le_without <- decade_life_expectancy(mx_without)
  le_with <- decade_life_expectancy(mx_with)
  closure_400 <- NA_real_
  result <- NULL

  for (n in horiuchi_steps) {
    contribution <- as.vector(horiuchi(
      func = decade_life_expectancy,
      pars1 = mx_without,
      pars2 = mx_with,
      N = n
    ))
    closure_months <- 12 * (sum(contribution) - (le_with - le_without))
    if (n == horiuchi_steps[1]) closure_400 <- closure_months
    result <- list(
      values = decade_grid(contribution),
      closure_months = closure_months,
      n = n,
      le_without = le_without,
      le_with = le_with,
      method = "DemoDecomp",
      quadrature_nodes = 0L
    )
    if (abs(closure_months) <= 1e-8) break
  }

  numeric_last_closure <- result$closure_months
  if (abs(result$closure_months) > 1e-6) {
    for (n_nodes in c(8L, 16L, 32L, 64L)) {
      contribution <- horiuchi_analytic_gradient(mx_without, mx_with, n_nodes)
      closure_months <- 12 * (sum(contribution) - (le_with - le_without))
      result <- list(
        values = decade_grid(contribution),
        closure_months = closure_months,
        n = result$n,
        le_without = le_without,
        le_with = le_with,
        method = "analytic-gradient Gauss-Legendre",
        quadrature_nodes = n_nodes
      )
      if (abs(closure_months) <= 1e-10) break
    }
  }

  if (abs(result$closure_months) > 1e-6) {
    stop(
      context, " ", period, ": Horiuchi closure remains ",
      format(result$closure_months, scientific = TRUE),
      " months at N=", result$n, ".", call. = FALSE
    )
  }
  result$closure_400_months <- closure_400
  result$closure_numeric_last_months <- numeric_last_closure
  result
}

compute_period_metrics <- function(schedule, context) {
  start <- decompose_between_branches(schedule, period_names[["start"]], context)
  end <- decompose_between_branches(schedule, period_names[["end"]], context)
  g0 <- end$le_without - start$le_without
  g1 <- end$le_with - start$le_with
  loss_by_group <- function(parts, causes) {
    parts$values[cause %in% causes, sum(contribution)]
  }

  heat_loss <- 12 * (loss_by_group(start, heat_causes) - loss_by_group(end, heat_causes))
  cold_loss <- 12 * (loss_by_group(start, cold_causes) - loss_by_group(end, cold_causes))
  rest_loss <- 12 * (loss_by_group(start, "rest") - loss_by_group(end, "rest"))
  all_loss <- 12 * (g0 - g1)
  additivity_gap <- heat_loss + cold_loss + rest_loss - all_loss
  closure_start <- start$closure_months
  closure_end <- end$closure_months
  if (abs(rest_loss) >= 1e-9) {
    stop(context, ": rest loss is ", format(rest_loss, scientific = TRUE),
         " months although rest rates match.", call. = FALSE)
  }
  if (abs(additivity_gap) >= 1e-6) {
    stop(context, ": additivity gap is ", format(additivity_gap, scientific = TRUE),
         " months.", call. = FALSE)
  }

  data.table(
    g0_months = 12 * g0,
    g1_months = 12 * g1,
    all_loss_months = all_loss,
    heat_loss_months = heat_loss,
    cold_loss_months = cold_loss,
    rest_loss_months = rest_loss,
    all_ratio_pct = 100 * all_loss / (12 * g0),
    additivity_gap_months = additivity_gap,
    closure_start_n400_months = start$closure_400_months,
    closure_end_n400_months = end$closure_400_months,
    closure_start_numeric_last_months = start$closure_numeric_last_months,
    closure_end_numeric_last_months = end$closure_numeric_last_months,
    closure_start_final_months = closure_start,
    closure_end_final_months = closure_end,
    n_start = start$n,
    n_end = end$n,
    method_start = start$method,
    method_end = end$method,
    quadrature_nodes_start = start$quadrature_nodes,
    quadrature_nodes_end = end$quadrature_nodes
  )
}

run_case <- function(ssp_value, gcm_name, dem, members, out_file) {
  t0 <- proc.time()[["elapsed"]]
  pooled <- pooled_case(ssp_value, gcm_name, dem, members)
  schedule <- period_mean_schedule(pooled)
  res <- schedule[, compute_period_metrics(
    .SD, sprintf("SSP%d %s/%s %s", ssp_value, level, id, gcm_name)
  ), by = .(level, id)]
  res[, `:=`(ssp = ssp_value, gcm = gcm_name, seconds = proc.time()[["elapsed"]] - t0)]
  setcolorder(res, c("ssp", "gcm", "level", "id"))
  fwrite(res, out_file)
  res
}

run_gcm_stage <- function(timing_only = FALSE) {
  dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- city_meta()
  members <- geo_members(meta)
  cities <- unique(meta$URAU_CODE)
  ssps <- if (timing_only) 3L else 1:3
  dem <- setNames(lapply(ssps, function(s) {
    f <- file.path(scratch_dir, paste0("demography_ssp", s, ".rds"))
    if (file.exists(f)) return(readRDS(f))
    d <- read_demography(s, cities)
    saveRDS(d, f)
    d
  }), as.character(ssps))
  cases <- if (timing_only) data.table(ssp = 3L, gcm = gcm_names[8]) else CJ(ssp = ssps, gcm = gcm_names)
  work <- function(i) {
    s <- cases$ssp[i]
    g <- cases$gcm[i]
    f <- file.path(scratch_dir, sprintf("percase_periodmean_gradient_ssp%d_%s.csv", s, g))
    prior_cache <- file.path(
      scratch_dir, sprintf("percase_periodmean_adaptive_ssp%d_%s.csv", s, g)
    )
    if (!timing_only && file.exists(f)) return(fread(f))
    if (!timing_only && file.exists(prior_cache)) {
      cached <- fread(prior_cache)
      cache_valid <- nrow(cached) == 35L &&
        all(abs(cached$additivity_gap_months) < 1e-6) &&
        all(abs(cached$rest_loss_months) < 1e-9) &&
        all(pmax(abs(cached$closure_start_final_months),
                 abs(cached$closure_end_final_months)) < 1e-6)
      if (cache_valid) {
        cached[, `:=`(
          closure_start_numeric_last_months = closure_start_final_months,
          closure_end_numeric_last_months = closure_end_final_months,
          method_start = "DemoDecomp",
          method_end = "DemoDecomp",
          quadrature_nodes_start = 0L,
          quadrature_nodes_end = 0L
        )]
        fwrite(cached, f)
        return(cached)
      }
    }
    out_file <- if (timing_only) tempfile(fileext = ".csv") else f
    result <- run_case(s, g, dem[[as.character(s)]], members, out_file)
    if (timing_only) unlink(out_file)
    result
  }
  t0 <- proc.time()[["elapsed"]]
  case_results <- mclapply(
    seq_len(nrow(cases)), work,
    mc.cores = if (timing_only) 1L else n_cores,
    mc.preschedule = FALSE
  )
  failures <- vapply(case_results, inherits, logical(1), what = "try-error")
  if (any(failures)) {
    failed_cases <- cases[failures]
    details <- vapply(case_results[failures], as.character, character(1))
    stop("Per-GCM decomposition failed: ",
         paste(sprintf("SSP%s %s: %s", failed_cases$ssp, failed_cases$gcm, details),
               collapse = "; "), call. = FALSE)
  }
  out <- rbindlist(case_results, use.names = TRUE, fill = TRUE)
  elapsed <- proc.time()[["elapsed"]] - t0
  if (timing_only) {
    message(sprintf("One case (35 pooled geographies, SSP3, %s): %.1f s", cases$gcm, elapsed))
    message(sprintf("Projected total for 57 cases on %d cores: %.1f min",
                    n_cores, 57 * elapsed / 60 / n_cores))
    return(invisible(out))
  }
  expected_rows <- 3L * length(gcm_names) * 35L
  if (nrow(out) != expected_rows || anyDuplicated(out, by = c("ssp", "gcm", "level", "id"))) {
    stop("Per-GCM results do not cover all 35 geographies, 3 SSPs and 19 GCMs.",
         call. = FALSE)
  }
  if (any(abs(out$additivity_gap_months) >= 1e-6) ||
      any(abs(out$rest_loss_months) >= 1e-9) ||
      any(pmax(abs(out$closure_start_final_months),
               abs(out$closure_end_final_months)) > 1e-6) ||
      any(abs(out$all_ratio_pct -
                100 * (out$heat_loss_months + out$cold_loss_months) /
                  out$g0_months) >= 1e-6) ||
      any(abs(out$heat_loss_months + out$cold_loss_months +
                out$rest_loss_months - out$all_loss_months) >= 1e-6)) {
    stop("Per-GCM additivity, rest, or Horiuchi closure checks failed.",
         call. = FALSE)
  }
  fwrite(out, file.path(scratch_dir, "per_gcm_pooled_results.csv"))
  escalations <- rbindlist(list(
    out[n_start > 400L, .(
      ssp, gcm, level, id, period = period_names[["start"]],
      N = n_start, closure_n400_months = closure_start_n400_months,
      closure_numeric_last_months = closure_start_numeric_last_months,
      closure_final_months = closure_start_final_months,
      method = method_start, quadrature_nodes = quadrature_nodes_start
    )],
    out[n_end > 400L, .(
      ssp, gcm, level, id, period = period_names[["end"]],
      N = n_end, closure_n400_months = closure_end_n400_months,
      closure_numeric_last_months = closure_end_numeric_last_months,
      closure_final_months = closure_end_final_months,
      method = method_end, quadrature_nodes = quadrature_nodes_end
    )]
  ), use.names = TRUE)
  fwrite(escalations, file.path(scratch_dir, "horiuchi_escalations_per_gcm.csv"))
  message(sprintf(
    "Per-GCM stage finished in %.1f min (%d cases); %d decompositions needed N > 400.",
    elapsed / 60, nrow(cases), nrow(escalations)
  ))
  if (nrow(escalations)) {
    message("Escalated per-GCM decompositions are listed in ",
            file.path(scratch_dir, "horiuchi_escalations_per_gcm.csv"))
  }
  invisible(out)
}

#------------------------
# ENSEMBLE COLUMNS
#------------------------

geo_dir <- function(ssp_value, level, id) {
  file.path(batch_root, "geographies", paste0("ssp", ssp_value), level, id)
}

# Competing-risk lifetime risk from 65 (percent): all temperature causes and heat only.
lifetime_risk <- function(mx_by_cause) {
  m <- dcast(mx_by_cause, age ~ cause, value.var = "mx")[order(age)]
  m_tot <- rowSums(m[, ..decade_cause_levels])
  q <- 1 - exp(-m_tot * decade_nx)
  l <- head(cumprod(c(1, 1 - q)), -1)
  d_cause <- l * q * (m[, ..decade_cause_levels] / m_tot)
  c(all = 100 * sum(d_cause[, ..temp_causes]), heat = 100 * sum(d_cause[, ..heat_causes]))
}

first_negative_decade <- function(x) {
  hit <- x[value < 0, decade_to]
  if (length(hit)) paste0(min(hit), "s") else "not before 2100"
}

ensemble_row <- function(ssp_value, level, id, decade_with, decade_without) {
  d <- geo_dir(ssp_value, level, id)
  master <- fread(file.path(d, "03_master_table.csv"),
                  select = c("label", "branch", "year", "age", "range", "an", "deaths_component", "pop"))
  end <- master[year %in% end_years & age %in% decade_age_levels]
  with_deaths <- end[range %in% temp_causes & branch == "with_cc", .(d = sum(an)), by = range]
  with_total <- with_deaths[, sum(d)]
  without_total <- end[range %in% temp_causes & branch == "without_cc", sum(an)]
  life_schedule <- end[branch == "with_cc", .(year, age, cause = range, mx = deaths_component / pop)][
    , .(mx = mean(mx)), by = .(age, cause)]
  risk <- lifetime_risk(life_schedule)
  metric_input <- master[year %in% c(start_years, end_years), .(
    branch, year, age, cause = range, deaths = deaths_component, pop
  )]
  geo_level <- level
  geo_id <- id
  metric_input[, `:=`(level = geo_level, id = geo_id)]
  period_schedule <- period_mean_schedule(metric_input)
  metrics <- compute_period_metrics(
    period_schedule[, .(branch, period, age, cause, mx_cause)],
    sprintf("SSP%d %s/%s ENSEMBLE", ssp_value, geo_level, geo_id)
  )
  share <- function(causes) 100 * with_deaths[range %in% causes, sum(d)] / with_total
  by_decade <- function(x) x[ssp == ssp_value & cause %in% temp_causes,
                             .(value = sum(le_contribution_months)), by = decade_to]
  dw <- by_decade(decade_with)
  dn <- by_decade(decade_without)
  diff_dt <- merge(dw, dn, by = "decade_to", suffixes = c("", "_wo"))[, .(decade_to, value = value - value_wo)]
  data.table(
    ssp = ssp_value, level = level, id = id, geography = master$label[1],
    lifetime_risk_pct = risk[["all"]], lifetime_risk_heat_pct = risk[["heat"]],
    persistent_share_pct = 100 * without_total / with_total,
    general_pct = share(c("ModCold", "ModHeat")),
    sustained_pct = share("NormExtrHeat"),
    episodic_pct = share(c("ExtrCold", "ExtrHeat")),
    g0_months = metrics$g0_months,
    g1_months = metrics$g1_months,
    all_loss_months = metrics$all_loss_months,
    all_ratio_pct = metrics$all_ratio_pct,
    heat_loss_ensemble_months = metrics$heat_loss_months,
    cold_loss_ensemble_months = metrics$cold_loss_months,
    rest_loss_ensemble_months = metrics$rest_loss_months,
    additivity_gap_months = metrics$additivity_gap_months,
    closure_start_n400_months = metrics$closure_start_n400_months,
    closure_end_n400_months = metrics$closure_end_n400_months,
    closure_start_numeric_last_months = metrics$closure_start_numeric_last_months,
    closure_end_numeric_last_months = metrics$closure_end_numeric_last_months,
    closure_start_final_months = metrics$closure_start_final_months,
    closure_end_final_months = metrics$closure_end_final_months,
    n_start = metrics$n_start,
    n_end = metrics$n_end,
    method_start = metrics$method_start,
    method_end = metrics$method_end,
    quadrature_nodes_start = metrics$quadrature_nodes_start,
    quadrature_nodes_end = metrics$quadrature_nodes_end,
    net_drag_decade = first_negative_decade(dw[order(decade_to)]),
    net_drag_alt_decade = first_negative_decade(diff_dt[order(decade_to)])
  )
}

decade_tables <- function(geos) {
  dir.create(file.path(scratch_dir, "decade"), recursive = TRUE, showWarnings = FALSE)
  jobs <- CJ(i = seq_len(nrow(geos)), branch = c("with_cc", "without_cc"))
  invisible(mclapply(seq_len(nrow(jobs)), function(j) {
    g <- geos[jobs$i[j]]
    suffix <- if (jobs$branch[j] == "with_cc") "" else "_without_cc"
    f <- file.path(scratch_dir, "decade", sprintf("decade_decomposition_%s_%s%s.csv", g$level, g$id, suffix))
    if (!file.exists(f)) {
      compute_decade_decomposition(batch_root, g$level, g$id, 1:3, file.path(scratch_dir, "decade"),
                                   n_horiuchi, branch = jobs$branch[j])
    }
  }, mc.cores = n_cores))
  NULL
}

read_decade <- function(level, id, branch) {
  suffix <- if (branch == "with_cc") "" else "_without_cc"
  fread(file.path(scratch_dir, "decade", sprintf("decade_decomposition_%s_%s%s.csv", level, id, suffix)))
}

build_geo_index <- function(meta) {
  region_order <- c("Northern", "Western", "Eastern", "Southern")
  countries <- unique(meta[, .(id = CNTR_CODE, geography = cntr_name)])[order(geography)]
  rbindlist(list(
    data.table(level = "europe", id = "EUROPE", geography = "Europe"),
    data.table(level = "region", id = region_order, geography = paste(region_order, "Europe")),
    data.table(level = "country", id = countries$id, geography = countries$geography)
  ))
}

collect_table <- function() {
  geos <- build_geo_index(city_meta())
  decade_tables(geos)
  decades <- setNames(lapply(seq_len(nrow(geos)), function(i) {
    list(with = read_decade(geos$level[i], geos$id[i], "with_cc"),
         without = read_decade(geos$level[i], geos$id[i], "without_cc"))
  }), paste(geos$level, geos$id))
  rows <- rbindlist(lapply(1:3, function(s) {
    rbindlist(lapply(seq_len(nrow(geos)), function(i) {
      dd <- decades[[paste(geos$level[i], geos$id[i])]]
      r <- ensemble_row(s, geos$level[i], geos$id[i], dd$with, dd$without)
      r[, geography := geos$geography[i]]
      r
    }))
  }))
  per <- fread(file.path(scratch_dir, "per_gcm_pooled_results.csv"))
  per[, `:=`(
    heat_ratio_pct = 100 * heat_loss_months / g0_months,
    all_ratio_gcm_pct = 100 * all_loss_months / g0_months
  )]
  setorder(per, ssp, level, id, gcm)
  wide <- dcast(per, ssp + level + id ~ gcm, value.var = "heat_ratio_pct")
  setnames(wide, gcm_names, paste0("heat_ratio_", gcm_names))
  summ <- per[, .(
    heat_ratio_median_pct = median(heat_ratio_pct),
    heat_ratio_min_pct = min(heat_ratio_pct),
    heat_ratio_max_pct = max(heat_ratio_pct),
    all_ratio_gcm_median_pct = median(all_ratio_gcm_pct),
    max_additivity_gap_months = max(abs(additivity_gap_months)),
    max_rest_loss_months = max(abs(rest_loss_months)),
    max_closure_error_months = max(abs(c(
      closure_start_final_months, closure_end_final_months
    ))),
    max_closure_numeric_last_months = max(abs(c(
      closure_start_numeric_last_months, closure_end_numeric_last_months
    ))),
    analytic_gradient_fallbacks = sum(
      method_start == "analytic-gradient Gauss-Legendre" |
        method_end == "analytic-gradient Gauss-Legendre"
    ),
    max_n = max(n_start, n_end)
  ), by = .(ssp, level, id)]
  out <- Reduce(function(a, b) merge(a, b, by = c("ssp", "level", "id")), list(rows, summ, wide))
  if (any(abs(rows$additivity_gap_months) >= 1e-6) ||
      any(abs(rows$rest_loss_ensemble_months) >= 1e-9) ||
      any(abs(c(rows$closure_start_final_months, rows$closure_end_final_months)) > 1e-6)) {
    stop("Ensemble period-mean decomposition checks failed.", call. = FALSE)
  }

  per_escalations <- rbindlist(list(
    per[n_start > 400L, .(
      ssp, gcm, level, id, period = period_names[["start"]],
      N = n_start, closure_n400_months = closure_start_n400_months,
      closure_numeric_last_months = closure_start_numeric_last_months,
      closure_final_months = closure_start_final_months,
      method = method_start, quadrature_nodes = quadrature_nodes_start
    )],
    per[n_end > 400L, .(
      ssp, gcm, level, id, period = period_names[["end"]],
      N = n_end, closure_n400_months = closure_end_n400_months,
      closure_numeric_last_months = closure_end_numeric_last_months,
      closure_final_months = closure_end_final_months,
      method = method_end, quadrature_nodes = quadrature_nodes_end
    )]
  ), use.names = TRUE)
  ensemble_escalations <- rbindlist(list(
    rows[n_start > 400L, .(
      ssp, gcm = "ENSEMBLE", level, id, period = period_names[["start"]],
      N = n_start, closure_n400_months = closure_start_n400_months,
      closure_numeric_last_months = closure_start_numeric_last_months,
      closure_final_months = closure_start_final_months,
      method = method_start, quadrature_nodes = quadrature_nodes_start
    )],
    rows[n_end > 400L, .(
      ssp, gcm = "ENSEMBLE", level, id, period = period_names[["end"]],
      N = n_end, closure_n400_months = closure_end_n400_months,
      closure_numeric_last_months = closure_end_numeric_last_months,
      closure_final_months = closure_end_final_months,
      method = method_end, quadrature_nodes = quadrature_nodes_end
    )]
  ), use.names = TRUE)
  escalations <- rbindlist(list(per_escalations, ensemble_escalations), use.names = TRUE)
  fwrite(escalations, file.path(table_dir, "table1_n_escalation.csv"))
  out[, hidden := g0_months < 1]
  out[, all_ratio_shown_pct := fifelse(hidden, NA_real_, all_ratio_pct)]
  out[, gcm_minus_ensemble_pp := all_ratio_gcm_median_pct - all_ratio_pct]
  out[, geo_order := match(paste(level, id), paste(geos$level, geos$id))]
  setorder(out, ssp, geo_order)
  out
}

#------------------------
# FORMATTED TABLES
#------------------------

f1 <- function(x) {
  out <- sprintf("%.1f", x)
  out[out == "-0.0"] <- "0.0"
  ifelse(is.na(x), "n/a", minus(out))
}

heat_cell <- function(tab) {
  ifelse(tab$hidden, "n/a", sprintf("%s (%s to %s)", f1(tab$heat_ratio_median_pct),
                                    f1(tab$heat_ratio_min_pct), f1(tab$heat_ratio_max_pct)))
}

format_table1 <- function(tab, name_col = "Region") {
  d <- data.table(
    Scenario = scenario_labels[as.character(tab$ssp)],
    geography = tab$geography,
    c1 = f1(tab$lifetime_risk_pct), c2 = f1(tab$persistent_share_pct),
    c3 = f1(tab$general_pct), c4 = f1(tab$sustained_pct),
    c5 = heat_cell(tab), c6 = f1(tab$all_ratio_shown_pct), c7 = tab$net_drag_decade
  )
  setnames(d, "geography", name_col)
  d
}

table_title <- paste0(
  "Table 1 Temperature-related mortality among people aged 65+ in 854 European cities, ",
  "by scenario and region:\nthe whole burden in 2095\u201399, and erosion of LE65 gains, ",
  "2020\u201324 to 2095\u201399"
)

footnotes <- c(
  "a Probability that a person reaching 65 eventually dies from a temperature-attributable cause, under the mortality conditions of 2095\u201399, with climate change.",
  "b Share of all temperature-attributable deaths at 65+ that would occur even without warming beyond the 2010\u201314 temperature distribution. One minus this share is the part attributable to climate change.",
  "c,d Shares of all temperature-attributable deaths at 65+, by the response each temperature range mainly calls for: general public health: moderate cold and moderate heat, sustained heat-focused measures: normal extreme heat, the remainder (episodic responses) is extreme heat and extreme cold.",
  "e,f Percentage of the gain in LE65 between 2020\u201324 and 2095\u201399 that would be expected from social progress alone (without climate change) but is lost because of climate change. \u201ce\u201d counts heat only; \u201cf\u201d counts all temperature ranges, so falling cold mortality offsets rising heat. Negative values mean climate change adds to the gain. Ratios are not shown where the gain without climate change is below one month. GCM median, with range in brackets for column \u201ce\u201d. Losses in months, and the gains, are in Table S1.",
  "g First decade in which temperature-related mortality, all ranges together, contributes negatively to the gain in LE65.",
  "Regions: UN M49 groupings, as in Masselot et al. (2025). Central estimates; country results in Tables S2\u2013S4."
)

table_notes <- c(
  "Note 1. SSP1\u20132.6: Wittgenstein Centre survival ratios stored as exactly 1 were set to 0.9995 (sensitivity).",
  "Note 2. SSP2\u20134.5 includes 12 Swiss cities with within-age-band mortality-shape warnings (maximum decline 11.09%)."
)

make_flextable <- function(display, header, sup, title, notes, merge_scenario = TRUE,
                           n_left = 2L, footer_after_fit = FALSE) {
  suppressPackageStartupMessages(library(flextable))
  ft <- flextable(display)
  for (j in seq_along(header)) {
    ft <- compose(ft, j = j, part = "header", value = as_paragraph(
      as_chunk(header[j]), if (nzchar(sup[j])) as_sup(sup[j]) else as_chunk("")))
  }
  ft <- ft |> theme_booktabs(bold_header = TRUE) |>
    fontsize(size = 8, part = "all") |> font(fontname = "Arial", part = "all") |>
    align(j = seq_len(ncol(display)), align = "center", part = "all") |>
    align(j = seq_len(n_left), align = "left", part = "all") |>
    valign(valign = "center", part = "all") |> padding(padding.top = 1, padding.bottom = 1, part = "all")
  if (merge_scenario) ft <- merge_v(ft, j = "Scenario") |> valign(j = "Scenario", valign = "top", part = "body")
  if (!footer_after_fit) {
    ft <- add_footer_lines(ft, notes) |>
      fontsize(size = 7, part = "footer") |>
      font(fontname = "Arial", part = "footer") |>
      align(align = "left", part = "footer")
  }
  title_lines <- strsplit(title, "\n", fixed = TRUE)[[1]]
  ft <- add_header_lines(ft, title_lines) |>
    bold(i = seq_along(title_lines), part = "header") |>
    fontsize(i = seq_along(title_lines), size = 8.5, part = "header") |>
    align(i = seq_along(title_lines), align = "left", part = "header")
  ft <- autofit(ft) |> fit_to_width(max_width = 10.7)
  if (!footer_after_fit) return(ft)
  ft <- add_footer_lines(ft, notes)
  footnote_rows <- grep("^[a-g](,[a-g])? ", notes)
  for (i in footnote_rows) {
    marker <- sub("^([a-g](,[a-g])?) .*", "\\1", notes[i])
    note_text <- sub("^[a-g](,[a-g])? ", "", notes[i])
    ft <- compose(
      ft, i = i, j = seq_len(ncol(display)), part = "footer",
      value = as_paragraph(as_sup(marker), as_chunk(paste0(" ", note_text)))
    )
  }
  ft |> fontsize(size = 7, part = "footer") |>
    font(fontname = "Arial", part = "footer") |>
    align(align = "left", part = "footer")
}

write_table <- function(ft, stem) {
  suppressPackageStartupMessages({ library(flextable); library(officer) })
  flextable::save_as_html(ft, path = paste0(stem, ".html"))
  docx <- paste0(stem, ".docx")
  pdf <- paste0(stem, ".pdf")
  prop <- prop_section(page_size = page_size(orient = "landscape", width = 8.27, height = 11.69),
                       page_margins = page_mar(0.4, 0.4, 0.4, 0.4, 0.2, 0.2))
  flextable::save_as_docx(ft, path = docx, pr_section = prop)
  if (file.exists(pdf) && !file.remove(pdf)) {
    stop("Could not remove the previous PDF before conversion: ", pdf, call. = FALSE)
  }
  status <- system2("env", c("-u", "LD_LIBRARY_PATH", "soffice", "--headless", "--convert-to", "pdf", "--outdir", shQuote(dirname(stem)), shQuote(docx)),
                    stdout = FALSE, stderr = FALSE)
  file.remove(docx)
  if (status != 0L || !file.exists(pdf) || file.info(pdf)$size <= 0) {
    stop("PDF conversion failed for ", stem, call. = FALSE)
  }
}

build_stage <- function() {
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  tab <- collect_table()
  fwrite(tab, file.path(table_dir, "table1_full.csv"))
  header <- c("Scenario", "Region", "Lifetime risk from 65 (%)", "Persistent share (%)", "General tier (%)",
              "Sustained tier (%)", "Heat erosion ratio (%)", "All-temperature erosion ratio (%)", "Net drag from")
  sup <- c("", "", "a", "b", "c", "d", "e", "f", "g")
  main <- tab[level %in% c("europe", "region")]
  write_table(make_flextable(
    format_table1(main), header, sup, table_title, c(table_notes, footnotes),
    footer_after_fit = TRUE
  ),
              file.path(table_dir, "table1"))

  s1 <- data.table(
    Scenario = scenario_labels[as.character(main$ssp)], Region = main$geography,
    `LE65 gain without climate change (months)` = f1(main$g0_months),
    `LE65 gain with climate change (months)` = f1(main$g1_months),
    `All-temperature loss (months)` = f1(main$all_loss_months),
    `Heat loss (months)` = f1(main$heat_loss_ensemble_months),
    `Cold loss (months)` = f1(main$cold_loss_ensemble_months),
    `Heat erosion ratio (%)` = heat_cell(main),
    `All-temperature erosion ratio (%)` = f1(main$all_ratio_shown_pct)
  )
  s1_note <- paste(
    "All quantities use period-mean schedules for the pooled ensemble.",
    "Heat loss + cold loss = all-temperature loss; the remaining-cause contribution is zero because rest rates match between branches.",
    "Ratios use the same period-mean G0 and G1 as Table 1.")
  write_table(make_flextable(s1, names(s1), rep("", ncol(s1)),
                             "Table S1 Losses and gains of LE65, 2020\u201324 to 2095\u201399, in months", s1_note,
                             footer_after_fit = TRUE),
              file.path(table_dir, "tableS1"))
  for (s in 1:3) {
    d <- format_table1(tab[ssp == s & level == "country"], "Country")[, Scenario := NULL]
    write_table(
      make_flextable(d, c("Country", header[3:9]), sup[2:9],
                     sprintf("Table S%d Temperature-related mortality among people aged 65+ by country, %s",
                             s + 1L, scenario_labels[[as.character(s)]]),
                     "Columns as in Table 1. Column e uses per-GCM pooled decompositions for each country.",
                     merge_scenario = FALSE, n_left = 1L),
      file.path(table_dir, sprintf("tableS%d", s + 1L)))
  }
  message("Saved Table 1 and supplements to ", table_dir)
  invisible(tab)
}

#------------------------
# CHECKS
#------------------------

test_stage <- function() {
  tab <- fread(file.path(table_dir, "table1_full.csv"))
  per <- fread(file.path(scratch_dir, "per_gcm_pooled_results.csv"))
  check_results <- list()
  check <- function(cond, msg) {
    if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
    message("PASS: ", msg)
    check_results[[length(check_results) + 1L]] <<- data.table(
      check = msg,
      status = "PASS"
    )
  }
  gradient_test_rates <- rep(0.01, length(decade_age_levels) * length(decade_cause_levels))
  gradient_test <- life_expectancy_gradient(gradient_test_rates)
  gradient_test_indices <- c(1L, 35L, 36L, 100L, length(gradient_test_rates))
  gradient_test_errors <- vapply(gradient_test_indices, function(i) {
    plus <- gradient_test_rates
    minus <- gradient_test_rates
    plus[i] <- plus[i] + 1e-7
    minus[i] <- minus[i] - 1e-7
    numeric <- (decade_life_expectancy(plus) - decade_life_expectancy(minus)) / 2e-7
    abs(gradient_test[i] - numeric)
  }, numeric(1))
  check(max(gradient_test_errors) < 1e-6,
        "analytic LE65 gradient matches finite differences for representative age-cause rates")
  check(nrow(tab) == 105L &&
          nrow(tab[level %in% c("europe", "region")]) == 15L &&
          nrow(tab[level == "country"]) == 90L &&
          uniqueN(tab[level == "country", id]) == 30L,
        "105 rows cover 15 Table 1 rows and 30 countries across three scenarios")
  gcm_columns <- paste0("heat_ratio_", gcm_names)
  check(all(gcm_columns %in% names(tab)) &&
          all(vapply(tab[, ..gcm_columns], function(x) all(is.finite(x)), logical(1))),
        "all 19 per-GCM heat ratios are present and finite")
  eu <- tab[ssp == 3 & level == "europe"]
  check(abs(eu$all_ratio_pct - 4.3) < 0.05 && abs(eu$g0_months - 58) < 0.5,
        "SSP3 Europe erosion ratio about 4.3% with G0 about 58 months")
  check(abs(100 * (eu$g0_months - eu$g1_months) / eu$g0_months - eu$all_ratio_pct) < 1e-9,
        "ratio recomputed by hand from G0 and G1")
  ensemble_gap <- max(abs(tab$heat_loss_ensemble_months +
                            tab$cold_loss_ensemble_months +
                            tab$rest_loss_ensemble_months -
                            tab$all_loss_months))
  gcm_gap <- max(abs(per$additivity_gap_months))
  check(nrow(per) == 3L * 19L * 35L &&
          all(abs(per$additivity_gap_months) < 1e-6) &&
          ensemble_gap < 1e-6,
        sprintf("heat + cold + rest = G0 - G1 for all GCM and ensemble rows (largest gap %.3e months)",
                max(gcm_gap, ensemble_gap)))
  check(max(abs(c(per$rest_loss_months, tab$rest_loss_ensemble_months))) < 1e-9,
        sprintf("rest-cause loss is zero (largest absolute value %.3e months)",
                max(abs(c(per$rest_loss_months, tab$rest_loss_ensemble_months)))))
  col6_gap <- max(abs(tab$all_ratio_pct -
                        100 * (tab$heat_loss_ensemble_months +
                                 tab$cold_loss_ensemble_months) / tab$g0_months))
  gcm_col6_gap <- max(abs(per$all_ratio_pct -
                            100 * (per$heat_loss_months + per$cold_loss_months) /
                              per$g0_months))
  check(max(col6_gap, gcm_col6_gap) < 1e-6,
        sprintf("all-temperature ratios equal 100 * (heat loss + cold loss) / G0 (ensemble gap %.3e points; GCM gap %.3e points)",
                col6_gap, gcm_col6_gap))
  hidden_cells <- tab[hidden == TRUE, .(scenario = scenario_labels[as.character(ssp)], geography)]
  check(all(is.na(tab[hidden == TRUE, all_ratio_shown_pct])) &&
          all(!is.na(tab[hidden == FALSE, all_ratio_shown_pct])),
        "column 6 is hidden only when G0 is below one month")
  message("Column 6 hidden cells: ",
          if (nrow(hidden_cells)) paste(paste(hidden_cells$scenario, hidden_cells$geography, sep = " / "),
                                        collapse = "; ") else "none")
  closure_before <- max(abs(c(
    per$closure_start_n400_months, per$closure_end_n400_months,
    tab$closure_start_n400_months, tab$closure_end_n400_months
  )))
  closure_after <- max(abs(c(
    per$closure_start_final_months, per$closure_end_final_months,
    tab$closure_start_final_months, tab$closure_end_final_months
  )))
  closure_numeric_last <- max(abs(c(
    per$closure_start_numeric_last_months, per$closure_end_numeric_last_months,
    tab$closure_start_numeric_last_months, tab$closure_end_numeric_last_months
  )))
  check(closure_after < 1e-6,
        sprintf("final Horiuchi closure is below 1e-6 months (%.3e)", closure_after))
  message(sprintf(
    "Largest closure: N=400 %.3e months; last DemoDecomp N %.3e months; final %.3e months.",
    closure_before, closure_numeric_last, closure_after
  ))
  fallback_count <- per[
    method_start == "analytic-gradient Gauss-Legendre" |
      method_end == "analytic-gradient Gauss-Legendre", .N
  ]
  message(sprintf("Per-GCM geography rows with analytic-gradient fallback: %d", fallback_count))
  escalations <- fread(file.path(table_dir, "table1_n_escalation.csv"))
  check(all(c("ssp", "gcm", "level", "id", "period", "N",
              "closure_n400_months", "closure_numeric_last_months",
              "closure_final_months", "method", "quadrature_nodes") %in% names(escalations)),
        "ensemble and per-GCM N escalation report is present")
  message(sprintf("Decompositions with N > 400: %d, including %d ensemble decompositions.",
                  nrow(escalations), escalations[gcm == "ENSEMBLE", .N]))
  if (nrow(escalations)) print(escalations)
  simon <- c(EUROPE = 4.3, Northern = -0.1, Western = 3.4, Eastern = 4.8, Southern = 8.9)
  got <- tab[ssp == 3 & level %in% c("europe", "region"), setNames(all_ratio_pct, id)][names(simon)]
  print(round(rbind(simon, got, diff = got - simon), 3))
  check(all(is.finite(got)) && all(abs(got - simon) <= 0.5),
        "SSP3 column 6 is within 0.5 percentage points of Simon's provisional values")
  check(all(abs(tab$general_pct + tab$sustained_pct + tab$episodic_pct - 100) < 1e-9),
        "general + sustained + episodic = 100% (tolerance 1e-9)")
  meta <- city_meta()[region == "Northern"]
  cities <- head(meta$URAU_CODE, 5L)
  dem <- read_demography(3L, cities)
  pooling_members <- data.table(geo_id = cities, level = "region", id = "Northern")
  city_members <- data.table(geo_id = cities, level = "city", id = cities)
  pooled_schedule <- period_mean_schedule(pooled_case(3L, gcm_names[8], dem, pooling_members))
  city_schedule <- period_mean_schedule(pooled_case(3L, gcm_names[8], dem, city_members))
  pooled_metric <- pooled_schedule[, compute_period_metrics(.SD, "pooling test pooled"),
                                   by = .(level, id)]
  city_metrics <- city_schedule[, compute_period_metrics(
    .SD, paste("pooling test city", level, id)
  ), by = .(level, id)]
  pooled_ratio <- pooled_metric$all_ratio_pct
  mean_city_ratio <- mean(city_metrics$all_ratio_pct)
  message(sprintf("Northern SSP3, five-city check: mean city ratio %.3f%%, pooled ratio %.3f%%",
                  mean_city_ratio, pooled_ratio))
  check(abs(mean_city_ratio - pooled_ratio) > 1e-3,
        "average of five city ratios differs from the pooled period-mean ratio")
  fixture_path <- file.path(
    batch_root, "ssp3", "ES001C", "GFDL_ESM4", "checks", "01_attribution_checks.csv"
  )
  fixture <- fread(fixture_path)
  check(fixture[check_name == "masselot2023_fixture", status] == "PASS",
        "Masselot 2023 fixture passes")
  message("Largest |GCM median - ensemble| per row group (points):")
  comparison <- tab[level %in% c("europe", "region"),
                    .(ssp, id, gcm_minus_ensemble_pp = round(gcm_minus_ensemble_pp, 3))]
  print(comparison)
  material_diff <- comparison[abs(gcm_minus_ensemble_pp) > 0.3]
  if (nrow(material_diff)) {
    message("GCM median versus ensemble differences above 0.3 points:")
    print(material_diff)
  }
  fwrite(rbindlist(check_results), file.path(table_dir, "table1_checks.csv"))
}

if (!nzchar(Sys.getenv("TABLE1_SOURCE_ONLY"))) {
  if (stage == "gcm_time") run_gcm_stage(TRUE)
  if (stage %in% c("gcm", "all")) run_gcm_stage()
  if (stage %in% c("build", "all")) build_stage()
  if (stage %in% c("test", "all")) test_stage()
}
