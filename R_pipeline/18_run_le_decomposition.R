#!/usr/bin/env Rscript
# ==============================================================================
# 18_run_le_decomposition.R
#
# Runs the Horiuchi life-expectancy decomposition on Lloyd-compatible
# future life-table inputs produced by script 16.
#
# For each scenario combination (geo × SSP × GCM × sim × sex),
# decomposes changes in remaining LE at 65 into contributions by
# age and cause (extr_cold, mod_cold, mod_heat, extr_heat, rest)
# for each consecutive year pair.
#
# Based on Code_1.R from Lloyd et al. (2024), Environment International 193.
# doi: 10.1016/j.envint.2024.109050
# Life table functions are preserved exactly as published.
#
# Inputs:
#   Wide-format CSV from script 16, with columns:
#     geo_id, geo_name, year, ssp, gcm, sim, sex, reth, age, pop,
#     total_deaths, extr_cold, mod_cold, mod_heat, extr_heat, rest
#
#   IMPORTANT: years must be consecutive within each scenario group.
#   The Horiuchi decomposition compares year y-1 to year y.
#
# Outputs:
#   decomposition_results.csv  — contributions to LE change at 65 by
#                                 age, cause, and year-step
#   le_by_year.csv             — computed LE at 65 for each scenario × year
#   decomposition_checks.csv   — accuracy checks per year-step
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(DemoDecomp)
  library(parallel)
})

message("\n[18] Running LE decomposition on future Lloyd-compatible inputs...")

# --- Configuration via environment variables ----------------------------------

input_file     <- trimws(Sys.getenv("INPUT_FILE",     unset = "results/le_li_input/future_temperature_deaths_wide.csv"))
output_file    <- trimws(Sys.getenv("OUTPUT_FILE",    unset = "results/le_li_decomposition/decomposition_results.csv"))
le_file        <- trimws(Sys.getenv("LE_FILE",        unset = "results/le_li_decomposition/le_by_year.csv"))
check_file     <- trimws(Sys.getenv("CHECK_FILE",     unset = "results/le_li_decomposition/decomposition_checks.csv"))
country_filter <- trimws(Sys.getenv("COUNTRY_FILTER", unset = ""))
ssp_filter     <- trimws(Sys.getenv("SSP_FILTER",     unset = ""))
gcm_filter     <- trimws(Sys.getenv("GCM_FILTER",     unset = ""))
sim_filter     <- trimws(Sys.getenv("SIM_FILTER",     unset = ""))
year_filter    <- trimws(Sys.getenv("YEAR_FILTER",    unset = ""))
n_cores_env    <- suppressWarnings(as.integer(trimws(Sys.getenv("N_CORES", unset = ""))))
n_horiuchi_env <- suppressWarnings(as.integer(trimws(Sys.getenv("N_HORIUCHI", unset = ""))))

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(le_file),     recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(check_file),  recursive = TRUE, showWarnings = FALSE)

# --- Constants ----------------------------------------------------------------

COND_SURV  <- 65L
CAUSE_COLS <- c("extr_cold", "mod_cold", "mod_heat", "extr_heat", "rest")
N_HORIUCHI <- if (is.na(n_horiuchi_env) || n_horiuchi_env < 1L) 50L else n_horiuchi_env

# --- Filter helpers -----------------------------------------------------------

parse_int_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

parse_chr_filter <- function(value) {
  if (!nzchar(value)) return(NULL)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

wanted_countries <- parse_chr_filter(country_filter)
wanted_ssps      <- parse_int_filter(ssp_filter)
wanted_gcms      <- parse_chr_filter(gcm_filter)
wanted_sims      <- parse_int_filter(sim_filter)
wanted_years     <- parse_int_filter(year_filter)

# ==============================================================================
# Life table functions — from Lloyd et al. (2024), preserved as published
# ==============================================================================

life_expectancy_from_mx <- function(mx, x,
                                    nx = c(rep(1, 100 - COND_SURV), Inf),
                                    age = 0) {
  px <- exp(-mx * nx)
  qx <- 1 - px
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  return(ex[age + 1])
}

life_expectancy_cod <- function(mx.cod, x,
                                nx = c(rep(1, 100 - COND_SURV), Inf),
                                cond_age = 0) {
  dim(mx.cod) <- c(length(x), length(mx.cod) / length(x))
  mx <- rowSums(mx.cod)
  life_expectancy_from_mx(mx, x, nx, cond_age)
}

# ==============================================================================
# Load and prepare data
# ==============================================================================

message("  Reading input: ", input_file)
wide <- fread(input_file)

required_cols <- c("geo_id", "geo_name", "year", "ssp", "gcm", "sim",
                   "sex", "reth", "age", "pop", "total_deaths", CAUSE_COLS)
missing <- setdiff(required_cols, names(wide))
if (length(missing)) {
  stop("Input missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Apply filters
if (!is.null(wanted_countries)) wide <- wide[geo_id %in% wanted_countries]
if (!is.null(wanted_years))     wide <- wide[year %in% wanted_years]
if (!is.null(wanted_ssps))      wide <- wide[ssp %in% wanted_ssps]
if (!is.null(wanted_gcms))      wide <- wide[gcm %in% wanted_gcms]
if (!is.null(wanted_sims))      wide <- wide[sim %in% wanted_sims]

if (!nrow(wide)) {
  stop("All rows filtered out; nothing to process.", call. = FALSE)
}

# Filter to ages 65+
wide <- wide[age >= COND_SURV]

# Melt to long format
# data.table::melt groups by measure variable: all ages for cause 1 first,
# then all ages for cause 2, etc. This ordering is required by the
# Horiuchi decomposition's vector ↔ matrix reshape.
long <- melt(
  wide,
  id.vars = c("geo_id", "geo_name", "year", "ssp", "gcm", "sim",
              "sex", "reth", "age", "pop", "total_deaths"),
  measure.vars = CAUSE_COLS,
  variable.name = "cause",
  value.name = "deaths_cause"
)

# Compute cause-specific and total mortality rates
long[, mx_cause := deaths_cause / pop]
long[, mx_total := total_deaths / pop]

# Explicit sort: within each scenario, order by year → cause → age
# This ensures the mx_cause vector is correctly structured for horiuchi()
setorder(long, geo_id, ssp, gcm, sim, sex, reth, year, cause, age)

# ==============================================================================
# Validate consecutive years within each scenario
# ==============================================================================

scenario_keys <- c("geo_id", "ssp", "gcm", "sim", "sex", "reth")
year_check <- long[, .(years = list(sort(unique(year)))),
                   by = mget(scenario_keys)]

for (i in seq_len(nrow(year_check))) {
  yrs <- year_check$years[[i]]
  if (length(yrs) < 2L) {
    stop(sprintf(
      "Scenario %s/SSP%s/%s/sim%s/%s has only %d year(s); need >= 2.",
      year_check$geo_id[i], year_check$ssp[i], year_check$gcm[i],
      year_check$sim[i], year_check$sex[i], length(yrs)
    ), call. = FALSE)
  }
  gaps <- diff(yrs)
  if (any(gaps != 1L)) {
    stop(sprintf(
      paste0("Scenario %s/SSP%s/%s/sim%s/%s has non-consecutive years: %s.\n",
             "  The Horiuchi decomposition requires consecutive annual data.\n",
             "  Re-run scripts 17 + 16 without YEAR_FILTER or with consecutive years."),
      year_check$geo_id[i], year_check$ssp[i], year_check$gcm[i],
      year_check$sim[i], year_check$sex[i],
      paste(yrs, collapse = ",")
    ), call. = FALSE)
  }
}

n_scenarios <- nrow(year_check)
message(sprintf("  Scenarios: %d | Years: %d–%d | Consecutive: ✓",
                n_scenarios, min(long$year), max(long$year)))

# ==============================================================================
# Decomposition function
# ==============================================================================

decompose_group <- function(DT) {
  x     <- sort(unique(DT$age))
  years <- sort(unique(DT$year))
  nx    <- c(rep(1, 100 - COND_SURV), Inf)
  causes <- levels(DT$cause)
  if (is.null(causes)) causes <- unique(as.character(DT$cause))

  n_steps <- length(years) - 1L
  results <- vector("list", n_steps)

  for (j in seq_len(n_steps)) {
    y0 <- years[j]
    y1 <- years[j + 1L]

    # Extract mx vectors: ordered by cause then age (from setorder above)
    mx1 <- DT[year == y0, mx_cause]
    mx2 <- DT[year == y1, mx_cause]

    # Horiuchi decomposition (N = number of integration steps)
    hor <- horiuchi(
      func = life_expectancy_cod,
      pars1 = mx1,
      pars2 = mx2,
      N = N_HORIUCHI,
      x = x, nx = nx, cond_age = 0
    )

    # Reshape vector → matrix (ages × causes), R fills by column
    dim(hor) <- c(length(x), length(hor) / length(x))

    results[[j]] <- data.table(
      age          = rep(x, length(causes)),
      cause        = rep(causes, each = length(x)),
      contribution = as.vector(hor),
      year_from    = y0,
      year_to      = y1
    )
  }

  rbindlist(results)
}

# ==============================================================================
# Run decomposition across all scenarios
# ==============================================================================

message("  Running Horiuchi decomposition (N = ", N_HORIUCHI, ")...")

num_cores <- if (is.na(n_cores_env) || n_cores_env < 1L) {
  min(4L, max(1L, detectCores() - 1L))
} else {
  n_cores_env
}

# Build scenario groups
group_keys <- unique(long[, mget(c(scenario_keys, "geo_name"))])
group_list <- lapply(seq_len(nrow(group_keys)), function(i) {
  k <- group_keys[i]
  long[geo_id == k$geo_id & ssp == k$ssp & gcm == k$gcm &
       sim == k$sim & sex == k$sex & reth == k$reth]
})

results_list <- mclapply(seq_along(group_list), function(i) {
  k <- group_keys[i]
  res <- decompose_group(group_list[[i]])
  if (nrow(res)) {
    res[, `:=`(
      geo_id   = k$geo_id,
      geo_name = k$geo_name,
      ssp      = k$ssp,
      gcm      = k$gcm,
      sim      = k$sim,
      sex      = k$sex,
      reth     = k$reth
    )]
  }
  res
}, mc.cores = num_cores)

decomp <- rbindlist(results_list, use.names = TRUE)

setcolorder(decomp, c(
  "geo_id", "geo_name", "ssp", "gcm", "sim", "sex", "reth",
  "year_from", "year_to", "age", "cause", "contribution"
))
setorder(decomp, geo_id, ssp, gcm, sim, sex, reth, year_from, cause, age)

# ==============================================================================
# Compute LE at 65 for each scenario × year (validation + downstream use)
# ==============================================================================

le_table <- long[, {
  x  <- sort(unique(age))
  nx <- c(rep(1, 100 - COND_SURV), Inf)
  # Sum mx across causes to get total mx, ordered by age
  mx_total_vec <- .SD[, .(mx = sum(mx_cause)), by = age][order(age)]$mx
  le_65 <- life_expectancy_from_mx(mx_total_vec, x, nx, age = 0)
  .(le_65 = le_65)
}, by = .(geo_id, geo_name, ssp, gcm, sim, sex, reth, year)]

setorder(le_table, geo_id, ssp, gcm, sim, sex, reth, year)

# ==============================================================================
# Checks: verify decomposition sums match direct LE changes
# ==============================================================================

checks <- decomp[, .(
  total_contribution = sum(contribution),
  n_ages             = uniqueN(age),
  n_causes           = uniqueN(cause),
  any_na             = anyNA(contribution)
), by = .(geo_id, geo_name, ssp, gcm, sim, sex, reth, year_from, year_to)]

# Merge LE for year_from
checks <- merge(
  checks,
  le_table[, .(geo_id, ssp, gcm, sim, sex, reth, year, le_65)],
  by.x = c("geo_id", "ssp", "gcm", "sim", "sex", "reth", "year_from"),
  by.y = c("geo_id", "ssp", "gcm", "sim", "sex", "reth", "year"),
  all.x = TRUE
)
setnames(checks, "le_65", "le_65_from")

# Merge LE for year_to
checks <- merge(
  checks,
  le_table[, .(geo_id, ssp, gcm, sim, sex, reth, year, le_65)],
  by.x = c("geo_id", "ssp", "gcm", "sim", "sex", "reth", "year_to"),
  by.y = c("geo_id", "ssp", "gcm", "sim", "sex", "reth", "year"),
  all.x = TRUE
)
setnames(checks, "le_65", "le_65_to")

checks[, le_change   := le_65_to - le_65_from]
checks[, decomp_error := total_contribution - le_change]

setorder(checks, geo_id, ssp, gcm, sim, sex, reth, year_from)

# ==============================================================================
# Save outputs
# ==============================================================================

fwrite(decomp,   output_file)
fwrite(le_table, le_file)
fwrite(checks,   check_file)

message("  Saved decomposition results to ", output_file)
message("  Saved LE by year to ", le_file)
message("  Saved checks to ", check_file)
message(sprintf(
  "  Decomp rows: %d | Scenarios: %d | Year-steps: %d | Any NA: %s",
  nrow(decomp),
  n_scenarios,
  nrow(checks),
  any(checks$any_na)
))
message(sprintf(
  "  Max decomp error: %g | LE at 65 range: %.2f to %.2f years",
  max(abs(checks$decomp_error), na.rm = TRUE),
  min(le_table$le_65),
  max(le_table$le_65)
))
