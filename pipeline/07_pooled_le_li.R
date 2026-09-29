#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# Pipeline Part 07: Pooled LE65 / LI65+ analysis above the city level
#   Implements the aggregation decision from the 25 Sep 2026 meeting: ANs are
#   estimated at city level, then population and deaths by cause are summed
#   across cities before the life table or decomposition is constructed. City
#   ERFs, city LE/LI values and city decomposition contributions are never
#   averaged to obtain country, regional or European results.
#
#   One invocation analyses one geography selected with GEO_LEVEL and GEO_ID:
#     GEO_LEVEL=europe GEO_ID=EUROPE
#     GEO_LEVEL=region GEO_ID=Southern
#     GEO_LEVEL=country GEO_ID=ES
#
################################################################################

source("pipeline/00_pkg_params.R")

geo_level <- tolower(Sys.getenv("GEO_LEVEL", "europe"))
geo_id <- Sys.getenv("GEO_ID", if (geo_level == "europe") "EUROPE" else "")
obj1_file <- Sys.getenv("OBJECT1_FILE", file.path("results/europe/collected", "object1_dataset.parquet"))

if (!geo_level %in% c("europe", "region", "country")) {
  stop("GEO_LEVEL must be europe, region or country.", call. = FALSE)
}
if (geo_level != "europe" && !nzchar(geo_id)) stop("GEO_ID is required.", call. = FALSE)
if (!file.exists(obj1_file)) stop("Object 1 not found: ", obj1_file, call. = FALSE)

meta <- unique(fread("data/city_results.csv")[, .(
  city = URAU_CODE, city_name = LABEL, country = CNTR_CODE, region
)])
if (geo_level == "region" && !geo_id %in% meta$region) stop("Unknown region: ", geo_id, call. = FALSE)
if (geo_level == "country" && !geo_id %in% meta$country) stop("Unknown country: ", geo_id, call. = FALSE)

geo_label <- switch(geo_level,
  europe = "Europe",
  region = paste(geo_id, "Europe"),
  country = geo_id
)
message(sprintf("\n[07] Pooled %s analysis for %s, %s (N = %d)...",
  geo_level, geo_label, ssplabs[ssp_name], N_HORIUCHI))

pooled_file <- file.path(out_dir, "07_pooled_master.csv")
levels_file <- file.path(out_dir, "07_le_li_levels.csv")
between_file <- file.path(out_dir, "07_between_branch_decomposition.csv")
within_file <- file.path(out_dir, "07_within_branch_period_decomposition.csv")
age98_file <- file.path(out_dir, "07_age98_share.csv")
checks_file <- file.path(check_dir, "07_pooled_checks.csv")
failures_file <- file.path(check_dir, "07_pooled_failures.csv")

#----- Sum city counts before constructing mortality rates

ds <- open_dataset(obj1_file)
q <- ds %>% filter(ssp == !!as.integer(ssp_name))
if (geo_level == "region") q <- q %>% filter(region == !!geo_id)
if (geo_level == "country") q <- q %>% filter(country == !!geo_id)

cities_included <- q %>% select(city) %>% distinct() %>% collect()
if (!nrow(cities_included)) stop("No completed cities found for this geography.", call. = FALSE)
expected_cities <- switch(geo_level,
  europe = meta$city,
  region = meta[region == geo_id, city],
  country = meta[country == geo_id, city]
)
coverage_ok <- setequal(cities_included$city, expected_cities)
if (!coverage_ok) {
  missing_cities <- setdiff(expected_cities, cities_included$city)
  extra_cities <- setdiff(cities_included$city, expected_cities)
  stop(sprintf(
    "City coverage mismatch for %s: %d missing (%s); %d unexpected (%s).",
    geo_label, length(missing_cities), paste(missing_cities, collapse = ", "),
    length(extra_cities), paste(extra_cities, collapse = ", ")
  ), call. = FALSE)
}

deaths <- q %>%
  group_by(scenario, year, age, cause) %>%
  summarise(deaths_cause = sum(deaths, na.rm = TRUE)) %>%
  collect() %>%
  as.data.table()

# Population is repeated for every cause in Object 1. Select one cause before
# summing so each city-age-year contributes exactly once.
pop <- q %>%
  filter(cause == "rest") %>%
  group_by(scenario, year, age) %>%
  summarise(pop = sum(pop, na.rm = TRUE)) %>%
  collect() %>%
  as.data.table()

setnames(deaths, "scenario", "branch")
setnames(pop, "scenario", "branch")
cause_dt <- merge(deaths, pop, by = c("branch", "year", "age"), all = TRUE)
setorder(cause_dt, branch, year, cause, age)
cause_dt[, `:=`(
  geo_level = geo_level,
  geo_id = geo_id,
  label = geo_label,
  ssp = as.integer(ssp_name)
)]
setcolorder(cause_dt, c("geo_level", "geo_id", "label", "ssp", "branch", "year", "age", "cause", "pop", "deaths_cause"))
fwrite(cause_dt, pooled_file)

expected_rows <- length(branch_levels) * length(future_years) * length(age_levels) * length(cause_levels)
if (nrow(cause_dt) != expected_rows) stop("Pooled input grid is incomplete or duplicated.", call. = FALSE)
if (anyNA(cause_dt) || any(!is.finite(cause_dt$pop)) || any(!is.finite(cause_dt$deaths_cause))) {
  stop("Pooled input contains missing or non-finite values.", call. = FALSE)
}
if (any(cause_dt$pop <= 0) || any(cause_dt$deaths_cause < -1e-9)) {
  stop("Pooled input contains nonpositive population or negative deaths.", call. = FALSE)
}

cause_dt[, mx_cause := deaths_cause / pop]
cause_dt[, mx_total := sum(mx_cause), by = .(branch, year, age)]

#----- Lloyd et al. (2024) / Aburto et al. (2022) functions (unchanged)

life_expectancy_from_mx_65plus <- function(mx, x, nx = c(rep(1, 100 - 65), Inf), age = 0) {
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  ex[age + 1]
}

sd_from_mx_fun_65_plus <- function(mx, x, nx = c(rep(1, 100 - 65), Inf), age = 0) {
  x_conditional <- x - 65
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  sqrt(sum(dx * (x_conditional + 0.5 - ex[age + 1])^2))
}

life_expectancy_cod <- function(mx.cod, x, nx = c(rep(1, 100 - 65), Inf), cond_age = 0) {
  dim(mx.cod) <- c(length(x), length(mx.cod) / length(x))
  life_expectancy_from_mx_65plus(rowSums(mx.cod), x, nx, cond_age)
}

sd.cod.fun.65plus <- function(mx.cod, x, nx, cond_age = 0) {
  dim(mx.cod) <- c(length(x), length(mx.cod) / length(x))
  sd_from_mx_fun_65_plus(rowSums(mx.cod), x, nx, cond_age)
}

cause_vector <- function(dt) {
  ordered <- dt[order(match(cause, cause_levels), age)]
  if (nrow(ordered) != length(age_levels) * length(cause_levels)) stop("Cause-age grid mismatch.", call. = FALSE)
  ordered$mx_cause
}

horiuchi_pair <- function(mx1, mx2) {
  list(
    le = as.vector(horiuchi(func = life_expectancy_cod, pars1 = mx1, pars2 = mx2,
      N = N_HORIUCHI, x = age_levels, nx = nx, cond_age = 0)),
    li = as.vector(horiuchi(func = sd.cod.fun.65plus, pars1 = mx1, pars2 = mx2,
      N = N_HORIUCHI, x = age_levels, nx = nx, cond_age = 0))
  )
}

input_hash <- unname(tools::md5sum(pooled_file))
cache_dir <- file.path(out_dir, "07_cache", input_hash)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cached <- function(key, expr) {
  f <- file.path(cache_dir, paste0(key, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  val <- force(expr)
  saveRDS(val, f)
  val
}

grid_dt <- function(values) data.table(
  age = rep(age_levels, length(cause_levels)),
  cause = rep(cause_levels, each = length(age_levels)),
  contribution = values
)

#----- Annual pooled LE65 and LI65+ levels

levels_dt <- cause_dt[cause == cause_levels[1], .(
  LE65 = life_expectancy_from_mx_65plus(mx_total[order(age)], age_levels, nx, 0),
  LI65 = sd_from_mx_fun_65_plus(mx_total[order(age)], age_levels, nx, 0)
), by = .(branch, year)]
levels_dt[, `:=`(geo_level = geo_level, geo_id = geo_id, label = geo_label, ssp = as.integer(ssp_name))]
setcolorder(levels_dt, c("geo_level", "geo_id", "label", "ssp", "branch", "year", "LE65", "LI65"))

#----- Decomposition on 5-year-period mean mortality schedules

period_len <- perlen
cause_dt[, period := 2020L + ((year - 2020L) %/% period_len) * period_len]
period_mx <- cause_dt[, .(mx_cause = mean(mx_cause)), by = .(branch, period, age, cause)]
periods_all <- sort(unique(period_mx$period))
f_le <- function(v) life_expectancy_cod(v, x = age_levels, nx = nx)
f_li <- function(v) sd.cod.fun.65plus(v, x = age_levels, nx = nx)

between <- list(); between_checks <- list()
for (p in periods_all) {
  mx_wo <- cause_vector(period_mx[branch == "without_cc" & period == p])
  mx_w <- cause_vector(period_mx[branch == "with_cc" & period == p])
  h <- cached(sprintf("between_%d_N%d", p, N_HORIUCHI), horiuchi_pair(mx_wo, mx_w))
  g <- grid_dt(h$le)
  setnames(g, "contribution", "le_contribution")
  g[, li_contribution := h$li]
  between[[length(between) + 1L]] <- cbind(data.table(
    geo_level = geo_level, geo_id = geo_id, label = geo_label, ssp = as.integer(ssp_name),
    period = sprintf("%d-%d", p, p + period_len - 1L)
  ), g)
  between_checks[[length(between_checks) + 1L]] <- data.table(
    period = p,
    le_closure_error = sum(h$le) - (f_le(mx_w) - f_le(mx_wo)),
    li_closure_error = sum(h$li) - (f_li(mx_w) - f_li(mx_wo)),
    max_abs_rest = max(abs(g[cause == "rest", c(le_contribution, li_contribution)]))
  )
}
between <- rbindlist(between)
between_checks <- rbindlist(between_checks)

within <- list(); within_checks <- list()
for (b in branch_levels) {
  for (i in seq_len(length(periods_all) - 1L)) {
    p0 <- periods_all[i]; p1 <- periods_all[i + 1L]
    m0 <- cause_vector(period_mx[branch == b & period == p0])
    m1 <- cause_vector(period_mx[branch == b & period == p1])
    h <- cached(sprintf("within_%s_%d_%d_N%d", b, p0, p1, N_HORIUCHI), horiuchi_pair(m0, m1))
    g <- grid_dt(h$le)
    setnames(g, "contribution", "le_contribution")
    g[, li_contribution := h$li]
    within[[length(within) + 1L]] <- cbind(data.table(
      geo_level = geo_level, geo_id = geo_id, label = geo_label, ssp = as.integer(ssp_name), branch = b,
      period_from = sprintf("%d-%d", p0, p0 + period_len - 1L),
      period_to = sprintf("%d-%d", p1, p1 + period_len - 1L)
    ), g)
    within_checks[[length(within_checks) + 1L]] <- data.table(
      branch = b, period_from = p0,
      le_closure_error = sum(h$le) - (f_le(m1) - f_le(m0)),
      li_closure_error = sum(h$li) - (f_li(m1) - f_li(m0))
    )
  }
}
within <- rbindlist(within)
within_checks <- rbindlist(within_checks)

#----- Simon's ages 98+ diagnostic

age98 <- rbind(
  between[cause != "rest" & age >= 90L, .(
    contribution_90plus = sum(le_contribution),
    contribution_98plus = sum(le_contribution[age >= 98L]),
    abs_contribution_90plus = sum(abs(le_contribution)),
    abs_contribution_98plus = sum(abs(le_contribution[age >= 98L]))
  ), by = .(period, cause)][, measure := "LE65"],
  between[cause != "rest" & age >= 90L, .(
    contribution_90plus = sum(li_contribution),
    contribution_98plus = sum(li_contribution[age >= 98L]),
    abs_contribution_90plus = sum(abs(li_contribution)),
    abs_contribution_98plus = sum(abs(li_contribution[age >= 98L]))
  ), by = .(period, cause)][, measure := "LI65"]
)
age98[, `:=`(
  signed_share_98plus = fifelse(abs(contribution_90plus) > 0, contribution_98plus / contribution_90plus, NA_real_),
  absolute_share_98plus = fifelse(abs_contribution_90plus > 0, abs_contribution_98plus / abs_contribution_90plus, NA_real_),
  geo_level = geo_level, geo_id = geo_id, label = geo_label, ssp = as.integer(ssp_name)
)]
setcolorder(age98, c("geo_level", "geo_id", "label", "ssp", "period", "measure", "cause",
  "contribution_90plus", "contribution_98plus", "signed_share_98plus", "absolute_share_98plus"))

#----- Invariant checks

pop_branch <- dcast(unique(cause_dt[cause == "rest", .(branch, year, age, pop)]), year + age ~ branch, value.var = "pop")
rest_branch <- dcast(cause_dt[cause == "rest", .(branch, year, age, deaths_cause)], year + age ~ branch, value.var = "deaths_cause")
max_between_error <- max(abs(c(between_checks$le_closure_error, between_checks$li_closure_error)))
max_within_error <- max(abs(c(within_checks$le_closure_error, within_checks$li_closure_error)))

checks <- data.table(
  check_name = c("city_coverage", "pooled_grid_complete", "no_negative_deaths", "population_same_between_branches",
    "rest_same_between_branches", "levels_complete", "between_branch_closure", "between_branch_rest_zero",
    "within_period_closure", "contributions_finite"),
  status = c(
    if (coverage_ok) "PASS" else "FAIL",
    if (nrow(cause_dt) == expected_rows) "PASS" else "FAIL",
    if (all(cause_dt$deaths_cause >= -1e-9)) "PASS" else "FAIL",
    if (max(abs(pop_branch$with_cc - pop_branch$without_cc)) <= 1e-6) "PASS" else "FAIL",
    if (max(abs(rest_branch$with_cc - rest_branch$without_cc)) <= 1e-6) "PASS" else "FAIL",
    if (nrow(levels_dt) == length(branch_levels) * length(future_years) && all(is.finite(c(levels_dt$LE65, levels_dt$LI65)))) "PASS" else "FAIL",
    if (max_between_error <= closure_tol) "PASS" else "FAIL",
    if (max(between_checks$max_abs_rest) <= 1e-12) "PASS" else "FAIL",
    if (max_within_error <= closure_tol) "PASS" else "FAIL",
    if (all(is.finite(c(between$le_contribution, between$li_contribution, within$le_contribution, within$li_contribution)))) "PASS" else "FAIL"
  ),
  value = c(
    sprintf("%d/%d cities", nrow(cities_included), length(expected_cities)),
    sprintf("%d rows", nrow(cause_dt)),
    sprintf("min=%0.3e", min(cause_dt$deaths_cause)),
    sprintf("max_abs_diff=%0.3e", max(abs(pop_branch$with_cc - pop_branch$without_cc))),
    sprintf("max_abs_diff=%0.3e", max(abs(rest_branch$with_cc - rest_branch$without_cc))),
    sprintf("%d rows", nrow(levels_dt)),
    sprintf("max_abs_error=%0.3e", max_between_error),
    sprintf("max_abs_rest=%0.3e", max(between_checks$max_abs_rest)),
    sprintf("max_abs_error=%0.3e", max_within_error),
    "all finite"
  ),
  threshold = c(sprintf("= %d", length(expected_cities)), sprintf("= %d", expected_rows), ">= 0", "<= 1e-6", "<= 1e-6", "160 finite rows",
    sprintf("<= %g", closure_tol), "<= 1e-12", sprintf("<= %g", closure_tol), "finite")
)

fwrite(checks, checks_file)
failed <- checks[status == "FAIL"]
if (nrow(failed)) {
  fwrite(failed, failures_file)
  stop(sprintf("07_pooled_le_li.R failed %d invariant(s); see %s", nrow(failed), failures_file), call. = FALSE)
}
if (file.exists(failures_file)) file.remove(failures_file)
invisible(file.create(failures_file))

fwrite(levels_dt, levels_file)
fwrite(between, between_file)
fwrite(within, within_file)
fwrite(age98, age98_file)

message(sprintf("Saved pooled results for %s (%d cities).", geo_label, nrow(cities_included)))
message("Saved checks to ", checks_file)

