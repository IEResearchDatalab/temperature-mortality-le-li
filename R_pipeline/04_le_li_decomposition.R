#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# R Pipeline Part 04: Life tables, LE65 / LI65+ levels and decompositions
#   Life-table, LE and SD functions are those of Lloyd et al. (2024) Code_1.R /
#   Code_2.R (modified from Aburto et al. 2022). Outputs (the three data objects
#   agreed at the 10 Sep 2026 meeting, Objects 2 and 3 written here):
#     04_le_li_levels.csv                LE65 and LI65+ by branch x year   (Object 3)
#     04_le_decomposition.csv            Horiuchi contributions to year-on-year
#     04_li_decomposition.csv              change, by branch x age x cause (Object 2)
#     04_between_branch_decomposition.csv  with-CC minus without-CC difference,
#                                        by 5-year period x age x cause
#     04_within_branch_period_decomposition.csv  change between consecutive
#                                        5-year-period means, by branch x age x cause
#   Horiuchi steps are cached in <out_dir>/04_cache/<md5 of the master
#   table>/ so an interrupted run resumes where it stopped.
#
################################################################################

source("R_pipeline/00_pkg_params.R")

message(sprintf("\n[04] Running %s LE/LI decomposition (N = %d)...", city_name, N_HORIUCHI))

# The cache is keyed on the md5 of the master table, so any upstream change
# starts a fresh cache instead of reusing stale contributions.
master_file <- file.path(out_dir, "03_master_table.csv")
cache_dir <- file.path(out_dir, "04_cache", unname(tools::md5sum(master_file)))
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

le_file <- file.path(out_dir, "04_le_decomposition.csv")
li_file <- file.path(out_dir, "04_li_decomposition.csv")
levels_file <- file.path(out_dir, "04_le_li_levels.csv")
between_file <- file.path(out_dir, "04_between_branch_decomposition.csv")
within_period_file <- file.path(out_dir, "04_within_branch_period_decomposition.csv")
checks_file <- file.path(check_dir, "04_le_li_decomposition_checks.csv")
steps_file <- file.path(check_dir, "04_le_li_decomposition_steps.csv")
failures_file <- file.path(check_dir, "04_le_li_decomposition_failures.csv")
fig_file <- file.path(fig_dir, "04_le_li_decomposition_diagnostic.png")

# Between-branch decomposition is run on mean mortality rates over the 5-year
# periods of the demographic projections (annual with - without differences
# mix the climate signal with single-year weather).
period_len <- perlen

master <- fread(master_file)
master <- master[geo_id == city_id & gcm == gcm_name]

if (!nrow(master)) stop("Master table for the city is missing; run 03_master_table.R first.", call. = FALSE)

master <- master[branch %in% branch_levels & age %in% age_levels]

if (any(!is.finite(master$pop)) || any(!is.finite(master$death)) || any(!is.finite(master$an)) || any(!is.finite(master$rest))) {
  stop("Master table contains non-finite required values.", call. = FALSE)
}

if (any(master$rest < -1e-9)) {
  stop("Master table contains negative rest mortality.", call. = FALSE)
}

cause_dt <- master[, .(
  deaths_cause = deaths_component,
  pop = first(pop),
  total_deaths = first(death)
), by = .(branch, year, age, cause = range)]

setorder(cause_dt, branch, year, cause, age)
cause_dt[, mx_cause := deaths_cause / pop]
# Total mortality must equal the sum of cause-specific rates so the LE/LI
# closure target matches what the Horiuchi decomposition of mx_cause measures.
# Using raw `death` here would disagree with `deaths_cause` for the with_cc
# branch, whose components sum to `adjusted_death`, not `death`.
cause_dt[, mx_total := sum(mx_cause), by = .(branch, year, age)]

if (any(!cause_dt$cause %in% cause_levels)) {
  stop("Unexpected cause labels in decomposition input.", call. = FALSE)
}
if (anyDuplicated(cause_dt, by = c("branch", "year", "age", "cause"))) {
  stop("Decomposition input has duplicate branch/year/age/cause keys.", call. = FALSE)
}
if (nrow(cause_dt) != length(branch_levels) * length(sort(unique(master$year))) * length(age_levels) * length(cause_levels)) {
  stop("Decomposition input grid is incomplete or oversized.", call. = FALSE)
}

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
  mx <- rowSums(mx.cod)
  life_expectancy_from_mx_65plus(mx, x, nx, cond_age)
}

sd.cod.fun.65plus <- function(mx.cod, x, nx, cond_age = 0) {
  dim(mx.cod) <- c(length(x), length(mx.cod) / length(x))
  mx <- rowSums(mx.cod)
  sd_from_mx_fun_65_plus(mx, x, nx, cond_age)
}

# mx vector ordered cause-major, age-minor (as in Lloyd Code_1.R)
cause_vector <- function(dt) {
  ordered <- dt[order(match(cause, cause_levels), age)]
  if (!identical(as.integer(unique(ordered$age)), age_levels)) stop("Age grid mismatch.", call. = FALSE)
  ordered$mx_cause
}

horiuchi_pair <- function(mx1, mx2) {
  hor_le <- horiuchi(func = life_expectancy_cod, pars1 = mx1, pars2 = mx2, N = N_HORIUCHI, x = age_levels, nx = nx, cond_age = 0)
  hor_li <- horiuchi(func = sd.cod.fun.65plus, pars1 = mx1, pars2 = mx2, N = N_HORIUCHI, x = age_levels, nx = nx, cond_age = 0)
  list(le = as.vector(hor_le), li = as.vector(hor_li))
}

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

#----- Object 3: LE65 and LI65+ levels

levels_dt <- cause_dt[cause == cause_levels[1], .(
  LE65 = life_expectancy_from_mx_65plus(mx_total[order(age)], x = age_levels, nx = nx, age = 0),
  LI65 = sd_from_mx_fun_65_plus(mx_total[order(age)], x = age_levels, nx = nx, age = 0)
), by = .(branch, year)]
levels_dt[, `:=`(geo_id = city_id, label = city_name, gcm = gcm_name)]
setcolorder(levels_dt, c("geo_id", "label", "gcm", "branch", "year", "LE65", "LI65"))
setorder(levels_dt, branch, year)

if (levels_only) {
  fwrite(levels_dt, levels_file)
  message("Levels-only run: saved LE/LI levels to ", levels_file)
  quit(save = "no")
}

#----- Object 2: year-on-year decomposition within each branch

decomp_le <- list(); decomp_li <- list(); step_checks <- list()
for (b in if (decomp_annual) branch_levels else character(0)) {
  years <- sort(unique(cause_dt[branch == b]$year))
  for (i in seq_len(length(years) - 1L)) {
    y0 <- years[i]; y1 <- years[i + 1L]
    h <- cached(sprintf("within_%s_%d_%d_N%d", b, y0, y1, N_HORIUCHI),
      horiuchi_pair(cause_vector(cause_dt[branch == b & year == y0]), cause_vector(cause_dt[branch == b & year == y1])))
    le0 <- levels_dt[branch == b & year == y0]; le1 <- levels_dt[branch == b & year == y1]
    base <- data.table(geo_id = city_id, label = city_name, gcm = gcm_name, branch = b, year_from = y0, year_to = y1)
    decomp_le[[length(decomp_le) + 1L]] <- cbind(base, grid_dt(h$le))
    decomp_li[[length(decomp_li) + 1L]] <- cbind(base, grid_dt(h$li))
    step_checks[[length(step_checks) + 1L]] <- data.table(
      branch = b, year_from = y0, year_to = y1,
      le_change = le1$LE65 - le0$LE65, li_change = le1$LI65 - le0$LI65,
      le_closure_error = sum(h$le) - (le1$LE65 - le0$LE65),
      li_closure_error = sum(h$li) - (le1$LI65 - le0$LI65)
    )
  }
}
decomp_le <- rbindlist(decomp_le); decomp_li <- rbindlist(decomp_li); step_checks <- rbindlist(step_checks)
if (!decomp_annual) step_checks <- data.table(branch = character(0), year_from = integer(0), year_to = integer(0),
  le_change = numeric(0), li_change = numeric(0), le_closure_error = numeric(0), li_closure_error = numeric(0))

# Simon's test (11 Sep): life-table change over the whole period must equal the
# sum of all age x cause contributions
period_checks <- rbindlist(lapply(if (decomp_annual) branch_levels else character(0), function(b) {
  l <- levels_dt[branch == b]
  data.table(branch = b,
    le_error = decomp_le[branch == b, sum(contribution)] - (l[year == max(year)]$LE65 - l[year == min(year)]$LE65),
    li_error = decomp_li[branch == b, sum(contribution)] - (l[year == max(year)]$LI65 - l[year == min(year)]$LI65))
}))

#----- Between-branch decomposition: with-CC vs without-CC, by 5-year period

cause_dt[, period := (year %/% period_len) * period_len]
period_mx <- cause_dt[, .(mx_cause = mean(mx_cause)), by = .(branch, period, age, cause)]
between <- list(); between_checks <- list()
for (p in sort(unique(period_mx$period))) {
  mx_wo <- cause_vector(period_mx[branch == "without_cc" & period == p])
  mx_w <- cause_vector(period_mx[branch == "with_cc" & period == p])
  h <- cached(sprintf("between_%d_N%d", p, N_HORIUCHI), horiuchi_pair(mx_wo, mx_w))
  f_le <- function(v) life_expectancy_cod(v, x = age_levels, nx = nx)
  f_li <- function(v) sd.cod.fun.65plus(v, x = age_levels, nx = nx)
  g <- grid_dt(h$le); setnames(g, "contribution", "le_contribution"); g[, li_contribution := h$li]
  between[[length(between) + 1L]] <- cbind(data.table(geo_id = city_id, label = city_name, gcm = gcm_name,
    period = sprintf("%d-%d", p, p + period_len - 1L)), g)
  between_checks[[length(between_checks) + 1L]] <- data.table(period = p,
    LE65_without = f_le(mx_wo), LE65_with = f_le(mx_w), LI65_without = f_li(mx_wo), LI65_with = f_li(mx_w),
    le_closure_error = sum(h$le) - (f_le(mx_w) - f_le(mx_wo)),
    li_closure_error = sum(h$li) - (f_li(mx_w) - f_li(mx_wo)),
    max_abs_rest = max(abs(g[cause == "rest", c(le_contribution, li_contribution)])))
}
between <- rbindlist(between); between_checks <- rbindlist(between_checks)

#----- Within-branch decomposition between consecutive 5-year-period means
# Summing year-on-year contributions over a block telescopes to the change
# between the block's first and last single years, which is dominated by those
# years' weather. Lloyd et al. (2024) decomposed changes in multi-year average
# ANs; here consecutive 5-year-period mean schedules are used, so block sums
# reflect the change between period means (used by 05_figures.R, Fig 2).

periods_all <- sort(unique(period_mx$period))
f_le_p <- function(v) life_expectancy_cod(v, x = age_levels, nx = nx)
f_li_p <- function(v) sd.cod.fun.65plus(v, x = age_levels, nx = nx)
within_p <- list(); within_p_checks <- list()
for (b in branch_levels) {
  for (i in seq_len(length(periods_all) - 1L)) {
    p0 <- periods_all[i]; p1 <- periods_all[i + 1L]
    m0 <- cause_vector(period_mx[branch == b & period == p0]); m1 <- cause_vector(period_mx[branch == b & period == p1])
    h <- cached(sprintf("withinperiod_%s_%d_%d_N%d", b, p0, p1, N_HORIUCHI), horiuchi_pair(m0, m1))
    g <- grid_dt(h$le); setnames(g, "contribution", "le_contribution"); g[, li_contribution := h$li]
    within_p[[length(within_p) + 1L]] <- cbind(data.table(geo_id = city_id, label = city_name, gcm = gcm_name, branch = b,
      period_from = sprintf("%d-%d", p0, p0 + period_len - 1L), period_to = sprintf("%d-%d", p1, p1 + period_len - 1L)), g)
    within_p_checks[[length(within_p_checks) + 1L]] <- data.table(branch = b, period_from = p0,
      le_closure_error = sum(h$le) - (f_le_p(m1) - f_le_p(m0)), li_closure_error = sum(h$li) - (f_li_p(m1) - f_li_p(m0)))
  }
}
within_p <- rbindlist(within_p); within_p_checks <- rbindlist(within_p_checks)

#----- Invariant checks

checks <- data.table(
  check_name = c("le_step_closure", "li_step_closure", "le_whole_period_closure", "li_whole_period_closure",
                 "between_branch_closure", "between_branch_rest_zero", "levels_complete", "within_period_closure"),
  status = c(
    if (!decomp_annual) "SKIPPED" else if (max(abs(step_checks$le_closure_error)) <= closure_tol) "PASS" else "FAIL",
    if (!decomp_annual) "SKIPPED" else if (max(abs(step_checks$li_closure_error)) <= closure_tol) "PASS" else "FAIL",
    if (!decomp_annual) "SKIPPED" else if (max(abs(period_checks$le_error)) <= closure_tol) "PASS" else "FAIL",
    if (!decomp_annual) "SKIPPED" else if (max(abs(period_checks$li_error)) <= closure_tol) "PASS" else "FAIL",
    if (max(abs(c(between_checks$le_closure_error, between_checks$li_closure_error))) <= closure_tol) "PASS" else "FAIL",
    if (max(between_checks$max_abs_rest) <= 1e-12) "PASS" else "FAIL",
    if (nrow(levels_dt) == length(branch_levels) * uniqueN(cause_dt$year) && all(is.finite(c(levels_dt$LE65, levels_dt$LI65)))) "PASS" else "FAIL",
    if (max(abs(c(within_p_checks$le_closure_error, within_p_checks$li_closure_error))) <= closure_tol) "PASS" else "FAIL"
  ),
  value = c(
    if (!decomp_annual) "annual decomposition not run" else sprintf("max_abs_error=%0.3e", max(abs(step_checks$le_closure_error))),
    if (!decomp_annual) "annual decomposition not run" else sprintf("max_abs_error=%0.3e", max(abs(step_checks$li_closure_error))),
    if (!decomp_annual) "annual decomposition not run" else sprintf("max_abs_error=%0.3e", max(abs(period_checks$le_error))),
    if (!decomp_annual) "annual decomposition not run" else sprintf("max_abs_error=%0.3e", max(abs(period_checks$li_error))),
    sprintf("max_abs_error=%0.3e", max(abs(c(between_checks$le_closure_error, between_checks$li_closure_error)))),
    sprintf("max_abs_rest=%0.3e", max(between_checks$max_abs_rest)),
    sprintf("%d rows", nrow(levels_dt)),
    sprintf("max_abs_error=%0.3e", max(abs(c(within_p_checks$le_closure_error, within_p_checks$li_closure_error))))
  ),
  threshold = c(rep(sprintf("<= %g", closure_tol), 5), "<= 1e-12 (rest identical in both branches)", "branch x year, finite", sprintf("<= %g", closure_tol))
)

fwrite(checks, checks_file)
fwrite(rbind(step_checks, fill = TRUE), steps_file)
failed <- checks[status == "FAIL"]
if (nrow(failed)) {
  fwrite(failed, failures_file)
  stop(sprintf("04_le_li_decomposition.R failed %d invariant(s); see %s", nrow(failed), failures_file), call. = FALSE)
} else {
  if (file.exists(failures_file)) file.remove(failures_file)
  invisible(file.create(failures_file))
}

# fwrite keeps 15 significant digits, so small single-age contributions (~1e-4)
# are not rounded to zero (10 Sep meeting, "small numbers")
fwrite(levels_dt, levels_file)
if (decomp_annual) {
  fwrite(decomp_le, le_file)
  fwrite(decomp_li, li_file)
}
fwrite(between, between_file)
fwrite(within_p, within_period_file)

if (decomp_annual) {
  plot_checks <- melt(step_checks[, .(branch, year_from, le_closure_error, li_closure_error)],
    id.vars = c("branch", "year_from"), variable.name = "measure", value.name = "error")
  p <- ggplot(plot_checks, aes(x = year_from, y = error, color = measure)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_line(linewidth = 0.7) +
    facet_wrap(~branch, scales = "free_y") +
    labs(title = sprintf("%s LE/LI decomposition closure diagnostics", city_name),
         subtitle = sprintf("Horiuchi N = %d; closure tolerance %g", N_HORIUCHI, closure_tol),
         x = "Year from", y = "Closure error") +
    theme_minimal(base_size = 11)
  ggsave(fig_file, p, width = 11, height = 6, dpi = 160)
}

message("Saved LE/LI levels to ", levels_file)
message("Saved LE decomposition to ", le_file)
message("Saved LI decomposition to ", li_file)
message("Saved between-branch decomposition to ", between_file)
message("Saved within-branch period decomposition to ", within_period_file)
message("Saved checks to ", checks_file)
