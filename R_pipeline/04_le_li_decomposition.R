#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality / life expectancy pipeline -- Madrid pilot
#
# R Pipeline Step 04: Life tables, LE65 / LI65+ levels and decompositions
#   Life-table, LE and SD functions are those of Lloyd et al. (2024) Code_1.R /
#   Code_2.R (modified from Aburto et al. 2022). Outputs (the three data objects
#   agreed at the 10 Sep 2026 meeting, Objects 2 and 3 written here):
#     04_le_li_levels.csv                LE65 and LI65+ by branch x year   (Object 3)
#     04_le_decomposition.csv            Horiuchi contributions to year-on-year
#     04_li_decomposition.csv              change, by branch x age x cause (Object 2)
#     04_between_branch_decomposition.csv  with-CC minus without-CC difference,
#                                        by 5-year period x age x cause
#   Horiuchi steps are cached in results/phase1_madrid/04_cache/<md5 of the
#   master table>/ so an interrupted run resumes where it stopped.
#
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(DemoDecomp)
  library(ggplot2)
})

message("\n[04] Running Madrid LE/LI decomposition (N = 400)...")

out_dir <- "results/phase1_madrid"
check_dir <- "results/checks"
fig_dir <- "results/figures"
for (d in c(out_dir, check_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
# The cache is keyed on the md5 of the master table, so any upstream change
# starts a fresh cache instead of reusing stale contributions.
master_file <- file.path(out_dir, "03_master_table.csv")
cache_dir <- file.path(out_dir, "04_cache", unname(tools::md5sum(master_file)))
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

le_file <- file.path(out_dir, "04_le_decomposition.csv")
li_file <- file.path(out_dir, "04_li_decomposition.csv")
levels_file <- file.path(out_dir, "04_le_li_levels.csv")
between_file <- file.path(out_dir, "04_between_branch_decomposition.csv")
checks_file <- file.path(check_dir, "04_le_li_decomposition_checks.csv")
steps_file <- file.path(check_dir, "04_le_li_decomposition_steps.csv")
failures_file <- file.path(check_dir, "04_le_li_decomposition_failures.csv")
fig_file <- file.path(fig_dir, "04_le_li_decomposition_diagnostic.png")

city_id <- "ES001C"
city_name <- "Madrid"
gcm_target <- "GFDL_ESM4"
branch_levels <- c("with_cc", "without_cc")
cause_levels <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat", "rest")
age_levels <- 65:100
nx <- c(rep(1, 100 - 65), Inf)
N_HORIUCHI <- 400L
# Between-branch decomposition is run on mean mortality rates over the 5-year
# periods of the demographic projections (annual with - without differences
# mix the climate signal with single-year weather).
period_len <- 5L
closure_tol <- 1e-6

master <- fread(master_file)
master <- master[geo_id == city_id & gcm == gcm_target]

if (!nrow(master)) stop("Master table for Madrid is missing; run 03_master_table.R first.", call. = FALSE)

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
levels_dt[, `:=`(geo_id = city_id, label = city_name, gcm = gcm_target)]
setcolorder(levels_dt, c("geo_id", "label", "gcm", "branch", "year", "LE65", "LI65"))
setorder(levels_dt, branch, year)

#----- Object 2: year-on-year decomposition within each branch

decomp_le <- list(); decomp_li <- list(); step_checks <- list()
for (b in branch_levels) {
  years <- sort(unique(cause_dt[branch == b]$year))
  for (i in seq_len(length(years) - 1L)) {
    y0 <- years[i]; y1 <- years[i + 1L]
    h <- cached(sprintf("within_%s_%d_%d_N%d", b, y0, y1, N_HORIUCHI),
      horiuchi_pair(cause_vector(cause_dt[branch == b & year == y0]), cause_vector(cause_dt[branch == b & year == y1])))
    le0 <- levels_dt[branch == b & year == y0]; le1 <- levels_dt[branch == b & year == y1]
    base <- data.table(geo_id = city_id, label = city_name, gcm = gcm_target, branch = b, year_from = y0, year_to = y1)
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

# Simon's test (11 Sep): life-table change over the whole period must equal the
# sum of all age x cause contributions
period_checks <- rbindlist(lapply(branch_levels, function(b) {
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
  between[[length(between) + 1L]] <- cbind(data.table(geo_id = city_id, label = city_name, gcm = gcm_target,
    period = sprintf("%d-%d", p, p + period_len - 1L)), g)
  between_checks[[length(between_checks) + 1L]] <- data.table(period = p,
    LE65_without = f_le(mx_wo), LE65_with = f_le(mx_w), LI65_without = f_li(mx_wo), LI65_with = f_li(mx_w),
    le_closure_error = sum(h$le) - (f_le(mx_w) - f_le(mx_wo)),
    li_closure_error = sum(h$li) - (f_li(mx_w) - f_li(mx_wo)),
    max_abs_rest = max(abs(g[cause == "rest", c(le_contribution, li_contribution)])))
}
between <- rbindlist(between); between_checks <- rbindlist(between_checks)

#----- Invariant checks

checks <- data.table(
  check_name = c("le_step_closure", "li_step_closure", "le_whole_period_closure", "li_whole_period_closure",
                 "between_branch_closure", "between_branch_rest_zero", "levels_complete"),
  status = c(
    if (max(abs(step_checks$le_closure_error)) <= closure_tol) "PASS" else "FAIL",
    if (max(abs(step_checks$li_closure_error)) <= closure_tol) "PASS" else "FAIL",
    if (max(abs(period_checks$le_error)) <= closure_tol) "PASS" else "FAIL",
    if (max(abs(period_checks$li_error)) <= closure_tol) "PASS" else "FAIL",
    if (max(abs(c(between_checks$le_closure_error, between_checks$li_closure_error))) <= closure_tol) "PASS" else "FAIL",
    if (max(between_checks$max_abs_rest) <= 1e-12) "PASS" else "FAIL",
    if (nrow(levels_dt) == length(branch_levels) * uniqueN(cause_dt$year) && all(is.finite(c(levels_dt$LE65, levels_dt$LI65)))) "PASS" else "FAIL"
  ),
  value = c(
    sprintf("max_abs_error=%0.3e", max(abs(step_checks$le_closure_error))),
    sprintf("max_abs_error=%0.3e", max(abs(step_checks$li_closure_error))),
    sprintf("max_abs_error=%0.3e", max(abs(period_checks$le_error))),
    sprintf("max_abs_error=%0.3e", max(abs(period_checks$li_error))),
    sprintf("max_abs_error=%0.3e", max(abs(c(between_checks$le_closure_error, between_checks$li_closure_error)))),
    sprintf("max_abs_rest=%0.3e", max(between_checks$max_abs_rest)),
    sprintf("%d rows", nrow(levels_dt))
  ),
  threshold = c(rep(sprintf("<= %g", closure_tol), 5), "<= 1e-12 (rest identical in both branches)", "branch x year, finite")
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
fwrite(decomp_le, le_file)
fwrite(decomp_li, li_file)
fwrite(between, between_file)

plot_checks <- melt(step_checks[, .(branch, year_from, le_closure_error, li_closure_error)],
  id.vars = c("branch", "year_from"), variable.name = "measure", value.name = "error")
p <- ggplot(plot_checks, aes(x = year_from, y = error, color = measure)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~branch, scales = "free_y") +
  labs(title = "Madrid LE/LI decomposition closure diagnostics",
       subtitle = sprintf("Horiuchi N = %d; closure tolerance %g", N_HORIUCHI, closure_tol),
       x = "Year from", y = "Closure error") +
  theme_minimal(base_size = 11)
ggsave(fig_file, p, width = 11, height = 6, dpi = 160)

message("Saved LE/LI levels to ", levels_file)
message("Saved LE decomposition to ", le_file)
message("Saved LI decomposition to ", li_file)
message("Saved between-branch decomposition to ", between_file)
message("Saved checks to ", checks_file)
