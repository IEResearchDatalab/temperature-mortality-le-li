#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DemoDecomp)
  library(parallel)
  library(ggplot2)
})

message("\n[04] Running Madrid LE/LI decomposition (N = 400)...")

out_dir <- "results/phase1_madrid"
check_dir <- "results/checks"
fig_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

le_file <- file.path(out_dir, "04_le_decomposition.csv")
li_file <- file.path(out_dir, "04_li_decomposition.csv")
checks_file <- file.path(check_dir, "04_le_li_decomposition_checks.csv")
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

master <- fread(file.path(out_dir, "03_master_table.csv"))
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
cause_dt[, mx_total := total_deaths / pop]

if (any(!cause_dt$cause %in% cause_levels)) {
  stop("Unexpected cause labels in decomposition input.", call. = FALSE)
}
if (anyDuplicated(cause_dt, by = c("branch", "year", "age", "cause"))) {
  stop("Decomposition input has duplicate branch/year/age/cause keys.", call. = FALSE)
}
if (nrow(cause_dt) != length(branch_levels) * length(sort(unique(master$year))) * length(age_levels) * length(cause_levels)) {
  stop("Decomposition input grid is incomplete or oversized.", call. = FALSE)
}

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

build_year_matrix <- function(dt, year_value) {
  x <- sort(unique(dt$age))
  if (!all(x == age_levels)) stop(sprintf("Age grid mismatch for year %s.", year_value), call. = FALSE)
  ordered <- dt[order(match(cause, cause_levels), age)]
  list(x = x, mx = ordered$mx_cause, mx_total = ordered$mx_total[match(x, ordered$age)])
}

decompose_branch <- function(branch_name) {
  branch_dt <- cause_dt[branch == branch_name]
  years <- sort(unique(branch_dt$year))
  results_le <- list()
  results_li <- list()
  checks <- list()

  for (i in seq_len(length(years) - 1L)) {
    y0 <- years[i]
    y1 <- years[i + 1L]
    dt0 <- branch_dt[year == y0]
    dt1 <- branch_dt[year == y1]

    mat0 <- build_year_matrix(dt0, y0)
    mat1 <- build_year_matrix(dt1, y1)

    le0 <- life_expectancy_from_mx_65plus(mat0$mx_total, x = age_levels, nx = nx, age = 0)
    le1 <- life_expectancy_from_mx_65plus(mat1$mx_total, x = age_levels, nx = nx, age = 0)
    li0 <- sd_from_mx_fun_65_plus(mat0$mx_total, x = age_levels, nx = nx, age = 0)
    li1 <- sd_from_mx_fun_65_plus(mat1$mx_total, x = age_levels, nx = nx, age = 0)

    hor_le <- horiuchi(
      func = life_expectancy_cod,
      pars1 = mat0$mx,
      pars2 = mat1$mx,
      N = N_HORIUCHI,
      x = age_levels,
      nx = nx,
      cond_age = 0
    )
    hor_li <- horiuchi(
      func = sd.cod.fun.65plus,
      pars1 = mat0$mx,
      pars2 = mat1$mx,
      N = N_HORIUCHI,
      x = age_levels,
      nx = nx,
      cond_age = 0
    )

    dim(hor_le) <- c(length(age_levels), length(cause_levels))
    dim(hor_li) <- c(length(age_levels), length(cause_levels))

    le_df <- data.table(
      geo_id = city_id,
      label = city_name,
      gcm = gcm_target,
      branch = branch_name,
      year_from = y0,
      year_to = y1,
      age = rep(age_levels, length(cause_levels)),
      cause = rep(cause_levels, each = length(age_levels)),
      contribution = as.vector(hor_le)
    )
    li_df <- copy(le_df)
    li_df[, contribution := as.vector(hor_li)]

    results_le[[length(results_le) + 1L]] <- le_df
    results_li[[length(results_li) + 1L]] <- li_df

    checks[[length(checks) + 1L]] <- data.table(
      geo_id = city_id,
      branch = branch_name,
      year_from = y0,
      year_to = y1,
      le_from = le0,
      le_to = le1,
      li_from = li0,
      li_to = li1,
      le_change = le1 - le0,
      li_change = li1 - li0,
      le_contrib_sum = sum(le_df$contribution),
      li_contrib_sum = sum(li_df$contribution),
      le_closure_error = sum(le_df$contribution) - (le1 - le0),
      li_closure_error = sum(li_df$contribution) - (li1 - li0),
      le_sign_plausible = abs(le1 - le0) < 1e-12 || any(sign(le_df$contribution[abs(le_df$contribution) > 0]) == sign(le1 - le0)),
      li_sign_plausible = abs(li1 - li0) < 1e-12 || any(sign(li_df$contribution[abs(li_df$contribution) > 0]) == sign(li1 - li0))
    )
  }

  list(
    le = rbindlist(results_le, use.names = TRUE),
    li = rbindlist(results_li, use.names = TRUE),
    checks = rbindlist(checks, use.names = TRUE)
  )
}

branch_results <- lapply(branch_levels, decompose_branch)

decomp_le <- rbindlist(lapply(branch_results, `[[`, "le"), use.names = TRUE)
decomp_li <- rbindlist(lapply(branch_results, `[[`, "li"), use.names = TRUE)
checks <- rbindlist(lapply(branch_results, `[[`, "checks"), use.names = TRUE)

setorder(decomp_le, branch, year_from, cause, age)
setorder(decomp_li, branch, year_from, cause, age)
setorder(checks, branch, year_from)

checks[, `:=`(
  le_closure_pass = abs(le_closure_error) <= 1e-6,
  li_closure_pass = abs(li_closure_error) <= 1e-6,
  sign_plausible = le_sign_plausible & li_sign_plausible
)]

summary_checks <- data.table(
  check_name = c("le_closure", "li_closure", "sign_plausibility", "branch_consistency"),
  status = c(
    if (all(checks$le_closure_pass)) "PASS" else "FAIL",
    if (all(checks$li_closure_pass)) "PASS" else "FAIL",
    if (all(checks$sign_plausible)) "PASS" else "FAIL",
    if (length(unique(c(unique(decomp_le$branch), unique(decomp_li$branch)))) == length(branch_levels)) "PASS" else "FAIL"
  ),
  value = c(
    sprintf("max_abs_error=%0.3e", max(abs(checks$le_closure_error))),
    sprintf("max_abs_error=%0.3e", max(abs(checks$li_closure_error))),
    sprintf("plausible_rows=%d/%d", sum(checks$sign_plausible), nrow(checks)),
    paste(branch_levels, collapse = ",")
  ),
  threshold = c(
    "<= 1e-6",
    "<= 1e-6",
    "all rows plausible",
    "both branches present"
  )
)

failures <- data.table()
if (any(summary_checks$status == "FAIL")) {
  failures <- rbindlist(list(
    if (summary_checks$status[summary_checks$check_name == "le_closure"] == "FAIL") {
      checks[le_closure_pass == FALSE, .(
        geo_id = city_id,
        branch,
        year_from,
        year_to,
        failing_check = "le_closure",
        observed_value = le_closure_error,
        expected_bound = "<= 1e-6"
      )]
    } else NULL,
    if (summary_checks$status[summary_checks$check_name == "li_closure"] == "FAIL") {
      checks[li_closure_pass == FALSE, .(
        geo_id = city_id,
        branch,
        year_from,
        year_to,
        failing_check = "li_closure",
        observed_value = li_closure_error,
        expected_bound = "<= 1e-6"
      )]
    } else NULL,
    if (summary_checks$status[summary_checks$check_name == "sign_plausibility"] == "FAIL") {
      checks[sign_plausible == FALSE, .(
        geo_id = city_id,
        branch,
        year_from,
        year_to,
        failing_check = "sign_plausibility",
        observed_value = sprintf("le_plausible=%s; li_plausible=%s", le_sign_plausible, li_sign_plausible),
        expected_bound = "at least one same-signed contribution exists for each measure"
      )]
    } else NULL
  ), fill = TRUE)
}

fwrite(decomp_le, le_file)
fwrite(decomp_li, li_file)
fwrite(summary_checks, checks_file)
if (nrow(failures)) {
  fwrite(failures, failures_file)
} else {
  if (file.exists(failures_file)) file.remove(failures_file)
  invisible(file.create(failures_file))
}

plot_checks <- melt(
  checks[, .(branch, year_from, le_closure_error, li_closure_error)],
  id.vars = c("branch", "year_from"),
  variable.name = "measure",
  value.name = "error"
)

p <- ggplot(plot_checks, aes(x = year_from, y = error, color = measure)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~branch, scales = "free_y") +
  labs(
    title = "Madrid LE/LI decomposition closure diagnostics",
    subtitle = sprintf("Horiuchi N = %d; closure tolerance 1e-6", N_HORIUCHI),
    x = "Year from",
    y = "Closure error"
  ) +
  theme_minimal(base_size = 11)

ggsave(fig_file, p, width = 11, height = 6, dpi = 160)

if (nrow(failures)) {
  stop(sprintf("04_le_li_decomposition.R failed %d invariant(s); see %s", nrow(failures), failures_file), call. = FALSE)
}

message("Saved LE decomposition to ", le_file)
message("Saved LI decomposition to ", li_file)
message("Saved checks to ", checks_file)
message("Saved diagnostic figure to ", fig_file)
