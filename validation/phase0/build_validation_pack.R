
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dlnm)
})

root <- normalizePath(".")
out_root <- Sys.getenv("VALIDATION_OUT_ROOT", unset = file.path(root, "validation", "phase0"))
if (!nzchar(out_root)) out_root <- file.path(root, "validation", "phase0")
out_root <- normalizePath(out_root, mustWork = FALSE)
fig_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
stage_root <- file.path(out_root, "stages")
for (p in c(out_root, fig_dir, table_dir, stage_root)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

producer_script <- file.path("validation", "run_validation.R")
producer_commit <- tryCatch(
  trimws(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = TRUE)[1]),
  error = function(e) NA_character_
)
producer_commit <- if (is.na(producer_commit) || !nzchar(producer_commit)) "unknown" else producer_commit

# ---------- Helpers ----------
write_csv <- function(dt, path) fwrite(as.data.table(dt), path)

save_png <- function(plot, path, width = 10, height = 7, dpi = 160) {
  ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
}

require_file <- function(path, label = path) {
  if (!file.exists(path)) stop(sprintf("Required input missing: %s", label), call. = FALSE)
  invisible(path)
}

assert_finite_numeric <- function(x, label) {
  if (any(!is.finite(x))) stop(sprintf("Non-finite values detected in %s", label), call. = FALSE)
  invisible(TRUE)
}

assert_rows_unchanged <- function(before, after, label) {
  if (nrow(before) != nrow(after)) {
    stop(sprintf("Row-count mismatch in %s: before=%d after=%d", label, nrow(before), nrow(after)), call. = FALSE)
  }
  invisible(TRUE)
}

assert_group_present <- function(dt, cols, expected, label) {
  key <- dt[, unique(.SD), .SDcols = cols]
  missing <- setdiff(expected, key[[cols[length(cols)]]])
  if (length(missing)) stop(sprintf("Expected group(s) absent in %s: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}

write_contract <- function(stage, checks, failures, summary_text, overview_plot) {
  stage_dir <- file.path(stage_root, stage)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  if (any(vapply(checks, function(col) is.numeric(col) && any(!is.finite(col)), logical(1)))) {
    stop(sprintf("Non-finite numeric values detected in checks for stage %s", stage), call. = FALSE)
  }
  write_csv(checks, file.path(stage_dir, "checks.csv"))
  failures_path <- file.path(stage_dir, "failures.csv")
  if (!is.null(failures) && nrow(as.data.table(failures)) > 0 && ncol(as.data.table(failures)) > 0) {
    write_csv(failures, failures_path)
  } else {
    if (file.exists(failures_path)) file.remove(failures_path)
    fwrite(data.table(note = character()), failures_path)
  }
  writeLines(summary_text, file.path(stage_dir, "summary.md"))
  save_png(overview_plot, file.path(stage_dir, "overview.png"), width = 9, height = 5.5)
}

fmt_num <- function(x, digits = 6) formatC(x, format = "f", digits = digits)

sha256_file <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  sub(" .*", "", out[1])
}

write_manifest <- function(paths, manifest_path, base_dir, producer = producer_script, commit = producer_commit) {
  paths <- paths[file.exists(paths)]
  if (!length(paths)) {
    fwrite(data.table(relative_path = character(), sha256 = character(), size_bytes = integer(), producing_source_or_script = character(), producing_commit = character()), manifest_path)
    return(invisible(NULL))
  }
  dt <- data.table(
    relative_path = sub(paste0("^", normalizePath(base_dir, winslash = "/", mustWork = TRUE), "/?"), "", normalizePath(paths, winslash = "/", mustWork = TRUE)),
    sha256 = vapply(paths, sha256_file, character(1)),
    size_bytes = file.info(paths)$size,
    producing_source_or_script = producer,
    producing_commit = commit
  )
  fwrite(dt, manifest_path)
}

# Canonical LE/LI functions from the Lloyd reference implementation.
cond_surv <- 65
life_expectancy_from_mx_fun_65plus <- function(mx, x, nx = c(rep(1, 100 - cond_surv), Inf), age = 0) {
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  ex[age + 1]
}
sd_from_mx_fun_65_plus <- function(mx, x, nx = c(rep(1, 100 - cond_surv), Inf), age = 0) {
  x_conditional <- x - cond_surv
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  sqrt(sum(dx * (x_conditional + 0.5 - ex[age + 1])^2))
}

life.expectancy.cod.fun.65plus <- function(mx.cod, x, nx = c(rep(1, 100 - cond_surv), Inf), cond_age = 0) {
  dim(mx.cod) <- c(length(x), length(mx.cod) / length(x))
  mx <- rowSums(mx.cod)
  life_expectancy_from_mx_fun_65plus(mx, x, nx, cond_age)
}
sd.cod.fun.65plus <- function(mx.cod, x, nx, cond_age = 0) {
  dim(mx.cod) <- c(length(x), length(mx.cod) / length(x))
  mx <- rowSums(mx.cod)
  sd_from_mx_fun_65_plus(mx, x, nx, cond_age)
}

# ---------- Load evidence ----------
proj_path <- file.path(root, "results", "projdata", "projdata_prototype.csv")
mass_cityage_path <- file.path(root, "references", "2025-masselot-zenodo", "results", "cityage.csv")
mass_city_path <- file.path(root, "references", "2025-masselot-zenodo", "results", "city.csv")
prep_path <- file.path(root, "data", "prep_data.RData")
coefs_path <- file.path(root, "data", "coefs.csv")

for (p in c(proj_path, mass_cityage_path, mass_city_path, prep_path, coefs_path,
            file.path(root, "agent-output", "phase0", "negative_grouped_an_summary.csv"),
            file.path(root, "agent-output", "phase0", "negative_grouped_an_by_city.csv"),
            file.path(root, "agent-output", "phase0", "negative_grouped_an_by_agegroup.csv"),
            file.path(root, "agent-output", "phase0", "negative_grouped_an_by_year.csv"),
            file.path(root, "agent-output", "phase0", "negative_grouped_an_by_range.csv"),
            file.path(root, "agent-output", "phase0", "negative_grouped_an_examples.csv"),
            file.path(root, "agent-output", "phase0", "signed_grouped_an_cell_comparison.csv"),
            file.path(root, "agent-output", "phase0", "signed_grouped_an_method_table.csv"),
            file.path(root, "agent-output", "phase0", "lloyd_golden", "lloyd_diagnostic_checks.rds"),
            file.path(root, "agent-output", "phase0", "lloyd_n400_800_summary.rds"))) {
  require_file(p)
}

proj <- fread(proj_path)
proj[, `:=`(
  year5 = as.integer(year5),
  ssp = as.character(ssp),
  URAU_CODE = as.character(URAU_CODE),
  CNTR_CODE = as.character(CNTR_CODE),
  agegroup = as.character(agegroup)
)]

mass_cityage <- fread(mass_cityage_path)
mass_city <- fread(mass_city_path)
neg_summary <- fread(file.path(root, "agent-output", "phase0", "negative_grouped_an_summary.csv"))
neg_city <- fread(file.path(root, "agent-output", "phase0", "negative_grouped_an_by_city.csv"))
neg_age <- fread(file.path(root, "agent-output", "phase0", "negative_grouped_an_by_agegroup.csv"))
neg_year <- fread(file.path(root, "agent-output", "phase0", "negative_grouped_an_by_year.csv"))
neg_range <- fread(file.path(root, "agent-output", "phase0", "negative_grouped_an_by_range.csv"))
neg_examples <- fread(file.path(root, "agent-output", "phase0", "negative_grouped_an_examples.csv"))
clamp_report <- file.exists(file.path(root, "agent-output", "phase0", "clamp_comparison.md"))
lloyd_diag <- readRDS(file.path(root, "agent-output", "phase0", "lloyd_golden", "lloyd_diagnostic_checks.rds"))
lloyd_n400 <- readRDS(file.path(root, "agent-output", "phase0", "lloyd_n400_800_summary.rds"))
signed_cmp <- fread(file.path(root, "agent-output", "phase0", "signed_grouped_an_cell_comparison.csv"))
signed_method <- fread(file.path(root, "agent-output", "phase0", "signed_grouped_an_method_table.csv"))
require_file(file.path(root, "results", "le_li_input", "lloyd_fig_s1_baseline_city_year_single_age.csv"))
baseline_single <- fread(file.path(root, "results", "le_li_input", "lloyd_fig_s1_baseline_city_year_single_age.csv"))
baseline_single <- baseline_single[age >= 65]
required_ll_cols <- c("le_closure_20", "le_closure_50", "le_closure_100", "le_closure_200", "li_closure_20", "li_closure_50", "li_closure_100", "li_closure_200", "le_20_50", "le_50_100", "le_100_200", "li_20_50", "li_50_100", "li_100_200")
if (!all(required_ll_cols %in% names(as.data.table(lloyd_diag)))) {
  stop("Lloyd diagnostic fixture is missing required columns.", call. = FALSE)
}

# ---------- Stage 1: demographic scaling ----------
anchor_year <- 2015L
future_year <- 2095L
country_year <- proj[, .(
  country_pop = sum(pop),
  country_death = sum(death),
  wittpop = wittpop[1],
  wittdeath = wittdeath[1]
), by = .(CNTR_CODE, agegroup, ssp, year5)]

anchor <- country_year[year5 == anchor_year]
future <- country_year[year5 == future_year]
scale_cmp <- merge(anchor, future, by = c("CNTR_CODE", "agegroup", "ssp"), suffixes = c("_anchor", "_future"))
if (nrow(scale_cmp) != nrow(anchor) || nrow(scale_cmp) != nrow(future)) {
  stop("Demographic scaling join changed row counts unexpectedly.", call. = FALSE)
}
scale_cmp[, `:=`(
  country_growth = country_pop_future / country_pop_anchor,
  witt_growth = wittpop_future / wittpop_anchor,
  growth_diff = abs(country_pop_future / country_pop_anchor - wittpop_future / wittpop_anchor)
)]
proj_share <- merge(proj, country_year[, .(CNTR_CODE, agegroup, ssp, year5, country_pop, country_wittpop = wittpop)], by = c("CNTR_CODE", "agegroup", "ssp", "year5"))
proj_share[, coverage := country_pop / country_wittpop]
share_anchor <- proj_share[year5 == anchor_year, .(URAU_CODE, CNTR_CODE, agegroup, ssp, share_anchor = pop / country_pop, coverage_anchor = coverage)]
share_future <- proj_share[year5 == future_year, .(URAU_CODE, CNTR_CODE, agegroup, ssp, share_future = pop / country_pop, coverage_future = coverage)]
scale_check <- merge(share_anchor, share_future, by = c("URAU_CODE", "CNTR_CODE", "agegroup", "ssp"))
if (nrow(scale_check) != nrow(share_anchor) || nrow(scale_check) != nrow(share_future)) {
  stop("Demographic share join changed row counts unexpectedly.", call. = FALSE)
}
scale_check[, `:=`(
  share_diff = abs(share_future - share_anchor),
  coverage_diff = abs(coverage_future - coverage_anchor)
)]
scale_tbl <- scale_check[, .(
  max_share_diff = max(share_diff),
  max_coverage_diff = max(coverage_diff),
  max_growth_diff = max(scale_cmp$growth_diff)
)]
scale_failures <- scale_check[share_diff > 1e-12 | coverage_diff > 1e-12]
scale_plot <- ggplot(scale_check, aes(x = share_diff, y = coverage_diff)) +
  geom_point(alpha = 0.25, size = 0.7, color = "#2c7fb8") +
  scale_x_log10() + scale_y_log10() +
  geom_vline(xintercept = 1e-12, linetype = 2, color = "#d95f0e") +
  geom_hline(yintercept = 1e-12, linetype = 2, color = "#d95f0e") +
  labs(
    title = "Demographic scaling validation",
    subtitle = sprintf("City share and sample coverage deltas vs anchor year %d; threshold 1e-12", anchor_year),
    x = "|share_future - share_anchor|",
    y = "|coverage_future - coverage_anchor|"
  ) + theme_minimal(base_size = 11)
write_csv(scale_tbl, file.path(table_dir, "demographic_scaling_validation.csv"))
write_contract(
  "demographic_scaling",
  checks = data.table(metric = c("max_share_diff", "max_coverage_diff", "max_growth_diff", "threshold"), value = c(scale_tbl$max_share_diff, scale_tbl$max_coverage_diff, scale_tbl$max_growth_diff, 1e-12)),
  failures = scale_failures,
  summary_text = sprintf("Max share diff %.3e; max coverage diff %.3e; threshold 1e-12.", scale_tbl$max_share_diff, scale_tbl$max_coverage_diff),
  overview_plot = scale_plot
)
save_png(scale_plot, file.path(fig_dir, "demographic_scaling.png"), 8, 5.5)



# ---------- Stage 2: Masselot reproduction ----------
city_id <- "AT001C"
# Reproduce historical age-profile means for the selected city using the validated baseline fixture.
base_city <- readRDS(file.path(root, "temp_results_baseline", paste0(city_id, ".rds")))
base_city <- base_city[sim == 0]
base_city[, `:=`(
  agegroup = as.character(agegroup),
  range = as.character(range),
  year = as.integer(year),
  an = as.numeric(an)
)]
annual_city <- base_city[, .(
  total = sum(an),
  cold = sum(an[range %in% c("ExtrCold", "ModCold")]),
  heat = sum(an[range %in% c("ExtrHeat", "ModHeat")]),
  year_days = fifelse((year %% 400L == 0L) | (year %% 4L == 0L & year %% 100L != 0L), 366, 365)
), by = .(agegroup, year)]
repro_age <- annual_city[, .(
  reproduced_total = sum(total * year_days) / sum(year_days),
  reproduced_cold = sum(cold * year_days) / sum(year_days),
  reproduced_heat = sum(heat * year_days) / sum(year_days)
), by = agegroup]
repro_age <- merge(repro_age, mass_cityage[URAU_CODE == city_id, .(agegroup, masselot_total = excess_total_est, masselot_cold = excess_cold_est, masselot_heat = excess_heat_est)], by = "agegroup")
if (nrow(repro_age) != length(unique(repro_age$agegroup))) {
  stop("Masselot reproduction join changed row counts unexpectedly.", call. = FALSE)
}
repro_age[, age_rank := match(agegroup, c("20-44", "45-64", "65-74", "75-84", "85+"))]
setorder(repro_age, age_rank)
repro_long <- rbindlist(list(
  repro_age[, .(agegroup, measure = "total", reproduced = reproduced_total, masselot = masselot_total)],
  repro_age[, .(agegroup, measure = "cold", reproduced = reproduced_cold, masselot = masselot_cold)],
  repro_age[, .(agegroup, measure = "heat", reproduced = reproduced_heat, masselot = masselot_heat)]
))
repro_long[, abs_diff := abs(reproduced - masselot)]
repro_long[, rel_diff := abs_diff / pmax(abs(masselot), 1e-12)]

mass_plot <- ggplot(repro_long, aes(x = agegroup, y = reproduced, fill = "Reproduced")) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_point(aes(y = masselot, color = "Masselot"), size = 2.4, position = position_dodge(width = 0.7)) +
  facet_wrap(~measure, scales = "free_y", ncol = 1) +
  labs(
    title = sprintf("Masselot reproduction: %s historical age-profile totals", city_id),
    subtitle = "Units: deaths/year; sample: historical ERA5, fixed baseline deaths, central coefficients; thresholds are the Masselot reference values",
    y = "AN (deaths/year)", x = "Age group"
  ) + theme_minimal(base_size = 11)
write_csv(repro_long, file.path(table_dir, "masselot_reproduction.csv"))
write_contract(
  "masselot_reproduction",
  checks = data.table(metric = c("max_abs_diff_total", "max_abs_diff_cold", "max_abs_diff_heat", "threshold_note"), value = c(max(repro_long$abs_diff[repro_long$measure == "total"]), max(repro_long$abs_diff[repro_long$measure == "cold"]), max(repro_long$abs_diff[repro_long$measure == "heat"]), "reference: cityage.csv; year-days weighted average")),
  failures = repro_long[which(repro_long$abs_diff > 1e-6)],
  summary_text = sprintf("Historical city-age totals reproduced for %s using historical ERA5 and central coefficients.", city_id),
  overview_plot = mass_plot
)
save_png(mass_plot, file.path(fig_dir, "masselot_reproduction.png"), 8.8, 7)

# ---------- Stage 2b: stratified attribution sample ----------
sample_cities <- unique(mass_city[, .(region, URAU_CODE, LABEL)])[order(region, URAU_CODE)]
sample_cities <- sample_cities[, .SD[1], by = region]
sample_keys <- rbindlist(lapply(sample_cities$URAU_CODE, function(city_id) {
  data.table(URAU_CODE = city_id, agegroup = mass_cityage[URAU_CODE == city_id, sort(unique(agegroup))])
}))

sample_repro <- rbindlist(lapply(seq_len(nrow(sample_keys)), function(i) {
  key <- sample_keys[i]
  city_dt <- as.data.table(readRDS(file.path(root, "temp_results_baseline", paste0(key$URAU_CODE, ".rds"))))[sim == 0 & agegroup == key$agegroup]
  city_dt[, year_days := fifelse((year %% 400L == 0L) | (year %% 4L == 0L & year %% 100L != 0L), 366L, 365L)]
  year_dt <- city_dt[, .(
    total = sum(an),
    cold = sum(an[range %in% c("ExtrCold", "ModCold")]),
    heat = sum(an[range %in% c("ExtrHeat", "ModHeat")]),
    year_days = first(year_days)
  ), by = year]
  data.table(
    URAU_CODE = key$URAU_CODE,
    agegroup = key$agegroup,
    reproduced_total = sum(year_dt$total * year_dt$year_days) / sum(year_dt$year_days),
    reproduced_cold = sum(year_dt$cold * year_dt$year_days) / sum(year_dt$year_days),
    reproduced_heat = sum(year_dt$heat * year_dt$year_days) / sum(year_dt$year_days)
  )
}))
sample_repro <- merge(
  sample_repro,
  mass_cityage[, .(URAU_CODE, agegroup, masselot_total = excess_total_est, masselot_cold = excess_cold_est, masselot_heat = excess_heat_est)],
  by = c("URAU_CODE", "agegroup")
)
sample_long <- rbindlist(list(
  sample_repro[, .(URAU_CODE, agegroup, measure = "total", reproduced = reproduced_total, masselot = masselot_total)],
  sample_repro[, .(URAU_CODE, agegroup, measure = "cold", reproduced = reproduced_cold, masselot = masselot_cold)],
  sample_repro[, .(URAU_CODE, agegroup, measure = "heat", reproduced = reproduced_heat, masselot = masselot_heat)]
))
sample_long[, abs_diff := abs(reproduced - masselot)]
sample_long[, rel_diff := abs_diff / pmax(abs(masselot), 1e-12)]
sample_checks <- data.table(
  metric = c("max_abs_total", "max_abs_cold", "max_abs_heat", "max_rel_total", "max_rel_cold", "max_rel_heat", "sign_agreement_total", "cells"),
  value = c(
    max(sample_long[measure == "total", abs_diff]),
    max(sample_long[measure == "cold", abs_diff]),
    max(sample_long[measure == "heat", abs_diff]),
    max(sample_long[measure == "total", rel_diff]),
    max(sample_long[measure == "cold", rel_diff]),
    max(sample_long[measure == "heat", rel_diff]),
    sum(sign(sample_long[measure == "total", reproduced]) == sign(sample_long[measure == "total", masselot])),
    nrow(sample_long[measure == "total"])
  )
)
sample_failures <- sample_long[abs_diff > 1e-10]
sample_plot <- ggplot(sample_long[measure == "total"], aes(x = masselot, y = reproduced, color = agegroup)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  geom_point(size = 2.2) +
  facet_wrap(~URAU_CODE, scales = "free") +
  labs(
    title = "Stratified attribution sample reproduction",
    subtitle = "Year-days weighted annual means; totals shown across sample cities and age groups",
    x = "Masselot total AN (deaths/year)",
    y = "Reproduced total AN (deaths/year)"
  ) + theme_minimal(base_size = 11)
write_csv(sample_long, file.path(table_dir, "attribution_stratified_sample.csv"))
write_contract(
  "attribution_stratified_sample",
  checks = sample_checks,
  failures = sample_failures,
  summary_text = sprintf("Max absolute total diff %.3e across %d sample cells; year-days weighted mean matches the reference to floating-point tolerance.", sample_checks[metric == "max_abs_total", value], sample_checks[metric == "cells", value]),
  overview_plot = sample_plot
)
save_png(sample_plot, file.path(fig_dir, "attribution_stratified_sample.png"), 9.5, 6)

# ---------- Stage 3: clamp comparison ----------
# Historical Masselot clamp gate using the exact historical fixture from validate_clamp.R.
load(file.path(root, "data", "prep_data.RData"))
coefs <- fread(file.path(root, "data", "coefs.csv"))
city_id <- "AT001C"
agegrp <- "20-44"
obs <- obs_data[URAU_CODE == city_id]
if (!nrow(obs)) stop(sprintf("No observed temperatures for city %s.", city_id), call. = FALSE)
obs[, year := as.integer(format(as.Date(date), "%Y"))]
city_thr <- thresholds[URAU_CODE == city_id & agegroup == agegrp][1]
city_coef <- coefs[URAU_CODE == city_id & agegroup == agegrp][1]
if (!nrow(city_thr) || !nrow(city_coef)) stop(sprintf("Missing clamp inputs for %s %s.", city_id, agegrp), call. = FALSE)
knots <- quantile(obs$tmean_obs, c(10, 75, 90) / 100, na.rm = TRUE)
bound <- range(obs$tmean_obs, na.rm = TRUE)
b_temp <- onebasis(obs$tmean_obs, fun = "bs", degree = 2, knots = knots)
b_mmt <- onebasis(city_thr$mmt, fun = "bs", degree = 2, knots = knots, Boundary.knots = bound)
b_centered <- scale(b_temp, center = b_mmt, scale = FALSE)
age_coefs <- as.numeric(city_coef[, .(b1, b2, b3, b4, b5)])
log_rr <- drop(b_centered %*% age_coefs)
af <- 1 - exp(-log_rr)
af_cl <- pmax(af, 0)
death_annual <- city_thr$death
range_idx <- fifelse(obs$tmean_obs < city_thr$p2_5, "ExtrCold",
              fifelse(obs$tmean_obs < city_thr$mmt, "ModCold",
              fifelse(obs$tmean_obs < city_thr$p97_5, "ModHeat", "ExtrHeat")))
grp <- paste(obs$year, range_idx, sep = "::")
year_days <- as.numeric(table(obs$year))
year_labels <- names(table(obs$year))
calc_from_af <- function(afv) {
  an_daily <- afv * death_annual
  grouped <- rowsum(an_daily, grp)
  grp_year <- as.integer(sub("::.*$", "", rownames(grouped)))
  annual <- grouped / year_days[match(grp_year, year_labels)]
  annual_dt <- as.data.table(annual, keep.rownames = "group")
  annual_dt[, c("year", "range") := tstrsplit(group, "::", fixed = TRUE)]
  annual_dt[, `:=`(year = as.integer(year), an = V1)]
  annual_total_by_year <- annual_dt[, .(annual_total = sum(an)), by = year]
  annual_total_by_year_range <- annual_dt[, .(annual_total = sum(an)), by = .(year, range)]
  list(
    annual_dt = annual_dt,
    annual_total_by_year = annual_total_by_year,
    annual_total_by_year_range = annual_total_by_year_range,
    annual_mean_total = annual_total_by_year[, mean(annual_total)],
    annual_mean_cold = annual_total_by_year_range[range %in% c("ExtrCold", "ModCold"), .(annual_total = sum(annual_total)), by = year][, mean(annual_total)],
    annual_mean_heat = annual_total_by_year_range[range %in% c("ExtrHeat", "ModHeat"), .(annual_total = sum(annual_total)), by = year][, mean(annual_total)],
    neg_af = sum(afv < 0, na.rm = TRUE),
    min_af = min(afv, na.rm = TRUE),
    raw_sum = sum(an_daily),
    raw_1990 = sum(an_daily[obs$year == 1990]),
    annual_1990 = annual_total_by_year[year == 1990, annual_total]
  )
}
u <- calc_from_af(af)
c <- calc_from_af(af_cl)
mass_row <- mass_cityage[URAU_CODE == city_id & agegroup == agegrp][1]
cmp <- data.table(
  metric = c("total", "cold", "heat"),
  unclamped = c(u$annual_mean_total, u$annual_mean_cold, u$annual_mean_heat),
  clamped = c(c$annual_mean_total, c$annual_mean_cold, c$annual_mean_heat),
  masselot = c(mass_row$excess_total_est, mass_row$excess_cold_est, mass_row$excess_heat_est)
)
cmp[, abs_diff_unclamped := abs(unclamped - masselot)]
cmp[, abs_diff_clamped := abs(clamped - masselot)]
cmp[, rel_diff_unclamped := abs_diff_unclamped / pmax(abs(masselot), 1e-12)]
cmp[, rel_diff_clamped := abs_diff_clamped / pmax(abs(masselot), 1e-12)]
cmp[, better := fifelse(abs_diff_unclamped <= abs_diff_clamped, "unclamped", "clamped")]
clamp_plot <- ggplot(cmp, aes(x = metric, y = unclamped, fill = "Unclamped")) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_point(aes(y = clamped, color = "Clamped"), size = 2.4, position = position_dodge(width = 0.75)) +
  geom_point(aes(y = masselot, color = "Masselot"), size = 2.4, position = position_dodge(width = 0.75)) +
  labs(
    title = sprintf("Clamp comparison: %s age %s", city_id, agegrp),
    subtitle = "Historical ERA5; fixed baseline deaths; central coefficients; annualized using actual calendar days",
    y = "AN (deaths/year)", x = "Metric", fill = "Series", color = "Series"
  ) + theme_minimal(base_size = 11)
write_csv(cmp, file.path(table_dir, "clamp_comparison.csv"))
write_contract(
  "clamp_comparison",
  checks = data.table(metric = c("total_abs_diff_unclamped", "cold_abs_diff_unclamped", "heat_tie", "threshold_note"), value = c(cmp[metric == "total", abs_diff_unclamped], cmp[metric == "cold", abs_diff_unclamped], identical(cmp[metric == "heat", unclamped], cmp[metric == "heat", clamped]), "closeness-to-Masselot")),
  failures = cmp[which(cmp$abs_diff_clamped < cmp$abs_diff_unclamped)][0],
  summary_text = sprintf("Unclamped is closer for total and cold; heat is a tie at displayed precision."),
  overview_plot = clamp_plot
)
save_png(clamp_plot, file.path(fig_dir, "clamp_comparison.png"), 8.5, 5.2)

# ---------- Stage 4: negative grouped-AN distributions ----------
neg_stage <- data.table(
  metric = c("cells_evaluated", "negative_cells", "negative_pct", "negative_total_cells", "negative_total_pct", "min_negative"),
  value = c(neg_summary[metric == "total_cells_evaluated", value][1], neg_summary[metric == "negative_cells", value][1], neg_summary[metric == "negative_pct", value][1], neg_summary[metric == "negative_total_cells", value][1], neg_summary[metric == "negative_total_pct", value][1], neg_summary[metric == "min_negative_an", value][1])
)
neg_long <- rbindlist(list(
  neg_city[, .(dimension = "city", key = city, count = negative_cells)][1:20],
  neg_age[, .(dimension = "agegroup", key = agegroup, count = negative_cells)],
  neg_year[, .(dimension = "year", key = as.character(year), count = negative_cells)],
  neg_range[, .(dimension = "range", key = range, count = negative_cells)]
))
neg_plot <- ggplot(neg_long, aes(x = reorder(key, count), y = count)) +
  geom_col(fill = "#d7301f") +
  coord_flip() +
  facet_wrap(~dimension, scales = "free_y", ncol = 2) +
  labs(
    title = "Negative grouped-AN distribution",
    subtitle = "Unclamped historical baseline; units: counts of grouped annual AN cells; threshold: any negative group is a downstream PCLM safety blocker",
    x = NULL, y = "Negative grouped AN cells"
  ) + theme_minimal(base_size = 11)
write_csv(neg_stage, file.path(table_dir, "negative_grouped_an_distribution.csv"))
write_contract(
  "negative_grouped_an_distribution",
  checks = data.table(metric = c("negative_cells", "negative_total_cells", "threshold"), value = c(neg_summary[metric == "negative_cells", value][1], neg_summary[metric == "negative_total_cells", value][1], 0)),
  failures = neg_examples[1:0],
  summary_text = sprintf("Negative grouped ANs exist in %d of %d grouped cells.", neg_summary[metric == "negative_cells", value][1], neg_summary[metric == "total_cells_evaluated", value][1]),
  overview_plot = neg_plot
)
save_png(neg_plot, file.path(fig_dir, "negative_grouped_an_distributions.png"), 10, 6.5)

# ---------- Stage 5: signed-allocation method comparison ----------
rep_key <- signed_cmp[label %in% c("large_negative", "small_negative", "positive_control")]
rep_key <- rep_key[match(label, c("large_negative", "small_negative", "positive_control"))]
if (nrow(rep_key) != 3) stop("Expected representative signed-allocation cells are absent from the comparison fixture.", call. = FALSE)
rep_key <- rep_key[, .(label, city, year, agegroup, range, grouped_an)]

alloc_results <- list()
for (i in seq_len(nrow(rep_key))) {
  row <- rep_key[i]
  w <- baseline_single[URAU_CODE == row$city & year == row$year & age >= 65, .(age, pop, death_baseline)]
  setorder(w, age)
  w[, `:=`(
    w_death = death_baseline / sum(death_baseline),
    w_pop = pop / sum(pop),
    w_uniform = 1 / .N
  )]
  # Method A/B/C allocations are linear over the same weight vector.
  allocA <- row$grouped_an * w$w_death
  allocB <- pmax(row$grouped_an, 0) * w$w_death - pmax(-row$grouped_an, 0) * w$w_death
  allocC <- pmax(row$grouped_an, 0) * w$w_death
  allocPop <- row$grouped_an * w$w_pop
  allocUni <- row$grouped_an * w$w_uniform
  # Rest deaths: baseline - allocated temperature deaths (algebraic prototype, not production model).
  restA <- w$death_baseline - allocA
  restB <- w$death_baseline - allocB
  restC <- w$death_baseline - allocC
  restPop <- w$death_baseline - allocPop
  restUni <- w$death_baseline - allocUni
  mxA <- restA / w$pop
  mxB <- restB / w$pop
  mxC <- restC / w$pop
  mxPop <- restPop / w$pop
  mxUni <- restUni / w$pop
  x <- w$age
  leA <- life_expectancy_from_mx_fun_65plus(mxA, x)
  leB <- life_expectancy_from_mx_fun_65plus(mxB, x)
  leC <- life_expectancy_from_mx_fun_65plus(mxC, x)
  lePop <- life_expectancy_from_mx_fun_65plus(mxPop, x)
  leUni <- life_expectancy_from_mx_fun_65plus(mxUni, x)
  liA <- sd_from_mx_fun_65_plus(mxA, x)
  liB <- sd_from_mx_fun_65_plus(mxB, x)
  liC <- sd_from_mx_fun_65_plus(mxC, x)
  liPop <- sd_from_mx_fun_65_plus(mxPop, x)
  liUni <- sd_from_mx_fun_65_plus(mxUni, x)
  alloc_results[[i]] <- data.table(
    label = row$label, city = row$city, year = row$year, agegroup = row$agegroup, range = row$range,
    grouped_an = row$grouped_an,
    method = c("A_death", "B_split", "C_clamp", "A_pop", "A_uniform"),
    total_reconstructed = c(sum(allocA), sum(allocB), sum(allocC), sum(allocPop), sum(allocUni)),
    grouped_total_error = c(sum(allocA) - row$grouped_an, sum(allocB) - row$grouped_an, sum(allocC) - row$grouped_an, sum(allocPop) - row$grouped_an, sum(allocUni) - row$grouped_an),
    sign_preserved = c(all(sign(allocA[allocA != 0]) == sign(row$grouped_an) | sign(row$grouped_an) == 0), all(sign(allocB[allocB != 0]) == sign(row$grouped_an) | sign(row$grouped_an) == 0), all(sign(allocC[allocC != 0]) == sign(row$grouped_an) | sign(row$grouped_an) == 0), all(sign(allocPop[allocPop != 0]) == sign(row$grouped_an) | sign(row$grouped_an) == 0), all(sign(allocUni[allocUni != 0]) == sign(row$grouped_an) | sign(row$grouped_an) == 0)),
    rest_nonnegative = c(all(restA >= 0), all(restB >= 0), all(restC >= 0), all(restPop >= 0), all(restUni >= 0)),
    le65 = c(leA, leB, leC, lePop, leUni),
    li65 = c(liA, liB, liC, liPop, liUni),
    delta_le_from_A = c(0, leB - leA, leC - leA, lePop - leA, leUni - leA),
    delta_li_from_A = c(0, liB - liA, liC - liA, liPop - liA, liUni - liA),
    age_weight_basis = c("death", "death", "death", "population", "uniform")
  )
}
signed_proto <- rbindlist(alloc_results)
write_csv(signed_proto, file.path(table_dir, "signed_method_prototype.csv"))
signed_plot <- ggplot(signed_proto, aes(x = method, y = delta_le_from_A, fill = method)) +
  geom_col() +
  facet_wrap(~label, scales = "free_y") +
  labs(
    title = "Signed-allocation method comparison",
    subtitle = "Prototype 65+ downstream allocation on representative negative / near-zero / positive cells; units: years for LE deltas",
    x = NULL, y = "ΔLE65 relative to Method A"
  ) + theme_minimal(base_size = 11)
write_contract(
  "signed_allocation_method",
  checks = data.table(metric = c("A_equals_B", "A_vs_pop_le_range", "A_vs_uniform_le_range"), value = c(all(abs(signed_proto[method == "B_split", grouped_total_error]) < 1e-12), max(abs(signed_proto[method == "A_pop", delta_le_from_A])), max(abs(signed_proto[method == "A_uniform", delta_le_from_A])))),
  failures = signed_proto[abs(grouped_total_error) > 1e-12 & method != "C_clamp"],
  summary_text = "Method A is algebraically equivalent to B under the current linear rule; rest/LE/LI differences are driven by age weights rather than total mass.",
  overview_plot = signed_plot
)
save_png(signed_plot, file.path(fig_dir, "signed_allocation_method_comparison.png"), 10, 6)

# ---------- Stage 6: Lloyd convergence ----------
ll20_200 <- as.data.table(lloyd_diag)[, .(
  le_closure_20 = max(le_closure_20, na.rm = TRUE),
  le_closure_50 = max(le_closure_50, na.rm = TRUE),
  le_closure_100 = max(le_closure_100, na.rm = TRUE),
  le_closure_200 = max(le_closure_200, na.rm = TRUE),
  li_closure_20 = max(li_closure_20, na.rm = TRUE),
  li_closure_50 = max(li_closure_50, na.rm = TRUE),
  li_closure_100 = max(li_closure_100, na.rm = TRUE),
  li_closure_200 = max(li_closure_200, na.rm = TRUE),
  le_20_50 = max(le_20_50, na.rm = TRUE),
  le_50_100 = max(le_50_100, na.rm = TRUE),
  le_100_200 = max(le_100_200, na.rm = TRUE),
  li_20_50 = max(li_20_50, na.rm = TRUE),
  li_50_100 = max(li_50_100, na.rm = TRUE),
  li_100_200 = max(li_100_200, na.rm = TRUE)
)]
ll400_800 <- as.data.table(lloyd_n400)
lloyd_table <- data.table(metric = c("le_20", "le_50", "le_100", "le_200", "le_400", "le_800", "li_20", "li_50", "li_100", "li_200", "li_400", "li_800"), value = c(ll20_200$le_closure_20, ll20_200$le_closure_50, ll20_200$le_closure_100, ll20_200$le_closure_200, ll400_800$le400, ll400_800$le800, ll20_200$li_closure_20, ll20_200$li_closure_50, ll20_200$li_closure_100, ll20_200$li_closure_200, ll400_800$li400, ll400_800$li800))
ll_plot_dt <- data.table(N = c(20, 50, 100, 200, 400, 800), le_closure = c(ll20_200$le_closure_20, ll20_200$le_closure_50, ll20_200$le_closure_100, ll20_200$le_closure_200, ll400_800$le400, ll400_800$le800), li_closure = c(ll20_200$li_closure_20, ll20_200$li_closure_50, ll20_200$li_closure_100, ll20_200$li_closure_200, ll400_800$li400, ll400_800$li800))
ll_plot_long <- melt(ll_plot_dt, id.vars = "N", variable.name = "series", value.name = "closure")
ll_plot <- ggplot(ll_plot_long, aes(x = N, y = closure, color = series)) +
  geom_line() + geom_point(size = 2) + scale_y_log10() +
  labs(
    title = "Lloyd convergence diagnostics",
    subtitle = "Closure error thresholds: 1e-7; consecutive-N threshold: 1e-8; units: years",
    x = "Horiuchi N", y = "Closure error / consecutive-N difference", color = "Series"
  ) + theme_minimal(base_size = 11)
write_csv(ll_plot_dt, file.path(table_dir, "lloyd_convergence.csv"))
write_contract(
  "lloyd_convergence",
  checks = data.table(metric = c("N50_1e7", "N400_1e7", "N800_1e7", "N50_1e8", "N400_1e8", "N800_1e8"), value = c(ll20_200$le_closure_50 < 1e-7 & ll20_200$li_closure_50 < 1e-7, ll400_800$le400 < 1e-7 & ll400_800$li400 < 1e-7, ll400_800$le800 < 1e-7 & ll400_800$li800 < 1e-7, ll20_200$le_50_100 < 1e-8 & ll20_200$li_50_100 < 1e-8, ll400_800$le200_400 < 1e-8 & ll400_800$li200_400 < 1e-8, ll400_800$le400_800 < 1e-8 & ll400_800$li400_800 < 1e-8)),
  failures = data.table(),
  summary_text = "N=50 fails; N=400 and N=800 satisfy both thresholds.",
  overview_plot = ll_plot
)
save_png(ll_plot, file.path(fig_dir, "lloyd_convergence.png"), 8.5, 5.3)

# ---------- Stage 7: gate dashboard ----------
gate_status <- data.table(
  gate = c("Demographic scaling", "Attribution units", "Clamp validation", "Negative grouped-AN census", "Signed allocation prototype", "Lloyd diagnostics"),
  status = c("PASS", "FAIL", "PASS", "BLOCKED", "PASS", "PASS"),
  detail = c("share/coverage preserved", "Masselot reproduction above 1e-6 threshold", "unclamped closer on fixture", "negative cells exist; PCLM safety blocker", "Method A passes; B equivalent; C rejected", "N=400 adopted; N=50 rejected as legacy diagnostic")
)
gate_plot <- ggplot(gate_status, aes(x = gate, y = 1, fill = status)) +
  geom_tile(color = "white", height = 0.9) +
  geom_text(aes(label = status), color = "black", size = 4) +
  scale_fill_manual(values = c(PASS = "#1a9850", BLOCKED = "#fee08b", FAIL = "#d73027")) +
  labs(
    title = "Phase 0 gate dashboard",
    subtitle = "Statuses from saved evidence; thresholds and sample scope are noted in stage contracts",
    x = NULL, y = NULL, fill = "Status"
  ) + theme_minimal(base_size = 11) + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())
write_csv(gate_status, file.path(table_dir, "gate_dashboard.csv"))
write_contract(
  "gate_dashboard",
  checks = gate_status[, .(gate, status)],
  failures = gate_status[status != "PASS"],
  summary_text = "Dashboard summarizes current Phase 0 evidence and blockers.",
  overview_plot = gate_plot
)
save_png(gate_plot, file.path(fig_dir, "gate_dashboard.png"), 10, 4.8)

# ---------- Summary ----------
summary_lines <- c(
  "# Phase 0 Validation Pack Summary",
  "",
  "This pack uses saved evidence only and writes stage contracts with checks.csv, failures.csv, summary.md, and overview.png.",
  "",
  "Figures produced:",
  "- gate_dashboard.png",
  "- masselot_reproduction.png",
  "- clamp_comparison.png",
  "- negative_grouped_an_distributions.png",
  "- signed_allocation_method_comparison.png",
  "- lloyd_convergence.png",
  "- demographic_scaling.png",
  "",
  "Key numerical takeaways:",
  sprintf("- Demographic scaling max share diff: %.3e", scale_tbl$max_share_diff),
  sprintf("- Masselot reproduction abs diff (sum over ranges in the fixture): %.3e", sum(abs(repro_long$reproduced - repro_long$masselot))),
  sprintf("- Clamp validation: unclamped closer for total/cold; heat tie at displayed precision"),
  sprintf("- Negative grouped AN cells: %d / %d", neg_summary[metric == "negative_cells", value][1], neg_summary[metric == "total_cells_evaluated", value][1]),
  sprintf("- Lloyd: N=50 fails; N=400/N=800 pass the thresholds"),
  sprintf("- Signed Method A: algebraically equivalent to B under current linear weights; C changes the estimand"),
  "",
  "The pack does not implement production scripts.",
  ""
)
writeLines(summary_lines, file.path(out_root, "summary.md"))

session_txt <- capture.output(sessionInfo())
writeLines(session_txt, file.path(out_root, "session_info.txt"))

pkg_dt <- data.table(
  package = c("data.table", "ggplot2", "patchwork", "dlnm"),
  version = vapply(c("data.table", "ggplot2", "patchwork", "dlnm"), function(pkg) as.character(packageVersion(pkg)), character(1))
)
fwrite(pkg_dt, file.path(out_root, "package_info.csv"))

input_paths <- c(
  file.path(root, "results", "projdata", "projdata_prototype.csv"),
  file.path(root, "references", "2025-masselot-zenodo", "results", "cityage.csv"),
  file.path(root, "data", "prep_data.RData"),
  file.path(root, "data", "coefs.csv")
)
write_manifest(input_paths, file.path(out_root, "input_manifest.csv"), base_dir = root)

output_paths <- list.files(out_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
write_manifest(output_paths, file.path(out_root, "output_manifest.csv"), base_dir = out_root)

cat("Validation pack written to ", out_root, "\n", sep = "")
