#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source("R_pipeline/functions/pclm_utils.R")

message("\n[00] Building Madrid SSP3-7.0 demographic tables...")

out_dir <- "results/phase1_madrid"
check_dir <- "results/checks"
fig_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

grouped_file <- file.path(out_dir, "00_demography_grouped.csv")
single_file <- file.path(out_dir, "00_demography_single_age.csv")
checks_file <- file.path(check_dir, "00_demography_checks.csv")
failures_file <- file.path(check_dir, "00_demography_failures.csv")
fig_file <- file.path(fig_dir, "00_demography_diagnostic.png")

city_id <- "ES001C"
city_name <- "Madrid"
ssp_target <- "3"
baseline_year <- 2020L
future_years <- 2020:2099
target_agegroups <- c("65-74", "75-84", "85+")
age_map <- data.table(
  age_band = c("65--69", "70--74", "75--79", "80--84", "85--89", "90--94", "95--99", "100+"),
  agegroup = c(rep("65-74", 2L), rep("75-84", 2L), rep("85+", 4L)),
  age_start = c(65L, 70L, 75L, 80L, 85L, 90L, 95L, 100L),
  nlast = c(rep(5L, 7L), 1L)
)

city_meta <- fread("data/city_results.csv")
city_meta <- unique(city_meta[URAU_CODE == city_id & LABEL == city_name & agegroup %in% target_agegroups])
if (nrow(city_meta) != 3L || !all(unique(city_meta$agegroup) %in% target_agegroups)) {
  stop(sprintf("Madrid identifier verification failed for %s; expected 3 age-group rows for %s.", city_id, city_name), call. = FALSE)
}

city_meta <- city_meta[match(target_agegroups, agegroup)]
if (anyNA(city_meta$agegroup)) {
  stop(sprintf("Madrid identifier verification failed: missing one of age groups %s.", paste(target_agegroups, collapse = ", ")), call. = FALSE)
}

pop_raw <- fread("data/wittgenstein_pop.csv")
assr_raw <- fread("data/wittgenstein_assr.csv")

pop_raw <- pop_raw[ssp == as.integer(ssp_target) & age %in% age_map$age_band]
assr_raw <- assr_raw[ssp == as.integer(ssp_target) & age %in% age_map$age_band]

if (!nrow(pop_raw) || !nrow(assr_raw)) {
  stop("SSP3 demographic source tables are empty after filtering.", call. = FALSE)
}

parse_year_start <- function(x) as.integer(sub("^([0-9]{4}).*$", "\\1", x))

pop_raw[, year_start := as.integer(year)]
assr_raw[, year_start := parse_year_start(period)]

if (anyNA(pop_raw$year_start) || anyNA(assr_raw$year_start)) {
  stop("Failed to parse SSP3 year anchors from the demographic sources.", call. = FALSE)
}

country_5y <- merge(
  pop_raw[, .(CNTR_CODE, cntr_name, ssp, year_start, age_band = as.character(age), sex, pop = as.numeric(pop))],
  assr_raw[, .(CNTR_CODE, age_band = as.character(age), sex, year_start, assr = as.numeric(assr))],
  by = c("CNTR_CODE", "age_band", "sex", "year_start"),
  all = FALSE,
  sort = FALSE
)

if (anyNA(country_5y$pop) || anyNA(country_5y$assr)) {
  stop("Merged SSP3 demographic source contains NA pop/assr values.", call. = FALSE)
}
if (any(country_5y$pop < 0) || any(!is.finite(country_5y$pop)) || any(!is.finite(country_5y$assr))) {
  stop("Invalid SSP3 demographic source values detected.", call. = FALSE)
}

country_5y[, death := pop * (1 - assr)]
country_5y[, `:=`(
  year_end = pmin(year_start + 4L, 2100L)
)]

expand_annual <- function(dt) {
  rbindlist(lapply(seq_len(nrow(dt)), function(i) {
    row <- dt[i]
    yrs <- seq(row$year_start, row$year_end)
    data.table(
      CNTR_CODE = row$CNTR_CODE,
      cntr_name = row$cntr_name,
      ssp = row$ssp,
      year = yrs,
      age_band = row$age_band,
      age_start = age_map$age_start[match(row$age_band, age_map$age_band)],
      sex = row$sex,
      pop = row$pop,
      death = row$death
    )
  }), use.names = TRUE, fill = TRUE)
}

country_annual_5y <- expand_annual(country_5y)

country_annual <- country_annual_5y[, .(
  pop = sum(pop),
  death = sum(death)
), by = .(CNTR_CODE, cntr_name, ssp, year, age_band, age_start)]

country_annual <- merge(country_annual, age_map[, .(age_band, agegroup, age_start)], by = c("age_band", "age_start"), all.x = TRUE, sort = FALSE)
if (anyNA(country_annual$agegroup)) {
  stop("Failed to map country age bands to Madrid age groups.", call. = FALSE)
}

country_agegroup <- country_annual[, .(
  country_pop = sum(pop),
  country_death = sum(death)
), by = .(year, agegroup)]

country_agegroup_2020 <- country_agegroup[year == baseline_year]
if (nrow(country_agegroup_2020) != length(target_agegroups)) {
  stop("Baseline year 2020 is missing one or more required age groups in the SSP3 source.", call. = FALSE)
}

city_base <- unique(city_meta[, .(agegroup, city_pop = agepop, city_death = death)])
share_tbl <- merge(city_base, country_agegroup_2020, by = "agegroup", all.x = TRUE, sort = FALSE)
if (anyNA(share_tbl$country_pop) || anyNA(share_tbl$country_death)) {
  stop("Could not compute Madrid baseline demographic shares from the SSP3 source.", call. = FALSE)
}
if (any(share_tbl$country_pop <= 0) || any(share_tbl$country_death <= 0)) {
  stop("Nonpositive baseline country demographic totals prevent share construction.", call. = FALSE)
}

share_tbl[, `:=`(
  pop_share = city_pop / (country_pop * 1000),
  death_share = city_death / (country_death * 1000)
)]

if (any(!is.finite(share_tbl$pop_share)) || any(!is.finite(share_tbl$death_share))) {
  stop("Madrid baseline shares are not finite.", call. = FALSE)
}

country_single <- rbindlist(lapply(sort(unique(country_annual$year)), function(yr) {
  yr_dt <- country_annual[year == yr]
  single_rows <- list()
  for (kind in c("pop", "death")) {
    grouped_vals <- yr_dt[, .(value = sum(get(kind))), by = .(age_start, agegroup)]
    grouped_vals <- grouped_vals[order(age_start)]
    x <- grouped_vals$age_start
    if (any(!is.finite(grouped_vals$value)) || any(grouped_vals$value < 0)) {
      stop(sprintf("Invalid country grouped %s values for year %d.", kind, yr), call. = FALSE)
    }
    fit <- if (sum(grouped_vals$value) == 0) {
      rep(0, pclm_expected_length(x, 1L))
    } else {
      pclm_disaggregate_nonnegative(x = x, y = grouped_vals$value, nlast = 1L)
    }
    if (length(fit) != 36L) {
      stop(sprintf("Unexpected PCLM output length for year %d and %s.", yr, kind), call. = FALSE)
    }
    fit_age <- data.table(year = yr, age = 65:100, value = as.numeric(fit), kind = kind)
    single_rows[[length(single_rows) + 1L]] <- fit_age
  }
  pop_dt <- single_rows[[1L]]
  death_dt <- single_rows[[2L]]
  out <- merge(pop_dt, death_dt, by = c("year", "age"), suffixes = c("_pop", "_death"))
  out[, `:=`(
    pop = value_pop,
    death = value_death
  )]
  out[, c("value_pop", "value_death") := NULL]
  out[, agegroup := fifelse(age <= 74L, "65-74", fifelse(age <= 84L, "75-84", "85+"))]
  out
}))

country_single_agegroup <- country_single[, .(
  country_pop = sum(pop),
  country_death = sum(death)
), by = .(year, agegroup, age)]

country_single_agegroup <- merge(
  country_single_agegroup,
  country_agegroup[, .(year, agegroup, country_agegroup_pop = country_pop, country_agegroup_death = country_death)],
  by = c("year", "agegroup"),
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(country_single_agegroup$country_agegroup_pop) || anyNA(country_single_agegroup$country_agegroup_death)) {
  stop("Failed to attach country agegroup totals to single-age country table.", call. = FALSE)
}

country_single_agegroup[, `:=`(
  pop_weight = fifelse(country_agegroup_pop > 0, country_pop / country_agegroup_pop, NA_real_),
  death_weight = fifelse(country_agegroup_death > 0, country_death / country_agegroup_death, NA_real_)
)]

country_single_agegroup[, `:=`(
  pop_weight = pop_weight / sum(pop_weight),
  death_weight = death_weight / sum(death_weight)
), by = .(year, agegroup)]

if (anyNA(country_single_agegroup$pop_weight) || anyNA(country_single_agegroup$death_weight)) {
  stop("Single-age country weights could not be constructed.", call. = FALSE)
}

city_grouped <- merge(country_agegroup, share_tbl[, .(agegroup, pop_share, death_share)], by = "agegroup", all.x = TRUE, sort = FALSE)
city_grouped[, `:=`(
  geo_id = city_id,
  label = city_name,
  country_code = "ES",
  cntr_name = "Spain",
  ssp = as.integer(ssp_target),
  city_pop = country_pop * pop_share * 1000,
  city_death = country_death * death_share * 1000
)]

city_grouped <- city_grouped[year %in% future_years, .(
  geo_id, label, country_code, cntr_name, ssp, year, agegroup,
  pop = as.numeric(city_pop),
  death = as.numeric(city_death),
  country_pop = as.numeric(country_pop * 1000),
  country_death = as.numeric(country_death * 1000),
  pop_share,
  death_share
)]

city_single <- merge(
  city_grouped,
  country_single_agegroup[, .(year, agegroup, age, pop_weight, death_weight)],
  by = c("year", "agegroup"),
  all.x = TRUE,
  sort = FALSE
)

city_single[, `:=`(
  geo_id = city_id,
  label = city_name,
  country_code = "ES",
  cntr_name = "Spain",
  ssp = as.integer(ssp_target),
  pop = pop * pop_weight,
  death = death * death_weight
)]

city_single <- city_single[year %in% future_years, .(
  geo_id, label, country_code, cntr_name, ssp, year, agegroup, age,
  pop = as.numeric(pop),
  death = as.numeric(death),
  country_pop = as.numeric(country_pop),
  country_death = as.numeric(country_death),
  pop_share, death_share,
  pop_weight, death_weight
)]

setorder(city_grouped, year, agegroup)
setorder(city_single, year, agegroup, age)

recon_check <- merge(
  city_grouped[, .(year, agegroup, grouped_pop = pop, grouped_death = death)],
  city_single[, .(single_pop = sum(pop), single_death = sum(death)), by = .(year, agegroup)],
  by = c("year", "agegroup"),
  all.x = TRUE,
  sort = FALSE
)
recon_check[, `:=`(
  pop_abs_diff = abs(grouped_pop - single_pop),
  death_abs_diff = abs(grouped_death - single_death)
)]

full_group_grid <- CJ(year = future_years, agegroup = target_agegroups)
full_single_grid <- CJ(year = future_years, age = 65:100)

grouped_keys <- unique(city_grouped[, .(year, agegroup)])
single_keys <- unique(city_single[, .(year, agegroup, age)])

checks <- data.table(
  check_name = c(
    "madrid_identifier_unique",
    "grouped_primary_keys_complete",
    "single_primary_keys_complete",
    "population_finite_nonnegative",
    "death_finite_nonnegative",
    "grouped_equals_single_sum",
    "baseline_shares_finite"
  ),
  status = c(
    if (nrow(city_meta) == 3L) "PASS" else "FAIL",
    if (nrow(grouped_keys) == nrow(full_group_grid) && nrow(unique(city_grouped, by = c("year", "agegroup"))) == nrow(city_grouped)) "PASS" else "FAIL",
    if (nrow(single_keys) == nrow(full_single_grid) && nrow(unique(city_single, by = c("year", "age"))) == nrow(city_single)) "PASS" else "FAIL",
    if (!any(!is.finite(city_grouped$pop)) && !any(city_grouped$pop < 0) && !any(!is.finite(city_single$pop)) && !any(city_single$pop < 0)) "PASS" else "FAIL",
    if (!any(!is.finite(city_grouped$death)) && !any(city_grouped$death < 0) && !any(!is.finite(city_single$death)) && !any(city_single$death < 0)) "PASS" else "FAIL",
    if (max(recon_check$pop_abs_diff) <= 1e-9 && max(recon_check$death_abs_diff) <= 1e-9) "PASS" else "FAIL",
    if (!any(!is.finite(share_tbl$pop_share)) && !any(!is.finite(share_tbl$death_share))) "PASS" else "FAIL"
  ),
  value = c(
    nrow(city_meta),
    nrow(grouped_keys),
    nrow(single_keys),
    sprintf("min_pop=%g; min_single_pop=%g", min(city_grouped$pop), min(city_single$pop)),
    sprintf("min_death=%g; min_single_death=%g", min(city_grouped$death), min(city_single$death)),
    sprintf("max_pop_diff=%0.3e; max_death_diff=%0.3e", max(recon_check$pop_abs_diff), max(recon_check$death_abs_diff)),
    sprintf("max_pop_share=%0.6f; max_death_share=%0.6f", max(share_tbl$pop_share), max(share_tbl$death_share))
  ),
  threshold = c(
    "unique city id with 3 age-group rows",
    sprintf("%d rows", nrow(full_group_grid)),
    sprintf("%d rows", nrow(full_single_grid)),
    "finite and >= 0",
    "finite and >= 0",
    "<= 1e-9",
    "finite shares"
  )
)

failures <- data.table()
if (any(checks$status == "FAIL")) {
  failures <- rbindlist(lapply(which(checks$status == "FAIL"), function(i) {
    data.table(
      failing_check = checks$check_name[i],
      observed_value = checks$value[i],
      expected_bound = checks$threshold[i]
    )
  }), fill = TRUE)
}

fwrite(city_grouped, grouped_file)
fwrite(city_single, single_file)
fwrite(checks, checks_file)
if (nrow(failures)) {
  fwrite(failures, failures_file)
} else {
  if (file.exists(failures_file)) file.remove(failures_file)
  file.create(failures_file)
}

plot_group <- city_grouped[, .(pop = sum(pop), death = sum(death)), by = .(year, agegroup)]
plot_group <- melt(plot_group, id.vars = c("year", "agegroup"), variable.name = "measure", value.name = "value")

p <- ggplot(plot_group, aes(x = year, y = value, color = agegroup)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~measure, scales = "free_y") +
  labs(
    title = "Madrid SSP3-7.0 demographic projection diagnostic",
    subtitle = sprintf("City %s (%s); piecewise-constant annualization from 5-year source snapshots", city_name, city_id),
    x = "Year",
    y = "Count"
  ) +
  theme_minimal(base_size = 11)

ggsave(fig_file, p, width = 10, height = 6, dpi = 160)

if (nrow(failures)) {
  stop(sprintf("00_demography.R failed %d invariant(s); see %s", nrow(failures), failures_file), call. = FALSE)
}

message("Saved grouped demography to ", grouped_file)
message("Saved single-age demography to ", single_file)
message("Saved checks to ", checks_file)
message("Saved diagnostic figure to ", fig_file)
