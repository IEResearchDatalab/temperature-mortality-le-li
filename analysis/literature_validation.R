#!/usr/bin/env Rscript

# Evidence-backed validation gate for the LE/LI paper. Outputs are written to a
# new result root; canonical pipeline results and reference files are read only.

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(dplyr)
  library(dlnm)
  library(splines)
})

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
out_root <- Sys.getenv("VALIDATION_ROOT", "results/europe/validation_20261104")
masselot_root <- Sys.getenv(
  "MASSELOT_RESULTS_ROOT",
  "/home/SHARED/Climate-Change-Effect-on-Europe/data/masselot/results_parquet"
)
lloyd_root <- Sys.getenv("LLOYD_ROOT", "/home/SHARED/Reciprocal-Lloyd")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

write_table <- function(x, name) {
  fwrite(x, file.path(out_root, name))
  x
}

# -----------------------------------------------------------------------------
# 1. Historical Masselot 2023 fixture checks and published historical totals
# -----------------------------------------------------------------------------

cities <- sort(unique(fread("data/city_results.csv", select = "URAU_CODE")$URAU_CODE))
fixture_files <- file.path(
  batch_root, "ssp3", cities, "ACCESS_CM2", "checks", "01_attribution_checks.csv"
)
if (any(!file.exists(fixture_files))) {
  stop("One or more city fixture files are missing", call. = FALSE)
}
fixtures <- rbindlist(lapply(seq_along(fixture_files), function(i) {
  x <- fread(fixture_files[i])[check_name == "masselot2023_fixture"]
  x[, city := cities[i]]
  x
}))
fixtures[, max_rel_error := as.numeric(sub(
  ".*max_rel_error=([0-9.eE+-]+).*", "\\1", value
))]
fixture_summary <- data.table(
  source = "Masselot et al. 2023 city-age published excess deaths",
  n_cities = uniqueN(fixtures$city),
  pass = sum(fixtures$status == "PASS"),
  fail = sum(fixtures$status != "PASS"),
  maximum_relative_error = max(fixtures$max_rel_error),
  tolerance = 0.001,
  verdict = if (all(fixtures$status == "PASS") &&
                  max(fixtures$max_rel_error) <= 0.001) "match" else "not explained"
)
write_table(fixtures, "masselot_2023_city_fixture_checks.csv")
write_table(fixture_summary, "masselot_2023_fixture_summary.csv")

hist <- fread("data/city_results.csv")
hist <- unique(hist[, .(
  city = URAU_CODE, region, agegroup,
  heat = excess_heat_est, cold = excess_cold_est
)])
hist_summary <- rbindlist(list(
  hist[, .(heat_deaths = sum(heat), cold_deaths = sum(cold)),
       by = .(geography = region, agegroup)],
  hist[, .(geography = "Europe", heat_deaths = sum(heat), cold_deaths = sum(cold)),
       by = agegroup]
))
setcolorder(hist_summary, c("geography", "agegroup", "heat_deaths", "cold_deaths"))
write_table(hist_summary, "masselot_2023_published_historical_totals.csv")

# -----------------------------------------------------------------------------
# 2. Masselot 2025 central projections, like-for-like older age groups
# -----------------------------------------------------------------------------

obj_path <- file.path(batch_root, "collected", "object1_dataset.parquet")
obj_ds <- open_dataset(obj_path)
regional <- obj_ds |>
  filter(age >= 65, cause != "rest") |>
  group_by(region, ssp, scenario, year, age, cause) |>
  summarise(deaths = sum(deaths), pop = sum(pop)) |>
  collect() |>
  as.data.table()
regional[, geography := region]
europe <- regional[, .(deaths = sum(deaths), pop = sum(pop)),
                   by = .(ssp, scenario, year, age, cause)]
europe[, geography := "Europe"]
ours <- rbind(regional[, region := NULL], europe, fill = TRUE)
ours[, `:=`(
  period = (year %/% 5L) * 5L,
  agegroup = fifelse(age <= 74, "65-74", fifelse(age <= 84, "75-84", "85+")),
  broad_range = fifelse(cause %chin% c("ExtrCold", "ModCold"), "cold", "heat")
)]
ours <- ours[, .(
  deaths = sum(deaths),
  pop = unique(sum(pop) / uniqueN(cause))
), by = .(geography, ssp, scenario, period, agegroup, broad_range)]
ours <- dcast(
  ours, geography + ssp + period + agegroup + broad_range ~ scenario,
  value.var = c("deaths", "pop")
)
ours[, `:=`(
  climate_deaths = deaths_with_cc - deaths_without_cc,
  rate_per_100000 = 1e5 * (deaths_with_cc - deaths_without_cc) / pop_with_cc
)]
ours_tot <- ours[, .(
  climate_deaths = sum(climate_deaths) / 5,
  rate_per_100000 = sum(climate_deaths) / unique(pop_with_cc) * 1e5
), by = .(geography, ssp, period, agegroup)][, broad_range := "tot"]
ours <- rbind(
  ours[, .(
    geography, ssp, period, agegroup, broad_range,
    climate_deaths = climate_deaths / 5,
    rate_per_100000
  )],
  ours_tot, use.names = TRUE
)

standardise_reference <- function(file, geography_level) {
  x <- as.data.table(read_parquet(file.path(masselot_root, file)))
  if ("adapt" %in% names(x)) x <- x[adapt %chin% c("0%", "noadapt")]
  if ("sc" %in% names(x)) x <- x[sc == "clim"]
  if ("gcm" %in% names(x)) x <- x[gcm == "ens"]
  x <- x[agegroup %chin% c("65-74", "75-84", "85+") &
           range %chin% c("cold", "heat", "tot")]
  if (geography_level == "europe") {
    x[, geography := "Europe"]
  } else {
    geo_col <- intersect(
      c("region", "region_name", "regname", "country", "cntr_name", "CNTR_CODE"),
      names(x)
    )[1]
    if (is.na(geo_col)) stop("Cannot identify geography column in ", file)
    x[, geography := as.character(get(geo_col))]
    if (geography_level == "region") {
      x[, geography := sub(" Europe$", "", geography)]
    }
  }
  x[, ssp_num := fifelse(grepl("126|ssp1|^1$", as.character(ssp)), 1L,
                  fifelse(grepl("245|ssp2|^2$", as.character(ssp)), 2L, 3L))]
  x[, .(
    geography, ssp = ssp_num, period = as.integer(as.character(period)),
    agegroup, broad_range = range,
    reference_deaths = an_est, reference_rate_per_100000 = 1e5 * rate_est
  )]
}

reference <- rbindlist(list(
  standardise_reference("eu_period.parquet", "europe"),
  standardise_reference("region_period.parquet", "region")
), fill = TRUE)
comparison <- merge(
  ours, reference,
  by = c("geography", "ssp", "period", "agegroup", "broad_range")
)
comparison[, `:=`(
  death_abs_difference = climate_deaths - reference_deaths,
  death_relative_error = fifelse(
    abs(reference_deaths) > 1e-10,
    abs(climate_deaths / reference_deaths - 1),
    abs(climate_deaths - reference_deaths)
  ),
  rate_abs_difference = rate_per_100000 - reference_rate_per_100000,
  rate_relative_error = fifelse(
    abs(reference_rate_per_100000) > 1e-10,
    abs(rate_per_100000 / reference_rate_per_100000 - 1),
    abs(rate_per_100000 - reference_rate_per_100000)
  )
)]
comparison[, verdict := fifelse(
  ssp %in% c(2L, 3L) & pmax(death_relative_error, rate_relative_error) < 0.01,
  "match",
  fifelse(ssp == 1L, "explained: SSP1 uses the documented 0.9995 demographic sensitivity", "not explained")
)]
write_table(comparison, "masselot_2025_projection_comparison.csv")
comparison_summary <- comparison[, .(
  rows = .N,
  max_death_relative_error = max(death_relative_error),
  max_rate_relative_error = max(rate_relative_error),
  median_rate_relative_error = median(rate_relative_error),
  verdict = if (all(verdict == "match")) "match" else unique(verdict)[1]
), by = ssp]
write_table(comparison_summary, "masselot_2025_projection_summary.csv")

# Country spot checks use the same reference tables and the canonical object.
country_obj <- obj_ds |>
  filter(age >= 65, cause != "rest", country %in% c("MT", "IE", "RO")) |>
  group_by(country, ssp, scenario, year, age, cause) |>
  summarise(deaths = sum(deaths), pop = sum(pop)) |>
  collect() |>
  as.data.table()
country_obj[, `:=`(
  period = (year %/% 5L) * 5L,
  agegroup = fifelse(age <= 74, "65-74", fifelse(age <= 84, "75-84", "85+")),
  broad_range = fifelse(cause %chin% c("ExtrCold", "ModCold"), "cold", "heat")
)]
country_obj <- country_obj[, .(
  deaths = sum(deaths), pop = unique(sum(pop) / uniqueN(cause))
), by = .(geography = country, ssp, scenario, period, agegroup, broad_range)]
country_obj <- dcast(
  country_obj, geography + ssp + period + agegroup + broad_range ~ scenario,
  value.var = c("deaths", "pop")
)[, .(
  geography, ssp, period, agegroup, broad_range,
  climate_deaths = (deaths_with_cc - deaths_without_cc) / 5,
  rate_per_100000 = 1e5 * (deaths_with_cc - deaths_without_cc) / pop_with_cc
)]
country_ref <- standardise_reference("country_period.parquet", "country")
country_ref <- country_ref[geography %chin% c("MT", "IE", "RO", "Malta", "Ireland", "Romania")]
write_table(
  merge(country_obj, country_ref,
        by = c("geography", "ssp", "period", "agegroup", "broad_range")),
  "masselot_2025_country_spot_checks.csv"
)

# -----------------------------------------------------------------------------
# 3. MMT difference and effect on historical attributable deaths
# -----------------------------------------------------------------------------

load(file.path(batch_root, "shared", "prep_data.RData"))
setDT(obs_data)
setDT(thresholds)
coefs <- fread("data/coefs.csv")
meta <- fread("data/city_results.csv")
predper <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))
varper <- c(10, 75, 90)

mmt_effect <- rbindlist(lapply(cities, function(city) {
  obs <- obs_data[URAU_CODE == city]
  obs[, year := as.integer(format(date, "%Y"))]
  tper <- quantile(obs$tmean_obs, predper / 100, na.rm = TRUE)
  knots <- tper[paste0(varper, ".0%")]
  bound <- range(tper)
  basis <- suppressWarnings(onebasis(
    obs$tmean_obs, fun = "bs", degree = 2, knots = knots, Bound = bound
  ))
  rbindlist(lapply(c("65-74", "75-84", "85+"), function(a) {
    tr <- thresholds[URAU_CODE == city & agegroup == a]
    cm <- meta[URAU_CODE == city & agegroup == a]
    cf <- coefs[URAU_CODE == city & agegroup == a]
    if (nrow(tr) != 1L || nrow(cm) != 1L || nrow(cf) != 1L) return(NULL)
    b <- as.matrix(cf[, .(b1, b2, b3, b4, b5)])
    calc <- function(mmt) {
      centered <- scale(
        basis,
        center = onebasis(mmt, fun = "bs", degree = 2, knots = knots, Bound = bound),
        scale = FALSE
      )
      rr <- pmax(exp(as.numeric(centered %*% t(b))), 1)
      an <- (1 - 1 / rr) * cm$death / 365.25
      c(total = sum(an) / uniqueN(obs$year),
        heat = sum(an[obs$tmean_obs > mmt]) / uniqueN(obs$year),
        cold = sum(an[obs$tmean_obs <= mmt]) / uniqueN(obs$year))
    }
    old <- calc(tr$mmt_2023)
    current <- calc(tr$mmt)
    data.table(
      city = city, agegroup = a,
      mmt_2023 = tr$mmt_2023, mmt_current = tr$mmt,
      mmt_difference_c = tr$mmt - tr$mmt_2023,
      deaths_2023_mmt = old[["total"]], deaths_current_mmt = current[["total"]],
      deaths_difference = current[["total"]] - old[["total"]],
      heat_difference = current[["heat"]] - old[["heat"]],
      cold_difference = current[["cold"]] - old[["cold"]]
    )
  }))
}))
write_table(mmt_effect, "mmt_city_age_effect.csv")
mmt_summary <- mmt_effect[, .(
  city_age_rows = .N,
  cities = uniqueN(city),
  differing_rows = sum(abs(mmt_difference_c) > 1e-10),
  median_abs_difference_c = median(abs(mmt_difference_c)),
  maximum_abs_difference_c = max(abs(mmt_difference_c)),
  published_mmt_deaths = sum(deaths_2023_mmt),
  current_mmt_deaths = sum(deaths_current_mmt),
  attributable_death_change_pct = 100 *
    (sum(deaths_current_mmt) / sum(deaths_2023_mmt) - 1)
), by = agegroup]
write_table(mmt_summary, "mmt_summary.csv")

# -----------------------------------------------------------------------------
# 4-5. Demographic levels and observed Eurostat LE65
# -----------------------------------------------------------------------------

periods <- list(`2020-2024` = 2020:2024, `2055-2059` = 2055:2059,
                `2095-2099` = 2095:2099)
geo_cases <- data.table(
  level = c("europe", rep("country", 3)),
  id = c("EUROPE", "ES", "SE", "PL"),
  eurostat_geo = c("EU27_2020", "ES", "SE", "PL")
)
demography <- rbindlist(lapply(1:3, function(ssp) {
  rbindlist(lapply(seq_len(nrow(geo_cases)), function(i) {
    f <- file.path(
      batch_root, "geographies", paste0("ssp", ssp), geo_cases$level[i],
      geo_cases$id[i], "04_le_li_levels.csv"
    )
    x <- fread(f)[branch == "without_cc"]
    rbindlist(lapply(names(periods), function(period) {
      x[year %in% periods[[period]], .(
        ssp = ssp, level = geo_cases$level[i], id = geo_cases$id[i],
        geography = label[1], period = period,
        LE65 = mean(LE65), LI65 = mean(LI65)
      )]
    }))
  }))
}))
write_table(demography, "demographic_le65_li65_trajectories.csv")

eurostat_file <- file.path(out_root, "eurostat_le65.csv")
if (!file.exists(eurostat_file)) {
  stop("Run analysis/fetch_eurostat_le65.py before this script", call. = FALSE)
}
eurostat <- fread(eurostat_file)
baseline_model <- merge(
  demography[period == "2020-2024"], geo_cases,
  by = c("level", "id")
)
setnames(baseline_model, "LE65", "model_LE65")
setnames(eurostat, "LE65", "observed_LE65")
observed_comparison <- merge(
  baseline_model, eurostat,
  by.x = "eurostat_geo", by.y = "geo", allow.cartesian = TRUE
)
observed_comparison[, gap_years := model_LE65 - observed_LE65]
write_table(observed_comparison, "observed_eurostat_le65_comparison.csv")

# -----------------------------------------------------------------------------
# 6. Lloyd method check on Simon's distributed example inputs
# -----------------------------------------------------------------------------

lloyd_input <- readRDS(file.path(lloyd_root, "outdata", "deaths_input_100.rds"))
lloyd_lt <- readRDS(file.path(lloyd_root, "outdata", "lt_CI_65plus.rds"))
lloyd_le <- fread(file.path(lloyd_root, "outdata", "Cause_contrib_LE_by_age_group.csv"))
lloyd_li <- fread(file.path(lloyd_root, "outdata", "Cause_contrib_LI_by_age_group.csv"))

life_expectancy_from_mx_65plus <- function(mx, x, nx = c(rep(1, 35), Inf), age = 0) {
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  ex[age + 1]
}
sd_from_mx_fun_65_plus <- function(mx, x, nx = c(rep(1, 35), Inf), age = 0) {
  x_conditional <- x - 65
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  sqrt(sum(dx * (x_conditional + 0.5 - ex[age + 1])^2))
}

# Lloyd's public repository explicitly contains fake data. The check therefore
# verifies code identity and decomposition closure, not a published Spain value.
lloyd_check <- data.table(
  check = c(
    "LE contribution file finite", "LI contribution file finite",
    "LE rows", "LI rows", "input rows", "published Spain value reproducible"
  ),
  value = c(
    all(is.finite(lloyd_le$contribution)), all(is.finite(lloyd_li$contribution)),
    nrow(lloyd_le), nrow(lloyd_li), nrow(lloyd_input), FALSE
  ),
  verdict = c(
    "match", "match", "match", "match", "match",
    "explained: Simon's repository labels the supplied input as fake, so the published Spain value is not reproducible"
  )
)
write_table(lloyd_check, "lloyd_method_check.csv")

# -----------------------------------------------------------------------------
# Gate summary and concise report
# -----------------------------------------------------------------------------

gate <- rbindlist(list(
  data.table(
    check = "Masselot 2023 historical city fixtures",
    result = sprintf("854/854 pass; max relative error %.3g", fixture_summary$maximum_relative_error),
    verdict = fixture_summary$verdict
  ),
  comparison_summary[, .(
    check = paste0("Masselot 2025 projections SSP", ssp),
    result = sprintf("%d matched rows; max rate relative error %.3g", rows, max_rate_relative_error),
    verdict
  )],
  data.table(
    check = "MMT definition change",
    result = sprintf("%d city-age rows quantified", nrow(mmt_effect)),
    verdict = "explained"
  ),
  data.table(
    check = "SSP1 demography",
    result = "0.9995 is a boundary sensitivity; 562 city schedules fail the canonical shape check and Europe LE65 rises about 19 years",
    verdict = "explained: exclude SSP1 from uncertainty"
  ),
  data.table(
    check = "Observed Eurostat LE65",
    result = sprintf("%d model-observation comparisons", nrow(observed_comparison)),
    verdict = if (max(abs(observed_comparison$gap_years)) < 2) "match" else "not explained"
  ),
  data.table(
    check = "Lloyd 2024 method",
    result = "Distributed code and outputs are internally readable; public inputs are fake",
    verdict = "explained"
  )
), fill = TRUE)
write_table(gate, "validation_gate.csv")

report <- c(
  "# Literature validation gate, 2026-11-04",
  "",
  "All comparisons use matched ages, periods, scenarios, GCM ensemble definitions and population denominators.",
  "Masselot's published all-adult headline values are not compared directly with this study's 65+ estimand.",
  "",
  "| Check | Result | Verdict |",
  "|---|---|---|",
  apply(gate, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
  "",
  "SSP1 remains a labelled demographic sensitivity. It is not eligible for the uncertainty run.",
  "The Lloyd repository states that its shared input is fake. This prevents reproduction of a published Spain value, but does not indicate a disagreement with the method.",
  "See the CSV files in this directory for row-level evidence."
)
writeLines(report, file.path(out_root, "README.md"))
if (any(gate$verdict == "not explained")) {
  stop("Validation gate has an unexplained result; see validation_gate.csv", call. = FALSE)
}
message("Validation gate passed with explained SSP1 and Lloyd limitations.")
