#!/usr/bin/env Rscript

# Pooled ERF-draw uncertainty for one SSP x GCM. The same draw IDs are used in
# every city, and city counts are pooled before life tables are constructed.

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(dplyr)
  library(dlnm)
  library(splines)
})

ssp <- as.integer(Sys.getenv("SSP", "3"))
ssp_text <- as.character(ssp)
gcm <- Sys.getenv("GCM", "GFDL_ESM4")
batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
uncertainty_root <- Sys.getenv("UNCERTAINTY_ROOT", "results/europe/uncertainty_20261104")
draw_text <- Sys.getenv("DRAW_IDS", "1-500")
city_text <- Sys.getenv("CITY_IDS", "all")

parse_ids <- function(x) {
  if (grepl("^[0-9]+-[0-9]+$", x)) {
    z <- as.integer(strsplit(x, "-", fixed = TRUE)[[1]])
    return(seq.int(z[1], z[2]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}
draw_ids <- parse_ids(draw_text)
if (anyNA(draw_ids) || anyDuplicated(draw_ids) || any(!draw_ids %in% 1:500)) {
  stop("DRAW_IDS must be unique IDs from 1 to 500", call. = FALSE)
}
draw_ids <- sort(draw_ids)

meta <- unique(fread("data/city_results.csv", select = c("URAU_CODE", "region")))
selected_cities <- if (city_text == "all") sort(meta$URAU_CODE) else sort(strsplit(city_text, ",", fixed = TRUE)[[1]])
if (any(!selected_cities %chin% meta$URAU_CODE)) stop("Unknown CITY_IDS", call. = FALSE)
meta <- meta[URAU_CODE %chin% selected_cities]
geographies <- c("Europe", paste(sort(unique(meta$region)), "Europe"))

out_dir <- file.path(uncertainty_root, "part01", paste0("ssp", ssp), gcm)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
done_file <- file.path(out_dir, ".done")
if (file.exists(done_file) && city_text == "all" && identical(draw_ids, 1:500)) {
  message("Already complete: ", out_dir)
  quit(save = "no")
}

years <- 2020:2099
endpoint_years <- c(2020:2024, 2095:2099)
agegroups <- c("65-74", "75-84", "85+")
ages <- 65:100
ranges <- c("ExtrCold", "ModCold", "ModHeat", "NormExtrHeat", "ExtrHeat")
branches <- c("with_cc", "without_cc")
heat_ranges <- c("ModHeat", "NormExtrHeat", "ExtrHeat")
cold_ranges <- c("ExtrCold", "ModCold")
nx <- c(rep(1, 35), Inf)
predper <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))
varper <- c(10, 75, 90)
projrange <- c(2015, seq(2030, 2100, by = 10))
histrange <- c(2000, 2014)

D <- length(draw_ids); G <- length(geographies); B <- 2L; Y <- 80L
A <- 3L; R <- 5L; E <- 10L; X <- 36L
dim_group <- c(D, G, B, Y, A, R)
dim_single <- c(D, G, B, E, X, R)
pool_group <- array(0, dim = dim_group)
pool_single <- array(0, dim = dim_single)
dimnames(pool_group) <- list(
  draw = draw_ids, geography = geographies, branch = branches, year = years,
  agegroup = agegroups, range = ranges
)
dimnames(pool_single) <- list(
  draw = draw_ids, geography = geographies, branch = branches,
  year = endpoint_years, age = ages, range = ranges
)
pool_pop <- array(0, dim = c(G, E, X))
pool_death <- array(0, dim = c(G, E, X))
pool_pop_scope <- array(0, dim = c(G, Y, 3L),
                        dimnames = list(geographies, years, c("65+", "85+", "95+")))
pool_95 <- array(0, dim = c(D, G, B, Y, R))

load(file.path(batch_root, "shared", "prep_data.RData"))
setDT(obs_data); setDT(thresholds)
draw_ds <- open_dataset("references/2025-masselot-data/coef_simu.gz.parquet")
temp_ds <- open_dataset("data/tmeanproj.gz.parquet")

isimip3 <- function(obshist, simhist, simfut, yearobshist, yearsimhist, yearsimfut, detrend = TRUE) {
  if (detrend) {
    obstrend <- lm(obshist ~ yearobshist, na.action = na.exclude) |> predict() |> scale(scale = FALSE)
    simhisttrend <- lm(simhist ~ yearsimhist, na.action = na.exclude) |> predict() |> scale(scale = FALSE)
    simfuttrend <- lm(simfut ~ yearsimfut, na.action = na.exclude) |> predict() |> scale(scale = FALSE)
    obshist <- obshist - obstrend
    simhist <- simhist - simhisttrend
    simfut <- simfut - simfuttrend
  }
  ecdfobs <- ecdf(obshist)(obshist)
  deltaadd <- quantile(simfut, ecdfobs, na.rm = TRUE) - quantile(simhist, ecdfobs, na.rm = TRUE)
  obsfut <- deltaadd + obshist
  simfutcdf <- pnorm(simfut, mean(simfut, na.rm = TRUE), sd(simfut, na.rm = TRUE))
  out <- qnorm(p = simfutcdf, mean = mean(obsfut, na.rm = TRUE), sd = sd(obsfut, na.rm = TRUE))
  if (detrend) out <- out + simfuttrend
  as.numeric(out)
}

life_table <- function(mx) {
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  li <- sqrt(sum(dx * (ages - 65 + 0.5 - ex[1])^2))
  c(LE65 = ex[1], LI65 = li)
}

lifetime_risk <- function(mx_cause, mx_rest) {
  total <- rowSums(mx_cause) + mx_rest
  q <- 1 - exp(-total * nx)
  l <- head(cumprod(c(1, 1 - q)), -1)
  fractions <- mx_cause / total
  fractions[!is.finite(fractions)] <- 0
  d <- fractions * (l * q)
  c(all = 100 * sum(d), heat = 100 * sum(d[, heat_ranges, drop = FALSE]))
}

add_geo <- function(city_region) match(c("Europe", paste(city_region, "Europe")), geographies)

for (ci in seq_along(selected_cities)) {
  city <- selected_cities[ci]
  region <- meta[URAU_CODE == city, region][1]
  geo_idx <- add_geo(region)
  if (anyNA(geo_idx)) {
    stop("Geography lookup failed for ", city, ": region=", region,
         "; available=", paste(geographies, collapse = ","), call. = FALSE)
  }
  obs <- copy(obs_data[URAU_CODE == city])
  obs[, `:=`(
    year = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    month_day = format(date, "%m-%d")
  )]
  obs_hist <- obs[year %between% histrange & month_day != "02-29"]
  tper <- quantile(obs$tmean_obs, predper / 100, na.rm = TRUE)
  knots <- tper[paste0(varper, ".0%")]
  bound <- range(tper)

  tcol <- paste0("tas_", gcm)
  tmean <- temp_ds |>
    filter(URAU_CODE == city, ssp %in% c("hist", ssp_text)) |>
    select(date, ssp, all_of(tcol)) |>
    collect() |>
    as.data.table()
  setnames(tmean, tcol, "tmean")
  tmean[, `:=`(
    year = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    month_day = format(date, "%m-%d")
  )]
  tmean <- tmean[month_day != "02-29"]
  if (gcm == "IITM_ESM" && ssp == 3L) {
    fill <- tmean[ssp == ssp_text & year == 2098L, .(month_day, fill = tmean)]
    tmean <- merge(tmean, fill, by = "month_day", all.x = TRUE, sort = FALSE)
    tmean[ssp == ssp_text & year == 2099L, tmean := fill]
    tmean[, fill := NULL]
    setorder(tmean, date)
  }
  hist_sim <- tmean[ssp == "hist" & year %between% histrange]
  cal <- rbind(
    tmean[ssp == "hist" & year %between% histrange],
    tmean[ssp == ssp_text & year >= 2015L]
  )
  cal[, calperiod := cut(year, c(histrange[1], projrange), right = FALSE)]
  cal[, full := {
    obs_m <- obs_hist[month == .BY$month]
    sim_m <- hist_sim[month == .BY$month]
    isimip3(obs_m$tmean_obs, sim_m$tmean, tmean,
            obs_m$year, sim_m$year, year)
  }, by = .(month, calperiod)]
  cal[, year5 := (year %/% 5L) * 5L]
  cal[, demo := {
    ref <- cal[month == .BY$month & year5 == 2010L]
    isimip3(ref$full, full, full, ref$year, year, year)
  }, by = .(month, year5)]
  static_p97 <- unique(thresholds[URAU_CODE == city, p97_5])
  model_base <- quantile(cal[year %between% histrange, full], 0.975, na.rm = TRUE)
  moving <- rbindlist(lapply(2015:2099, function(target) {
    w <- cal[year %between% c(target - 15L, target - 1L), full]
    data.table(year = target, moving = max(
      static_p97,
      static_p97 + quantile(w, 0.975, na.rm = TRUE) - model_base
    ))
  }))

  dem <- fread(file.path(
    batch_root, paste0("ssp", ssp), city, "demography", "00_demography_single_age.csv"
  ))[year %in% years & age %in% ages]
  dem_group <- dem[, .(death = sum(death)), by = .(year, agegroup)]
  dem[, group_death := sum(death), by = .(year, agegroup)]
  dem[, weight := fifelse(group_death > 0, death / group_death, 0)]
  for (gi in geo_idx) {
    e <- dem[year %in% endpoint_years]
    pool_pop[gi, , ] <- pool_pop[gi, , ] +
      matrix(e[order(match(year, endpoint_years), age), pop], nrow = E, ncol = X, byrow = TRUE)
    pool_death[gi, , ] <- pool_death[gi, , ] +
      matrix(e[order(match(year, endpoint_years), age), death], nrow = E, ncol = X, byrow = TRUE)
    scope <- dem[, .(
      pop65 = sum(pop), pop85 = sum(pop[age >= 85]), pop95 = sum(pop[age >= 95])
    ), by = year][order(year)]
    pool_pop_scope[gi, , 1] <- pool_pop_scope[gi, , 1] + scope$pop65
    pool_pop_scope[gi, , 2] <- pool_pop_scope[gi, , 2] + scope$pop85
    pool_pop_scope[gi, , 3] <- pool_pop_scope[gi, , 3] + scope$pop95
  }

  coef_city <- draw_ds |>
    filter(URAU_CODE == city, sim %in% draw_ids, agegroup %in% agegroups) |>
    collect() |>
    as.data.table()
  if (nrow(coef_city) != D * A || uniqueN(coef_city$sim) != D) {
    stop("Incomplete coefficient draws for ", city, call. = FALSE)
  }

  branch_temp <- list(
    with_cc = merge(cal[year %in% years, .(year, temp = full)], moving, by = "year"),
    without_cc = cal[year %in% years, .(year, temp = demo, moving = static_p97)]
  )
  for (bi in seq_along(branches)) {
    branch <- branches[bi]
    daily <- branch_temp[[branch]]
    for (ai in seq_along(agegroups)) {
      ag <- agegroups[ai]
      tr <- thresholds[URAU_CODE == city & agegroup == ag]
      co <- coef_city[agegroup == ag][order(match(sim, draw_ids))]
      coef_matrix <- as.matrix(co[, .(b1, b2, b3, b4, b5)])
      basis <- onebasis(daily$temp, fun = "bs", degree = 2, knots = knots, Bound = bound)
      center <- onebasis(tr$mmt, fun = "bs", degree = 2, knots = knots, Bound = bound)
      log_rr <- scale(basis, center = center, scale = FALSE) %*% t(coef_matrix)
      af <- 1 - exp(-pmax(log_rr, 0))
      deaths <- dem_group[agegroup == ag][match(daily$year, year), death]
      weighted <- af * (deaths / 365)
      heat_idx <- if (tr$mmt > tr$p97_5) {
        rep("ExtrHeat", nrow(daily))
      } else {
        fifelse(
          daily$temp < tr$p97_5, "ModHeat",
          fifelse(daily$moving <= tr$p97_5 | daily$temp >= daily$moving,
                  "ExtrHeat", "NormExtrHeat")
        )
      }
      range_idx <- fifelse(
        daily$temp < tr$mmt,
        fifelse(daily$temp < tr$p2_5, "ExtrCold", "ModCold"),
        heat_idx
      )
      group_id <- (match(daily$year, years) - 1L) * R + match(range_idx, ranges)
      summed <- rowsum(weighted, group_id, reorder = FALSE)
      full <- matrix(0, Y * R, D)
      full[as.integer(rownames(summed)), ] <- summed
      arr <- array(t(full), c(D, Y, R))
      for (gi in geo_idx) {
        pool_group[, gi, bi, , ai, ] <- pool_group[, gi, bi, , ai, ] + arr
      }

      if (ag == "85+") {
        share95 <- dem[agegroup == "85+", .(
          share95 = sum(death[age >= 95]) / sum(death)
        ), by = year][order(year)]$share95
        arr95 <- sweep(arr, 2, share95, "*")
        for (gi in geo_idx) {
          pool_95[, gi, bi, , ] <- pool_95[, gi, bi, , ] + arr95
        }
      }

      ages_ag <- sort(unique(dem[agegroup == ag, age]))
      for (ei in seq_along(endpoint_years)) {
        yy <- endpoint_years[ei]
        yi <- match(yy, years)
        weights <- dem[year == yy & agegroup == ag][match(ages_ag, age), weight]
        for (xx in seq_along(ages_ag)) {
          xi <- match(ages_ag[xx], ages)
          add <- arr[, yi, ] * weights[xx]
          for (gi in geo_idx) {
            pool_single[, gi, bi, ei, xi, ] <-
              pool_single[, gi, bi, ei, xi, ] + add
          }
        }
      }
      rm(log_rr, af, weighted, summed, full, arr)
      gc(FALSE)
    }
  }
  message(sprintf("[%d/%d] %s", ci, length(selected_cities), city))
}

# Save pooled grouped counts with an explicit experiment dimension.
pooled_dt <- as.data.table(as.table(pool_group))
setnames(pooled_dt, "N", "an")
pooled_dt[, `:=`(ssp = ssp, gcm = gcm)]
write_parquet(pooled_dt, file.path(out_dir, "pooled_counts.parquet"), compression = "zstd")
rm(pooled_dt)
gc()

# Life-table and Table 1 metrics for each experiment.
metric_rows <- vector("list", D * G)
k <- 0L
max_rest_relative_roundoff <- 0
for (di in seq_len(D)) {
  for (gi in seq_len(G)) {
    k <- k + 1L
    schedules <- list()
    for (bi in seq_len(B)) {
      temp <- pool_single[di, gi, bi, , , , drop = TRUE]
      temp_without <- pool_single[di, gi, 2L, , , , drop = TRUE]
      rest <- pool_death[gi, , ] - apply(temp_without, c(1, 2), sum)
      rest_floor <- min(rest)
      rest_scale <- max(pool_death[gi, , ], 1)
      relative_roundoff <- max(0, -rest_floor / rest_scale)
      max_rest_relative_roundoff <- max(max_rest_relative_roundoff, relative_roundoff)
      if (relative_roundoff > 1e-8) {
        stop(sprintf(
          "Attributable mortality exceeds all-cause mortality: min rest=%g, relative excess=%g",
          rest_floor, relative_roundoff
        ))
      }
      rest <- pmax(rest, 0)
      make_period <- function(indices) {
        temp_rates <- sweep(
          temp[indices, , , drop = FALSE],
          c(1, 2), pool_pop[gi, indices, ], "/"
        )
        mx_temp <- apply(temp_rates, c(2, 3), mean)
        mx_rest <- apply(rest[indices, ] / pool_pop[gi, indices, ], 2, mean)
        list(temp = mx_temp, rest = mx_rest, total = rowSums(mx_temp) + mx_rest)
      }
      schedules[[branch <- branches[bi]]] <- list(
        start = make_period(1:5), end = make_period(6:10)
      )
    }
    stats <- lapply(schedules, function(z) lapply(z, function(p) life_table(p$total)))
    g0 <- stats$without_cc$end[1] - stats$without_cc$start[1]
    g1 <- stats$with_cc$end[1] - stats$with_cc$start[1]
    noheat_gain <- sapply(schedules, function(z) {
      life_table(z$end$rest + rowSums(z$end$temp[, cold_ranges, drop = FALSE]))[1] -
        life_table(z$start$rest + rowSums(z$start$temp[, cold_ranges, drop = FALSE]))[1]
    })
    all_loss <- g0 - g1
    cold_only_loss <- noheat_gain[2] - noheat_gain[1]
    heat_loss <- all_loss - cold_only_loss
    risk <- lifetime_risk(schedules$with_cc$end$temp, schedules$with_cc$end$rest)

    end_idx <- match(2095:2099, years)
    with_counts <- apply(pool_group[di, gi, 1L, end_idx, , , drop = FALSE], 6, sum)
    without_counts <- apply(pool_group[di, gi, 2L, end_idx, , , drop = FALSE], 6, sum)
    with_total <- sum(with_counts)
    metric_rows[[k]] <- data.table(
      ssp = ssp, gcm = gcm, draw = draw_ids[di], geography = geographies[gi],
      lifetime_risk_all_pct = risk[["all"]],
      lifetime_risk_heat_pct = risk[["heat"]],
      persistent_share_pct = 100 * sum(without_counts) / with_total,
      general_tier_pct = 100 * sum(with_counts[c("ModCold", "ModHeat")]) / with_total,
      sustained_tier_pct = 100 * with_counts[["NormExtrHeat"]] / with_total,
      episodic_tier_pct = 100 * sum(with_counts[c("ExtrCold", "ExtrHeat")]) / with_total,
      temperature_death_rate_65plus_end = 1e5 * with_total /
        sum(pool_pop_scope[gi, end_idx, 1]),
      climate_death_rate_65plus_end = 1e5 * sum(with_counts - without_counts) /
        sum(pool_pop_scope[gi, end_idx, 1]),
      g0_years = g0, g1_years = g1,
      heat_erosion_years = heat_loss, all_erosion_years = all_loss,
      heat_erosion_ratio_pct = 100 * heat_loss / g0,
      all_erosion_ratio_pct = 100 * all_loss / g0,
      delta_LE65_end = stats$with_cc$end[1] - stats$without_cc$end[1],
      delta_LI65_end = stats$with_cc$end[2] - stats$without_cc$end[2]
    )
  }
}
metrics <- rbindlist(metric_rows)
write_parquet(metrics, file.path(out_dir, "experiment_metrics.parquet"), compression = "zstd")

# Annual all-temperature risks at 85+ and 95+, absolute and climate increment.
risk_rows <- vector("list", G * 2L)
k <- 0L
for (gi in seq_len(G)) {
  for (si in 1:2) {
    k <- k + 1L
    source <- if (si == 1L) pool_group[, gi, , , 3L, , drop = TRUE] else pool_95[, gi, , , , drop = TRUE]
    denom <- pool_pop_scope[gi, , si + 1L]
    with <- apply(source[, 1L, , , drop = FALSE], c(1, 3), sum)
    without <- apply(source[, 2L, , , drop = FALSE], c(1, 3), sum)
    risk_table <- function(values, estimand_name) {
      z <- as.data.table(as.table(values))
      setnames(z, names(z)[1:2], c("draw_index", "year_index"))
      setnames(z, names(z)[3], "risk_per_100000")
      z[, .(
        draw = draw_ids[as.integer(factor(draw_index, levels = unique(draw_index)))],
        year = years[as.integer(factor(year_index, levels = unique(year_index)))],
        risk_per_100000,
        age_scope = c("85+", "95+")[si],
        estimand = estimand_name
      )]
    }
    risk_rows[[k]] <- rbindlist(list(
      risk_table(1e5 * with / rep(denom, each = D), "with_cc"),
      risk_table(1e5 * (with - without) / rep(denom, each = D), "climate_increment")
    ))[, `:=`(ssp = ssp, gcm = gcm, geography = geographies[gi])]
  }
}
annual_risk <- rbindlist(risk_rows)
write_parquet(annual_risk, file.path(out_dir, "annual_risk.parquet"), compression = "zstd")

checks <- data.table(
  check = c(
    "draws", "geographies", "tier_closure_max_abs_pp", "finite_metrics",
    "rest_nonnegative", "rest_roundoff_clamped_max_relative"
  ),
  value = c(
    as.character(uniqueN(metrics$draw)),
    as.character(uniqueN(metrics$geography)),
    as.character(max(abs(metrics$general_tier_pct + metrics$sustained_tier_pct +
                           metrics$episodic_tier_pct - 100))),
    if (all(vapply(metrics, function(x) !is.numeric(x) || all(is.finite(x)), logical(1))))
      "TRUE" else "FALSE",
    "TRUE",
    as.character(max_rest_relative_roundoff)
  ),
  expected = c(as.character(D), as.character(G), "<=1e-8", "TRUE", "TRUE", "<=1e-8")
)
fwrite(checks, file.path(out_dir, "checks.csv"))
if (as.numeric(checks[check == "tier_closure_max_abs_pp", value]) > 1e-8 ||
    checks[check == "finite_metrics", value] != "TRUE") {
  stop("Uncertainty checks failed", call. = FALSE)
}

manifest <- c(
  paste0("ssp=", ssp), paste0("gcm=", gcm),
  paste0("draw_ids=", paste(draw_ids, collapse = ",")),
  paste0("cities=", length(selected_cities)), paste0("geographies=", paste(geographies, collapse = ",")),
  paste0("coefficient_md5=", unname(tools::md5sum("references/2025-masselot-data/coef_simu.gz.parquet"))),
  paste0("git_commit=", system("git rev-parse HEAD", intern = TRUE)),
  paste0("finished=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
writeLines(manifest, file.path(out_dir, "manifest.txt"))
if (city_text == "all" && identical(draw_ids, 1:500)) writeLines("complete", done_file)
message("Saved uncertainty experiment outputs to ", out_dir)
