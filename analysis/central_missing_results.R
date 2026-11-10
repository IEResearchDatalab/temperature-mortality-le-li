#!/usr/bin/env Rscript

# Missing central results requested for the LE/LI paper. This script is
# presentation-only: it reads the completed pooled geography outputs and writes
# a new, dated result pack. It never changes the canonical batch.

suppressPackageStartupMessages(library(data.table))

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
out_root <- Sys.getenv("CENTRAL_EXTRAS_ROOT", "results/europe/central_extras_20261104")
table1_file <- Sys.getenv(
  "TABLE1_FILE",
  "results/figures_simon_20oct/table1/table1_full.csv"
)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

ages <- 65:100
nx <- c(rep(1, 35), Inf)
temp_causes <- c("ExtrCold", "ModCold", "ModHeat", "NormExtrHeat", "ExtrHeat")
heat_causes <- c("ModHeat", "NormExtrHeat", "ExtrHeat")
cold_causes <- c("ExtrCold", "ModCold")
all_causes <- c(temp_causes, "rest")
start_years <- 2020:2024
end_years <- 2095:2099

life_table <- function(mx) {
  stopifnot(length(mx) == length(ages), all(is.finite(mx)), all(mx >= 0))
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  li65 <- sqrt(sum(dx * (ages - 65 + 0.5 - ex[1])^2))
  list(LE65 = ex[1], LI65 = li65, ex = ex, lx = lx, dx = dx)
}

geo_index <- rbindlist(lapply(c("europe", "region", "country"), function(level) {
  dirs <- list.dirs(
    file.path(batch_root, "geographies", "ssp3", level),
    recursive = FALSE, full.names = TRUE
  )
  data.table(level = level, id = basename(dirs), directory = dirs)
}))
geo_index <- geo_index[file.exists(file.path(directory, "03_master_table.csv"))]
if (nrow(geo_index) != 35L) {
  stop("Expected 35 pooled geographies; found ", nrow(geo_index), call. = FALSE)
}

master_file <- function(ssp, level, id) file.path(
  batch_root, "geographies", paste0("ssp", ssp), level, id, "03_master_table.csv"
)
geo_file <- function(ssp, level, id, filename) file.path(
  batch_root, "geographies", paste0("ssp", ssp), level, id, filename
)

cause_deleted_one <- function(master, ssp, level, id, label) {
  x <- master[
    branch == "with_cc" & year %in% start_years &
      age %in% ages & range %in% all_causes,
    .(mx = deaths_component / pop),
    by = .(year, age, cause = range)
  ][, .(mx = mean(mx)), by = .(age, cause)]
  if (nrow(x) != length(ages) * length(all_causes)) {
    stop("Incomplete baseline cause grid for SSP", ssp, " ", level, "/", id)
  }
  base <- life_table(x[, sum(mx), by = age][order(age)]$V1)
  deletion_sets <- c(
    setNames(lapply(temp_causes, c), temp_causes),
    list(AllTemperature = temp_causes)
  )
  rbindlist(lapply(names(deletion_sets), function(deletion) {
    deleted_mx <- x[!cause %in% deletion_sets[[deletion]], sum(mx), by = age][order(age)]$V1
    deleted <- life_table(deleted_mx)
    data.table(
      ssp = ssp, level = level, id = id, geography = label,
      period = "2020-2024", deleted_cause = deletion,
      LE65_observed = base$LE65, LE65_cause_deleted = deleted$LE65,
      LE65_gain_if_deleted_years = deleted$LE65 - base$LE65,
      LI65_observed = base$LI65, LI65_cause_deleted = deleted$LI65,
      LI65_change_if_deleted = deleted$LI65 - base$LI65
    )
  }))
}

life_years_one <- function(master, ssp, level, id, label) {
  x <- master[
    year %between% c(2020L, 2099L) & age %in% ages & range %in% all_causes,
    .(branch, year, age, cause = range, deaths = deaths_component, pop)
  ]
  if (anyDuplicated(x, by = c("branch", "year", "age", "cause"))) {
    stop("Duplicate cumulative-burden rows for SSP", ssp, " ", level, "/", id)
  }
  ex_rows <- x[, {
    mx <- .SD[, sum(deaths / pop), by = age][order(age)]$V1
    data.table(age = ages, ex = life_table(mx)$ex)
  }, by = .(branch, year)]
  x <- merge(x, ex_rows, by = c("branch", "year", "age"))
  x[, life_years_lost := deaths * ex]
  denom <- unique(x[, .(year, age, pop)])[, sum(pop)]

  by_world <- x[cause %in% temp_causes, .(
    life_years_lost = sum(life_years_lost),
    attributable_deaths = sum(deaths)
  ), by = .(branch, cause)]
  totals <- by_world[, .(
    life_years_lost = sum(life_years_lost),
    attributable_deaths = sum(attributable_deaths)
  ), by = branch][, cause := "AllTemperature"]
  by_world <- rbind(by_world, totals, use.names = TRUE)

  wide <- dcast(
    by_world, cause ~ branch,
    value.var = c("life_years_lost", "attributable_deaths")
  )
  diff_worlds <- wide[, .(
    cause,
    life_years_lost = life_years_lost_with_cc - life_years_lost_without_cc,
    attributable_deaths = attributable_deaths_with_cc - attributable_deaths_without_cc
  )][, estimand := "climate_increment_difference_of_worlds"]

  # A fixed-ex counterfactual isolates the change in attributable deaths. It
  # weights both worlds with remaining life expectancy in the without-CC world.
  fixed <- merge(
    x[branch == "with_cc" & cause %in% temp_causes,
      .(year, age, cause, deaths_with = deaths)],
    x[branch == "without_cc" & cause %in% temp_causes,
      .(year, age, cause, deaths_without = deaths, ex_without = ex)],
    by = c("year", "age", "cause")
  )
  fixed[, life_years_lost := (deaths_with - deaths_without) * ex_without]
  fixed_result <- fixed[, .(
    life_years_lost = sum(life_years_lost),
    attributable_deaths = sum(deaths_with - deaths_without)
  ), by = cause]
  fixed_result <- rbind(
    fixed_result,
    fixed_result[, .(
      cause = "AllTemperature",
      life_years_lost = sum(life_years_lost),
      attributable_deaths = sum(attributable_deaths)
    )]
  )[, estimand := "climate_increment_fixed_without_cc_ex"]

  world_result <- rbindlist(list(
    by_world[branch == "with_cc", .(
      cause, life_years_lost, attributable_deaths,
      estimand = "whole_burden_with_cc"
    )],
    by_world[branch == "without_cc", .(
      cause, life_years_lost, attributable_deaths,
      estimand = "counterfactual_burden_without_cc"
    )],
    diff_worlds,
    fixed_result
  ), use.names = TRUE)
  world_result[, `:=`(
    ssp = ssp, level = level, id = id, geography = label,
    years = "2020-2099", population_person_years_65plus = denom,
    life_years_lost_per_100000_65plus = 1e5 * life_years_lost / denom,
    attributable_deaths_per_100000_65plus = 1e5 * attributable_deaths / denom
  )]
  setcolorder(world_result, c(
    "ssp", "level", "id", "geography", "years", "estimand", "cause",
    "life_years_lost", "life_years_lost_per_100000_65plus",
    "attributable_deaths", "attributable_deaths_per_100000_65plus",
    "population_person_years_65plus"
  ))
  world_result
}

tier_shares_one <- function(master, ssp, level, id, label) {
  rbindlist(lapply(c(65L, 85L), function(min_age) {
    x <- master[
      year %in% start_years & age >= min_age & range %in% temp_causes,
      .(deaths = sum(an)), by = .(branch, cause = range)
    ]
    with <- x[branch == "with_cc"]
    total <- with[, sum(deaths)]
    without_total <- x[branch == "without_cc", sum(deaths)]
    share <- function(causes) 100 * with[cause %in% causes, sum(deaths)] / total
    data.table(
      ssp = ssp, level = level, id = id, geography = label,
      period = "2020-2024", age_scope = paste0(min_age, "+"),
      temperature_deaths_with_cc = total,
      persistent_share_pct = 100 * without_total / total,
      climate_change_share_pct = 100 * (total - without_total) / total,
      general_tier_pct = share(c("ModCold", "ModHeat")),
      sustained_tier_pct = share("NormExtrHeat"),
      episodic_tier_pct = share(c("ExtrCold", "ExtrHeat"))
    )
  }))
}

age_blocks <- list(
  `65-74` = 65:74,
  `75-84` = 75:84,
  `85+` = 85:100,
  `85-94` = 85:94,
  `95+` = 95:100
)

age_erosion_one <- function(ssp, level, id, label) {
  within <- fread(geo_file(
    ssp, level, id, "04_within_branch_period_decomposition.csv"
  ))
  between <- fread(geo_file(
    ssp, level, id, "04_between_branch_decomposition.csv"
  ))
  # Summing consecutive-period Horiuchi contributions defines the age-block
  # contribution to the full 2020-2024 to 2095-2099 gain along the actual path.
  gain <- within[
    period_from >= "2020-2024" & period_to <= "2095-2099",
    .(gain_years = sum(le_contribution)), by = .(branch, age)
  ]
  gain <- dcast(gain, age ~ branch, value.var = "gain_years")
  start <- between[period == "2020-2024"]
  end <- between[period == "2095-2099"]
  cross <- merge(
    start[, .(age, cause, contribution_start = le_contribution)],
    end[, .(age, cause, contribution_end = le_contribution)],
    by = c("age", "cause")
  )
  cross[, loss_years := contribution_start - contribution_end]

  rbindlist(lapply(names(age_blocks), function(block) {
    aa <- age_blocks[[block]]
    g0 <- gain[age %in% aa, sum(without_cc)]
    g1 <- gain[age %in% aa, sum(with_cc)]
    heat_loss <- cross[age %in% aa & cause %in% heat_causes, sum(loss_years)]
    cross_all_loss <- cross[age %in% aa, sum(loss_years)]
    all_loss <- g0 - g1
    stable <- is.finite(g0) && abs(12 * g0) >= 0.1
    data.table(
      ssp = ssp, level = level, id = id, geography = label,
      period_from = "2020-2024", period_to = "2095-2099",
      age_group = block,
      demographic_gain_without_cc_months = 12 * g0,
      gain_with_cc_months = 12 * g1,
      heat_erosion_months = 12 * heat_loss,
      all_temperature_erosion_months = 12 * all_loss,
      heat_erosion_ratio_pct = 100 * heat_loss / g0,
      all_temperature_erosion_ratio_pct = 100 * all_loss / g0,
      cross_section_all_loss_months = 12 * cross_all_loss,
      path_vs_cross_section_gap_months = 12 * (all_loss - cross_all_loss),
      stable_denominator = stable
    )
  }))
}

cause_deleted <- list()
life_years <- list()
tiers <- list()
erosion <- list()
counter <- 0L
for (ssp in 1:3) {
  for (i in seq_len(nrow(geo_index))) {
    counter <- counter + 1L
    level <- geo_index$level[i]
    id <- geo_index$id[i]
    f <- master_file(ssp, level, id)
    if (!file.exists(f)) stop("Missing master table: ", f, call. = FALSE)
    master <- fread(
      f,
      select = c("label", "branch", "year", "age", "range", "pop", "an", "deaths_component")
    )
    label <- master$label[1]
    cause_deleted[[counter]] <- cause_deleted_one(master, ssp, level, id, label)
    life_years[[counter]] <- life_years_one(master, ssp, level, id, label)
    tiers[[counter]] <- tier_shares_one(master, ssp, level, id, label)
    erosion[[counter]] <- age_erosion_one(ssp, level, id, label)
    message(sprintf("[%03d/105] SSP%d %s/%s", counter, ssp, level, id))
  }
}

cause_deleted <- rbindlist(cause_deleted)
life_years <- rbindlist(life_years)
tiers <- rbindlist(tiers)
erosion <- rbindlist(erosion)

fwrite(cause_deleted, file.path(out_root, "cause_deleted_2020_2024.csv"))
fwrite(life_years, file.path(out_root, "life_years_lost_2020_2099.csv"))
fwrite(tiers, file.path(out_root, "baseline_tier_shares_2020_2024.csv"))
fwrite(erosion, file.path(out_root, "age_specific_erosion_2020_2099.csv"))

table1 <- fread(table1_file)
consistency <- table1[, .(
  ssp, level, id, geography,
  heat_ratio_ensemble_pct = 100 * heat_loss_ensemble_months / g0_months,
  heat_ratio_gcm_median_pct = heat_ratio_median_pct,
  heat_median_minus_ensemble_pp = heat_ratio_median_pct -
    100 * heat_loss_ensemble_months / g0_months,
  all_ratio_ensemble_pct = all_ratio_pct,
  all_ratio_gcm_median_pct = all_ratio_gcm_median_pct,
  all_median_minus_ensemble_pp = all_ratio_gcm_median_pct - all_ratio_pct
)]
fwrite(consistency, file.path(out_root, "erosion_ensemble_vs_gcm_median.csv"))

checks <- data.table(
  check = c(
    "geographies", "cause_deleted_rows", "life_year_rows", "tier_rows",
    "age_erosion_rows", "tier_closure_max_abs_pp", "cause_deleted_le_nonnegative",
    "all_temperature_equals_sum_causes_yll_max_abs"
  ),
  value = c(
    nrow(geo_index), nrow(cause_deleted), nrow(life_years), nrow(tiers),
    nrow(erosion),
    max(abs(tiers$general_tier_pct + tiers$sustained_tier_pct +
              tiers$episodic_tier_pct - 100)),
    min(cause_deleted$LE65_gain_if_deleted_years),
    life_years[cause != "AllTemperature", .(v = sum(life_years_lost)),
               by = .(ssp, level, id, estimand)][
      life_years[cause == "AllTemperature"],
      max(abs(v - i.life_years_lost)),
      on = .(ssp, level, id, estimand)
    ]
  ),
  expected = c(
    "35", "630", "2520", "210", "525", "<=1e-9", ">=0", "<=1e-6"
  )
)
fwrite(checks, file.path(out_root, "checks.csv"))
if (checks[check == "tier_closure_max_abs_pp", as.numeric(value)] > 1e-9 ||
    checks[check == "cause_deleted_le_nonnegative", as.numeric(value)] < -1e-10 ||
    checks[check == "all_temperature_equals_sum_causes_yll_max_abs", as.numeric(value)] > 1e-6) {
  stop("One or more central-extra checks failed; see checks.csv", call. = FALSE)
}

readme <- c(
  "# Central missing results, 2026-11-04",
  "",
  "Source: completed pooled ENSEMBLE geography outputs under `heatsplit_20261013`.",
  "The canonical batch was read only.",
  "",
  "## Estimands",
  "",
  "- Cause-deleted LE65 is the gain after setting the selected cause-specific mortality rate to zero while retaining every competing risk. LI65+ reports the signed change after deletion.",
  "- Life-years lost equal attributable deaths at age x multiplied by period remaining life expectancy at x. The pack reports both climate worlds, their difference, and a fixed-without-CC-ex increment that isolates changed attributable deaths.",
  "- An age block's gain is the sum of its single-age Horiuchi contributions over every consecutive five-year transition from 2020-2024 to 2095-2099. Its all-temperature erosion ratio is (gain without CC - gain with CC) / gain without CC. Heat erosion uses the change from start to end in the heat components of the between-world Horiuchi decomposition, divided by the same denominator.",
  "- The 95+ age block is retained but marked unstable whenever the absolute denominator is below 0.1 month.",
  "- Baseline tier shares use attributable deaths in 2020-2024. General is moderate cold plus moderate heat. Sustained is normal extreme heat. Episodic is extreme cold plus extreme heat. Persistent is the without-CC total divided by the with-CC total.",
  "- The consistency table reports heat and all-temperature erosion using both the pooled ENSEMBLE calculation and the median of 19 per-GCM ratios.",
  "",
  "All rates per 100,000 use summed population person-years at ages 65+ over 2020-2099."
)
writeLines(readme, file.path(out_root, "README.md"))
message("Wrote central missing results to ", out_root)
