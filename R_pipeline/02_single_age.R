#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality / life expectancy pipeline -- Madrid pilot
#
# R Pipeline Step 02: Allocate grouped attributable numbers to single ages
#   Allocates Step 01 grouped attributable numbers to ages 65:100 using the
#   canonical single-age all-cause death schedule produced by Step 00.
#
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

message("\n[02] Allocating Madrid grouped AN to single ages...")

out_dir <- "results/phase1_madrid"
check_dir <- "results/checks"
fig_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

single_file <- file.path(out_dir, "02_single_age_an.csv")
checks_file <- file.path(check_dir, "02_single_age_checks.csv")
failures_file <- file.path(check_dir, "02_single_age_failures.csv")
fig_file <- file.path(fig_dir, "02_single_age_diagnostic.png")

city_id <- "ES001C"
city_name <- "Madrid"
future_years <- 2020:2099
age_groups <- c("65-74", "75-84", "85+")
age_slices <- list(
  "65-74" = 65:74,
  "75-84" = 75:84,
  "85+" = 85:100
)
#----- Load Step 00/01 outputs and restrict to the Madrid/SSP3 domain

grouped_dem <- fread(file.path(out_dir, "00_demography_grouped.csv"))
single_dem <- fread(file.path(out_dir, "00_demography_single_age.csv"))
grouped_an <- fread(file.path(out_dir, "01_attribution_grouped.csv"))

grouped_dem <- grouped_dem[geo_id == city_id & ssp == 3]
single_dem <- single_dem[geo_id == city_id & ssp == 3 & age %in% 65:100]
grouped_an <- grouped_an[geo_id == city_id & ssp == 3]

if (!nrow(grouped_dem)) stop("Grouped demography for Madrid is missing; run 00_demography.R first.", call. = FALSE)
if (!nrow(single_dem)) stop("Single-age demography for Madrid is missing; run 00_demography.R first.", call. = FALSE)
if (!nrow(grouped_an)) stop("Grouped AN for Madrid is missing; run 01_attribution.R first.", call. = FALSE)

grouped_dem <- grouped_dem[agegroup %in% age_groups]
grouped_an <- grouped_an[agegroup %in% age_groups]

if (any(!is.finite(grouped_dem$death)) || any(grouped_dem$death < 0)) {
  stop("Invalid grouped demographic deaths detected.", call. = FALSE)
}
if (any(!is.finite(single_dem$death)) || any(single_dem$death < 0)) {
  stop("Invalid single-age demographic deaths detected.", call. = FALSE)
}
if (any(!is.finite(grouped_an$an))) {
  stop("Grouped AN contains non-finite values.", call. = FALSE)
}

full_dem_grid <- CJ(year = future_years, agegroup = age_groups)
full_an_grid <- CJ(year = future_years, agegroup = age_groups, range = c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat"), branch = c("with_cc", "without_cc"))

dem_keys <- grouped_dem[, .(year, agegroup)]
an_keys <- grouped_an[, .(year, agegroup, range, branch)]
if (nrow(fsetdiff(full_dem_grid, unique(dem_keys))) || nrow(fsetdiff(unique(dem_keys), full_dem_grid)) || anyDuplicated(dem_keys)) {
  stop("Grouped demographic domain is incomplete for Madrid.", call. = FALSE)
}
if (nrow(fsetdiff(full_an_grid, unique(an_keys))) || nrow(fsetdiff(unique(an_keys), full_an_grid)) || anyDuplicated(an_keys)) {
  stop("Grouped AN domain is incomplete for Madrid.", call. = FALSE)
}

grouped_dem <- merge(full_dem_grid, grouped_dem, by = c("year", "agegroup"), all.x = TRUE, sort = FALSE)
grouped_an <- merge(full_an_grid, grouped_an, by = c("year", "agegroup", "range", "branch"), all.x = TRUE, sort = FALSE)

if (anyNA(grouped_dem$death) || anyNA(grouped_an$an)) {
  stop("Merged demographic or AN inputs contain NA values.", call. = FALSE)
}

#----- Derive Method A weights from Step 00's canonical single-age deaths

single_dem[, source_agegroup := fifelse(age <= 74L, "65-74", fifelse(age <= 84L, "75-84", "85+"))]
if (any(single_dem$agegroup != single_dem$source_agegroup)) {
  stop("Step 00 single-age rows do not match the required source age bands.", call. = FALSE)
}

weights_dt <- single_dem[, .(death = sum(death)), by = .(year, source_agegroup, age)]
weights_dt[, group_death := sum(death), by = .(year, source_agegroup)]
weights_dt[, weight := fifelse(group_death > 0, death / group_death, 0)]

expected_weight_grid <- rbindlist(lapply(age_groups, function(grp) {
  CJ(year = future_years, source_agegroup = grp, age = age_slices[[grp]])
}))
if (nrow(fsetdiff(expected_weight_grid, weights_dt[, .(year, source_agegroup, age)])) ||
    nrow(fsetdiff(weights_dt[, .(year, source_agegroup, age)], expected_weight_grid)) ||
    anyDuplicated(weights_dt, by = c("year", "source_agegroup", "age"))) {
  stop("Step 00 single-age demographic domain is incomplete or duplicated.", call. = FALSE)
}

gcm_values <- unique(grouped_an$gcm)
if (length(gcm_values) != 1L || is.na(gcm_values)) {
  stop("Grouped AN must contain exactly one GCM.", call. = FALSE)
}

allocation_rows <- list()
reconstruction_rows <- list()

#----- Loop years, then age groups, then branch/range to allocate AN to single ages

for (yr in sort(unique(grouped_dem$year))) {
  for (grp in age_groups) {
    weight_rows <- weights_dt[year == yr & source_agegroup == grp][order(age)]
    ages <- age_slices[[grp]]
    if (!identical(as.integer(weight_rows$age), as.integer(ages))) {
      stop(sprintf("Age-band mapping is invalid for %s/%s.", yr, grp), call. = FALSE)
    }
    w_vec <- weight_rows$weight
    group_death <- unique(weight_rows$group_death)
    if (length(group_death) != 1L) stop(sprintf("Ambiguous group deaths for %s/%s.", yr, grp), call. = FALSE)
    if (group_death > 0 && abs(sum(w_vec) - 1) > 1e-12) {
      stop(sprintf("Weight vector for %s/%s does not sum to 1.", yr, grp), call. = FALSE)
    }

    for (br in c("with_cc", "without_cc")) {
      for (rg in c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")) {
        an_row <- grouped_an[year == yr & agegroup == grp & branch == br & range == rg]
        if (nrow(an_row) != 1L) stop(sprintf("Missing grouped AN row for %s/%s/%s/%s.", yr, br, grp, rg), call. = FALSE)
        group_an <- an_row$an
        if (group_an != 0 && group_death == 0) {
          stop(sprintf("Nonzero AN with zero group allocation weight for %s/%s/%s/%s.", yr, br, grp, rg), call. = FALSE)
        }
        single_an <- group_an * w_vec
        if (abs(sum(single_an) - group_an) > 1e-9) {
          stop(sprintf("Grouped AN reconstruction error exceeds tolerance for %s/%s/%s/%s.", yr, br, grp, rg), call. = FALSE)
        }

        allocation_rows[[length(allocation_rows) + 1L]] <- data.table(
          geo_id = city_id,
          label = city_name,
          ssp = 3L,
          gcm = gcm_values,
          branch = br,
          year = yr,
          source_agegroup = grp,
          age = ages,
          range = rg,
          group_an = group_an,
          weight = as.numeric(w_vec),
          an = as.numeric(single_an)
        )

        reconstruction_rows[[length(reconstruction_rows) + 1L]] <- data.table(
          year = yr,
          branch = br,
          source_agegroup = grp,
          range = rg,
          group_an = group_an,
          reconstructed_an = sum(single_an),
          abs_diff = abs(sum(single_an) - group_an),
          max_weight = max(w_vec),
          min_weight = min(w_vec),
          zero_weight_count = sum(w_vec == 0),
          group_death = group_death
        )
      }
    }
  }
}

single_an <- rbindlist(allocation_rows, use.names = TRUE)
recon_dt <- rbindlist(reconstruction_rows, use.names = TRUE)

setorder(single_an, branch, year, source_agegroup, age, range)
setorder(weights_dt, year, source_agegroup, age)

#----- Invariant checks (project convention: a failing check stops the run)

age_band_ok <- single_an[, {
  expected_ages <- age_slices[[.BY$source_agegroup]]
  .(ok = all(age %in% expected_ages))
}, by = .(year, branch, source_agegroup, range)]
weight_sums <- weights_dt[, .(weight_sum = sum(weight), group_death = unique(group_death)), by = .(year, source_agegroup)]

checks <- data.table(
  check_name = c(
    "grouped_demography_domain_complete",
    "grouped_an_domain_complete",
    "demographic_weight_source_consistent",
    "weight_vectors_sum_to_one",
    "single_age_grouped_total_preserved",
    "no_nonzero_an_zero_weight",
    "age_band_mapping_exact"
  ),
  status = c(
    if (nrow(unique(grouped_dem, by = c("year", "agegroup"))) == nrow(full_dem_grid)) "PASS" else "FAIL",
    if (nrow(unique(grouped_an, by = c("year", "agegroup", "range", "branch"))) == nrow(full_an_grid)) "PASS" else "FAIL",
    if (max(abs(weights_dt[, sum(death), by = .(year, source_agegroup)]$V1 - weights_dt[, unique(group_death), by = .(year, source_agegroup)]$V1)) <= 1e-9) "PASS" else "FAIL",
    if (all(abs(weight_sums[group_death > 0]$weight_sum - 1) <= 1e-12) && all(weight_sums[group_death == 0]$weight_sum == 0)) "PASS" else "FAIL",
    if (max(recon_dt$abs_diff) <= 1e-9) "PASS" else "FAIL",
    if (!nrow(recon_dt[group_an != 0 & group_death == 0])) "PASS" else "FAIL",
    if (all(age_band_ok$ok)) "PASS" else "FAIL"
  ),
  value = c(
    nrow(grouped_dem),
    nrow(grouped_an),
    sprintf("max_abs_diff=%0.3e", max(abs(weights_dt[, sum(death), by = .(year, source_agegroup)]$V1 - weights_dt[, unique(group_death), by = .(year, source_agegroup)]$V1))),
    sprintf("min_positive_weight_sum=%0.12f; max_positive_weight_sum=%0.12f", min(weight_sums[group_death > 0]$weight_sum), max(weight_sums[group_death > 0]$weight_sum)),
    sprintf("max_abs_diff=%0.3e", max(recon_dt$abs_diff)),
    sprintf("nonzero_an_zero_group_weight_rows=%d", nrow(recon_dt[group_an != 0 & group_death == 0])),
    sprintf("bad_band_rows=%d", sum(!age_band_ok$ok))
  ),
  threshold = c(
    sprintf("%d rows", nrow(full_dem_grid)),
    sprintf("%d rows", nrow(full_an_grid)),
    "<= 1e-9",
    "sum to 1 within 1e-12",
    "<= 1e-9",
    "0 groups",
    "all TRUE"
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

#----- Persist checks before exposing primary outputs

fwrite(checks, checks_file)
if (nrow(failures)) {
  fwrite(failures, failures_file)
  stop(sprintf("02_single_age.R failed %d invariant(s); see %s", nrow(failures), failures_file), call. = FALSE)
} else {
  if (file.exists(failures_file)) file.remove(failures_file)
  invisible(file.create(failures_file))
}

fwrite(single_an, single_file)

min_year <- min(weights_dt$year)
weights_plot <- weights_dt[year == min_year][, .(age, weight, source_agegroup)]
plot_weights <- ggplot(weights_plot, aes(x = age, y = weight)) +
  geom_line(linewidth = 0.6, color = "#2c7fb8") +
  facet_wrap(~source_agegroup, scales = "free_x") +
  labs(
    title = "Madrid single-age demographic weights",
    subtitle = "Method A weights from Step 00 single-age all-cause deaths",
    x = "Age",
    y = "Weight"
  ) +
  theme_minimal(base_size = 11)

plot_recon <- ggplot(recon_dt, aes(x = group_an, y = reconstructed_an, color = branch)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(
    title = "Grouped-to-single-age AN reconstruction",
    subtitle = sprintf("Max abs diff = %.3e", max(recon_dt$abs_diff)),
    x = "Grouped AN",
    y = "Reconstructed AN"
  ) +
  theme_minimal(base_size = 11)

p <- plot_weights / plot_recon
ggsave(fig_file, p, width = 11, height = 9, dpi = 160)

message("Saved single-age AN to ", single_file)
message("Saved checks to ", checks_file)
message("Saved diagnostic figure to ", fig_file)
