#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("R_pipeline/functions/pclm_utils.R")

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
age_groups <- c("65-74", "75-84", "85+")
age_slices <- list(
  "65-74" = 65:74,
  "75-84" = 75:84,
  "85+" = 85:100
)
age_starts <- c("65-74" = 65L, "75-84" = 75L, "85+" = 85L)
nlast <- 16L

grouped_dem <- fread(file.path(out_dir, "00_demography_grouped.csv"))
grouped_an <- fread(file.path(out_dir, "01_attribution_grouped.csv"))

grouped_dem <- grouped_dem[geo_id == city_id & ssp == 3]
grouped_an <- grouped_an[geo_id == city_id & ssp == 3]

if (!nrow(grouped_dem)) stop("Grouped demography for Madrid is missing; run 00_demography.R first.", call. = FALSE)
if (!nrow(grouped_an)) stop("Grouped AN for Madrid is missing; run 01_attribution.R first.", call. = FALSE)

grouped_dem <- grouped_dem[agegroup %in% age_groups]
grouped_an <- grouped_an[agegroup %in% age_groups]

if (any(!is.finite(grouped_dem$death)) || any(grouped_dem$death < 0)) {
  stop("Invalid grouped demographic deaths detected.", call. = FALSE)
}
if (any(!is.finite(grouped_an$an))) {
  stop("Grouped AN contains non-finite values.", call. = FALSE)
}

full_dem_grid <- CJ(year = sort(unique(grouped_dem$year)), agegroup = age_groups)
full_an_grid <- CJ(year = sort(unique(grouped_an$year)), agegroup = age_groups, range = c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat"), branch = c("with_cc", "without_cc"))

if (nrow(unique(grouped_dem, by = c("year", "agegroup"))) != nrow(full_dem_grid)) {
  stop("Grouped demographic domain is incomplete for Madrid.", call. = FALSE)
}
if (nrow(unique(grouped_an, by = c("year", "agegroup", "range", "branch"))) != nrow(full_an_grid)) {
  stop("Grouped AN domain is incomplete for Madrid.", call. = FALSE)
}

grouped_dem <- merge(full_dem_grid, grouped_dem, by = c("year", "agegroup"), all.x = TRUE, sort = FALSE)
grouped_an <- merge(full_an_grid, grouped_an, by = c("year", "agegroup", "range", "branch"), all.x = TRUE, sort = FALSE)

if (anyNA(grouped_dem$death) || anyNA(grouped_an$an)) {
  stop("Merged demographic or AN inputs contain NA values.", call. = FALSE)
}

build_weights <- function(death_vec, year, agegroup) {
  if (any(!is.finite(death_vec)) || any(death_vec < 0)) {
    stop(sprintf("Invalid demographic death input for %s/%s.", year, agegroup), call. = FALSE)
  }
  if (sum(death_vec) == 0) {
    return(rep(0, pclm_expected_length(c(65, 75, 85), nlast)))
  }
  fit <- pclm_disaggregate_nonnegative(x = c(65, 75, 85), y = death_vec, nlast = nlast)
  if (any(!is.finite(fit)) || any(fit < 0)) {
    stop(sprintf("PCLM output is invalid for %s/%s.", year, agegroup), call. = FALSE)
  }
  if (abs(sum(fit) - sum(death_vec)) > 1e-9) {
    stop(sprintf("Grouped reconstruction error exceeds tolerance for %s/%s.", year, agegroup), call. = FALSE)
  }
  slices <- list(
    "65-74" = fit[1:10],
    "75-84" = fit[11:20],
    "85+" = fit[21:36]
  )
  target_bins <- as.numeric(death_vec)
  for (i in seq_along(slices)) {
    s <- slices[[i]]
    if (sum(s) == 0) {
      slices[[i]] <- rep(0, length(s))
    } else {
      slices[[i]] <- s * (target_bins[i] / sum(s))
    }
  }
  if (max(abs(c(sum(slices[[1]]), sum(slices[[2]]), sum(slices[[3]])) - target_bins)) > 1e-9) {
    stop(sprintf("PCLM grouped-bin reconstruction failed for %s/%s.", year, agegroup), call. = FALSE)
  }
  weight_list <- lapply(names(slices), function(grp) {
    vec <- slices[[grp]]
    if (sum(vec) == 0) {
      rep(0, length(vec))
    } else {
      vec / sum(vec)
    }
  })
  names(weight_list) <- names(slices)
  list(fit = fit, weights = weight_list)
}

weights_rows <- list()
allocation_rows <- list()
reconstruction_rows <- list()

for (yr in sort(unique(grouped_dem$year))) {
  dem_year <- grouped_dem[year == yr][match(age_groups, agegroup)]
  if (nrow(dem_year) != length(age_groups)) stop(sprintf("Missing demographic rows for year %s.", yr), call. = FALSE)

  weight_obj <- build_weights(death_vec = dem_year$death, year = yr, agegroup = "65-74/75-84/85+")
  fit <- weight_obj$fit
  weight_slices <- list(
    "65-74" = fit[1:10],
    "75-84" = fit[11:20],
    "85+" = fit[21:36]
  )

  for (grp in age_groups) {
    w_vec <- weight_slices[[grp]]
    ages <- age_slices[[grp]]
    if (sum(w_vec) > 0) w_vec <- w_vec / sum(w_vec)
    weights_rows[[length(weights_rows) + 1L]] <- data.table(
      year = yr,
      source_agegroup = grp,
      age = ages,
      weight = as.numeric(w_vec)
    )
    if (abs(sum(w_vec) - 1) > 1e-12) {
      stop(sprintf("Weight vector for %s/%s does not sum to 1.", yr, grp), call. = FALSE)
    }

    for (br in c("with_cc", "without_cc")) {
      for (rg in c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")) {
        an_row <- grouped_an[year == yr & agegroup == grp & branch == br & range == rg]
        if (nrow(an_row) != 1L) stop(sprintf("Missing grouped AN row for %s/%s/%s/%s.", yr, br, grp, rg), call. = FALSE)
        group_an <- an_row$an
        if (group_an != 0 && any(w_vec == 0)) {
          stop(sprintf("Nonzero AN with zero allocation weight for %s/%s/%s/%s.", yr, br, grp, rg), call. = FALSE)
        }
        single_an <- group_an * w_vec
        if (abs(sum(single_an) - group_an) > 1e-9) {
          stop(sprintf("Grouped AN reconstruction error exceeds tolerance for %s/%s/%s/%s.", yr, br, grp, rg), call. = FALSE)
        }

        allocation_rows[[length(allocation_rows) + 1L]] <- data.table(
          geo_id = city_id,
          label = city_name,
          ssp = 3L,
          gcm = grouped_an$gcm[1],
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
          zero_weight_count = sum(w_vec == 0)
        )
      }
    }
  }
}

weights_dt <- rbindlist(weights_rows, use.names = TRUE)
single_an <- rbindlist(allocation_rows, use.names = TRUE)
recon_dt <- rbindlist(reconstruction_rows, use.names = TRUE)

setorder(single_an, branch, year, source_agegroup, age, range)
setorder(weights_dt, year, source_agegroup, age)

age_band_ok <- single_an[, all(age %in% age_slices[[source_agegroup]]), by = .(year, branch, source_agegroup, range)]

checks <- data.table(
  check_name = c(
    "grouped_demography_domain_complete",
    "grouped_an_domain_complete",
    "pclm_group_reconstruction",
    "weight_vectors_sum_to_one",
    "single_age_grouped_total_preserved",
    "no_nonzero_an_zero_weight",
    "age_band_mapping_exact"
  ),
  status = c(
    if (nrow(unique(grouped_dem, by = c("year", "agegroup"))) == nrow(full_dem_grid)) "PASS" else "FAIL",
    if (nrow(unique(grouped_an, by = c("year", "agegroup", "range", "branch"))) == nrow(full_an_grid)) "PASS" else "FAIL",
    if (max(recon_dt$abs_diff) <= 1e-9) "PASS" else "FAIL",
    if (all(abs(weights_dt[, sum(weight), by = .(year, source_agegroup)]$V1 - 1) <= 1e-12)) "PASS" else "FAIL",
    if (max(recon_dt$abs_diff) <= 1e-9) "PASS" else "FAIL",
    if (all(recon_dt[group_an != 0]$zero_weight_count == 0L)) "PASS" else "FAIL",
    if (all(age_band_ok$V1)) "PASS" else "FAIL"
  ),
  value = c(
    nrow(grouped_dem),
    nrow(grouped_an),
    sprintf("max_abs_diff=%0.3e", max(recon_dt$abs_diff)),
    sprintf("min_weight_sum=%0.12f; max_weight_sum=%0.12f", min(weights_dt[, sum(weight), by = .(year, source_agegroup)]$V1), max(weights_dt[, sum(weight), by = .(year, source_agegroup)]$V1)),
    sprintf("max_abs_diff=%0.3e", max(recon_dt$abs_diff)),
    sprintf("nonzero_zero_weight_rows=%d", nrow(recon_dt[group_an != 0 & zero_weight_count > 0])),
    sprintf("bad_band_rows=%d", sum(!age_band_ok$V1))
  ),
  threshold = c(
    sprintf("%d rows", nrow(full_dem_grid)),
    sprintf("%d rows", nrow(full_an_grid)),
    "<= 1e-9",
    "sum to 1 within 1e-12",
    "<= 1e-9",
    "0 rows",
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

fwrite(single_an, single_file)
fwrite(checks, checks_file)
if (nrow(failures)) {
  fwrite(failures, failures_file)
} else {
  if (file.exists(failures_file)) file.remove(failures_file)
  invisible(file.create(failures_file))
}

min_year <- min(weights_dt$year)
weights_plot <- weights_dt[year == min_year][, .(age, weight, source_agegroup)]
plot_weights <- ggplot(weights_plot, aes(x = age, y = weight)) +
  geom_line(linewidth = 0.6, color = "#2c7fb8") +
  facet_wrap(~source_agegroup, scales = "free_x") +
  labs(
    title = "Madrid single-age demographic weights",
    subtitle = "Corrected Method A weights derived from grouped deaths via PCLM",
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

if (nrow(failures)) {
  stop(sprintf("02_single_age.R failed %d invariant(s); see %s", nrow(failures), failures_file), call. = FALSE)
}

message("Saved single-age AN to ", single_file)
message("Saved checks to ", checks_file)
message("Saved diagnostic figure to ", fig_file)
