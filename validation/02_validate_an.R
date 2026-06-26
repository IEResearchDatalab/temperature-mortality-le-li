################################################################################
#
# Validation: Step 02 — Climate-only attributable number computation
#   Diagnose ANs over time, by temperature range, by age group.
#
# Plots:
#   02a_ans_timeseries.png  — Annual ANs by temp range over all projection years
#   02b_ans_by_age.png      — ANs by age group for a representative year
#   02c_mmt_check.png       — MMT verification: RR curves with MMT marked
#   02d_temp_dist.png       — Temperature distribution with thresholds marked
#
################################################################################

library(ggplot2); library(scales); library(data.table); library(dlnm)

validation_dir <- "../results/validation_plots"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Validation Step 02: Climate-only AN computation\n")
cat("========================================\n\n")

# If 02 hasn't been run yet, check for checkpoint files
ans_file <- file.path(path_out, "ans_annual_all_cities.csv")
if (!file.exists(ans_file)) {
  checkpoint_dir <- file.path(path_out, "an_checkpoint")
  if (dir.exists(checkpoint_dir)) {
    ans_files <- list.files(checkpoint_dir, pattern = "\\.csv$", full.names = TRUE)
    if (length(ans_files) > 0) {
      ans <- rbindlist(lapply(ans_files, fread), fill = TRUE)
      cat(sprintf("  Loaded from checkpoints: %d rows\n", nrow(ans)))
    } else {
      cat("  WARNING: No AN data found. Run 02_ComputeAN.R first.\n")
      ans <- data.table()
    }
  } else {
    cat("  WARNING: No AN data found. Run 02_ComputeAN.R first.\n")
    ans <- data.table()
  }
} else {
  ans <- fread(ans_file)
}

if (nrow(ans) > 0) {

  cat(sprintf("  Total AN rows: %s\n", format(nrow(ans), big.mark = ",")))
  cat(sprintf("  Cities: %d\n", uniqueN(ans$city_code)))
  cat(sprintf("  GCMs: %d\n", uniqueN(ans$gcm)))
  cat(sprintf("  SSPs: %s\n", paste(unique(ans$ssp), collapse = ", ")))
  cat(sprintf("  Temp ranges: %s\n", paste(unique(ans$temp_range), collapse = ", ")))
  cat(sprintf("  Years: %d-%d\n", min(ans$year), max(ans$year)))

  #---------------------------
  # 02a: AN time series by temperature range
  #---------------------------
  # Sum across age groups and GCMs for a clear view of the climate signal

  ans_ts <- ans[, .(total_AN = sum(an_est)), by = .(year, ssp, temp_range)]
  # Average across GCMs for display
  ans_ts <- ans[, .(total_AN = mean(an_est)), by = .(year, ssp, temp_range)]
  setorder(ans_ts, ssp, temp_range, year)

  p <- ggplot(ans_ts, aes(x = year, y = total_AN, color = temp_range)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ssp, scales = "free_y",
      labeller = labeller(ssp = ssp_labels)) +
    labs(title = "Madrid: Annual temperature-attributable deaths by range",
         subtitle = "Climate-only ANs under fixed baseline deaths (all GCMs averaged)",
         caption = "Attribution using Masselot ERFs with fixed (2022) baseline deaths.
Shaded ranges: extreme cold (<p2.5), moderate cold (p2.5–MMT),
moderate heat (MMT–p97.5), extreme heat (>p97.5).",
         x = "Year", y = "Attributable deaths (annual)") +
    scale_y_continuous(labels = comma) +
    theme_minimal() + theme(legend.position = "bottom")

  ggsave(file.path(validation_dir, "02a_ans_timeseries.png"), p, width = 12, height = 7)
  cat(sprintf("  02a: AN time series saved\n"))

  #---------------------------
  # 02b: AN by age group (representative year and SSP)
  #---------------------------
  # Pick a mid-century year to show age profile

  ans_age <- ans[year == 2050 & ssp == "2",
    .(total_AN = mean(an_est)), by = .(age_group, temp_range)]

  p <- ggplot(ans_age, aes(x = age_group, y = total_AN, fill = temp_range)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Madrid 2050 SSP3-7.0: ANs by age group and temperature range",
         subtitle = "Averaged across GCMs. Stacked bars show contribution of each range.",
         caption = "Age groups defined by Masselot (2023): 20-44, 45-64, 65-74, 75-84, 85+.",
         x = "Age group", y = "Attributable deaths") +
    scale_y_continuous(labels = comma) +
    theme_minimal() + theme(legend.position = "bottom")

  ggsave(file.path(validation_dir, "02b_ans_by_age.png"), p, width = 9, height = 6)
  cat(sprintf("  02b: AN by age saved\n"))

  #---------------------------
  # 02c: MMT check
  #---------------------------
  # For Madrid, recreate the RR curve and mark MMT

  coefs <- readRDS(file.path(path_out, "coefs_all.rds"))
  madrid_coefs <- coefs[URAU_CODE == "ES001C"]

  # Load city metadata for MMT reference
  city_res <- readRDS(file.path(path_out, "city_res.rds"))
  madrid_mmt_ref <- unique(city_res[URAU_CODE == "ES001C", .(agegroup, mmt)])

  temp_seq <- seq(-10, 40, by = 0.5)
  argvar <- list(fun = varfun, degree = vardegree,
    knots = quantile(temp_seq, varper / 100),
    Bound = range(temp_seq))
  basis <- do.call(onebasis, c(list(x = temp_seq), argvar))

  rr_data <- rbindlist(lapply(unique(madrid_coefs$agegroup), function(ag) {
    coef_row <- madrid_coefs[agegroup == ag]
    coef_vec <- as.numeric(coef_row[, .(b1, b2, b3, b4, b5)])
    log_rr <- as.numeric(basis %*% coef_vec)
    search_idx <- which(temp_seq >= -5 & temp_seq <= 35)
    mmt <- temp_seq[search_idx][which.min(log_rr[search_idx])]
    cenvec <- do.call(onebasis, c(list(x = mmt), argvar))
    rr <- as.numeric(exp(log_rr - drop(cenvec %*% coef_vec)))
    data.table(temp = temp_seq, agegroup = ag, RR = rr,
      mmt_computed = mmt,
      mmt_reference = madrid_mmt_ref[agegroup == ag, mmt])
  }))

  p <- ggplot(rr_data, aes(x = temp, y = RR, color = agegroup)) +
    geom_line(linewidth = 1) +
    geom_vline(aes(xintercept = mmt_reference, color = agegroup),
      linetype = "dashed", alpha = 0.4) +
    scale_y_log10() +
    labs(title = "Madrid: RR curves with MMT verification",
         subtitle = "Solid = RR curve. Dashed = MMT from city_results.csv reference.",
         caption = "Vertical dashed lines = MMT per age group from original Masselot estimates.
Alignment of computed and reference MMT confirms correct basis and search range.",
         x = "Temperature (°C)", y = "Relative Risk") +
    geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.4) +
    theme_minimal() + theme(legend.position = "bottom")

  ggsave(file.path(validation_dir, "02c_mmt_check.png"), p, width = 10, height = 6)
  cat(sprintf("  02c: MMT check saved\n"))

  #---------------------------
  # 02d: Temperature distribution with thresholds
  #---------------------------

  # Get temperature data for Madrid from the parquet
  library(arrow)
  temp_data <- open_dataset(file.path(path_data, "tmeanproj.gz.parquet")) |>
    filter(URAU_CODE == "ES001C", ssp == "hist") |>
    collect() |>
    as.data.table()

  gcm_cols <- setdiff(names(temp_data), c("URAU_CODE", "date", "ssp"))
  hist_temps <- unlist(temp_data[, ..gcm_cols], use.names = FALSE)
  hist_temps <- hist_temps[!is.na(hist_temps)]

  p025 <- quantile(hist_temps, 0.025, na.rm = TRUE)
  p975 <- quantile(hist_temps, 0.975, na.rm = TRUE)

  temp_df <- data.table(temp = hist_temps)
  p <- ggplot(temp_df, aes(x = temp)) +
    geom_histogram(bins = 60, fill = "grey70", color = "grey30", alpha = 0.7) +
    geom_vline(xintercept = c(p025, p975), color = c("#2c3e50", "#e74c3c"),
      linewidth = 1, linetype = "dashed") +
    annotate("text", x = p025 - 1.5, y = Inf, label = "p2.5", vjust = 2,
      color = "#2c3e50", size = 4) +
    annotate("text", x = p975 + 1.5, y = Inf, label = "p97.5", vjust = 2,
      color = "#e74c3c", size = 4) +
    labs(title = "Madrid: Historical temperature distribution (all GCMs, 1990-2019)",
         subtitle = paste0("Dashed lines = extreme cold (p2.5=", round(p025, 1),
           "°C) and extreme heat (p97.5=", round(p975, 1), "°C) thresholds"),
         caption = "Temperature distribution from CMIP6 multi-GCM ensemble.
Extreme cold: below p2.5, extreme heat: above p97.5 of historical temps.",
         x = "Temperature (°C)", y = "Days") +
    theme_minimal()

  ggsave(file.path(validation_dir, "02d_temp_dist.png"), p, width = 10, height = 6)
  cat(sprintf("  02d: Temperature distribution saved\n"))

} else {
  cat("  No AN data available. Skipping plots.\n")
}

cat("\n--- Step 02 validation complete ---\n")