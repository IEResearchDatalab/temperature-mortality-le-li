#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
output_dir <- Sys.getenv(
  "EXTREME_OUTPUT_DIR",
  "results/figures_simon_20oct/extreme_temperature"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

years <- 2015:2099
expected_year_count <- length(years)
gcm_names <- c(
  "ACCESS_CM2", "ACCESS_ESM1_5", "BCC_CSM2_MR", "CMCC_ESM2", "CanESM5",
  "EC_Earth3", "EC_Earth3_Veg_LR", "GFDL_ESM4", "IITM_ESM", "INM_CM4_8",
  "INM_CM5_0", "IPSL_CM6A_LR", "KACE_1_0_G", "MIROC6", "MPI_ESM1_2_HR",
  "MPI_ESM1_2_LR", "MRI_ESM2_0", "NorESM2_LM", "NorESM2_MM"
)
ssp_labels <- c(`1` = "SSP1–2.6", `2` = "SSP2–4.5", `3` = "SSP3–7.0")
ssp_order <- c("1", "2", "3")
scenario_colours <- c(
  "SSP1–2.6" = "#24558A",
  "SSP2–4.5" = "#F06B18",
  "SSP3–7.0" = "#AD1A2B"
)
region_order <- c("Northern", "Western", "Eastern", "Southern")
footnote_for <- function(version) {
  summary <- if (version == "A") {
    "Mean of 19 GCMs; mean, maximum and minimum across cities."
  } else {
    "Mean, maximum and minimum across GCMs."
  }
  text <- paste(
    "The extreme-heat threshold for each year is the historical 97.5th",
    "percentile (ERA5-Land, 2000–2014) plus the change in the calibrated GCM",
    "97.5th percentile over the 15 preceding years; it never falls below the",
    "historical value. Without climate change the threshold stays at the",
    "historical value. For SSP1, Wittgenstein Centre survival ratios stored as",
    "exactly 1 were set to 0.9995 (this does not affect temperatures).",
    summary
  )
  paste(strwrap(text, width = 145L), collapse = "\n")
}

metadata_file <- file.path("data", "city_results.csv")
if (!file.exists(metadata_file)) stop("Missing city metadata: ", metadata_file, call. = FALSE)
city_metadata <- unique(
  fread(metadata_file, select = c("URAU_CODE", "region"))[, .(city = URAU_CODE, region)]
)
if (anyNA(city_metadata) || anyDuplicated(city_metadata$city)) {
  stop("City metadata must map each city to exactly one region.", call. = FALSE)
}
expected_cities <- sort(city_metadata$city)
if (length(expected_cities) != 854L) {
  stop("Expected 854 cities in data/city_results.csv; found ", length(expected_cities), call. = FALSE)
}
if (!setequal(unique(city_metadata$region), region_order)) {
  stop("City metadata does not contain the expected four regions.", call. = FALSE)
}

read_scenario <- function(ssp) {
  ssp_dir <- file.path(batch_root, paste0("ssp", ssp))
  files <- list.files(
    ssp_dir, pattern = "^01_moving_heat_thresholds\\.csv$",
    recursive = TRUE, full.names = TRUE
  )
  if (!length(files)) stop("No moving-threshold checks found under ", ssp_dir, call. = FALSE)

  paths <- data.table(file = files)
  paths[, city := basename(dirname(dirname(dirname(file))))]
  paths[, gcm := basename(dirname(dirname(file)))]
  missing_pairs <- CJ(city = expected_cities, gcm = gcm_names)[
    !paths[, .(city, gcm)], on = .(city, gcm)
  ]
  unexpected_pairs <- paths[!city %chin% expected_cities | !gcm %chin% gcm_names]
  pair_counts <- paths[, .N, by = .(city, gcm)]
  duplicate_pairs <- pair_counts[N != 1L]
  if (nrow(missing_pairs) || nrow(unexpected_pairs) || nrow(duplicate_pairs)) {
    print(list(
      missing_city_gcm_pairs = missing_pairs,
      unexpected_files = unexpected_pairs,
      duplicate_city_gcm_pairs = duplicate_pairs
    ))
    stop("Incomplete city/GCM threshold-table domain for SSP", ssp, ".", call. = FALSE)
  }
  cat(sprintf(
    "SSP%s: %d cities × %d GCMs; missing city/GCM pairs: 0\n",
    ssp, uniqueN(paths$city), uniqueN(paths$gcm)
  ))

  tables <- lapply(seq_len(nrow(paths)), function(i) {
    x <- fread(
      paths$file[i],
      select = c(
        "year", "static_p97_5", "moving_p97_5",
        "moving_p97_5_anchored_raw", "clamped"
      )
    )
    x[, `:=`(city = paths$city[i], gcm = paths$gcm[i], ssp = as.character(ssp))]
    x
  })
  d <- rbindlist(tables)
  expected_years <- data.table(year = years)
  year_domain <- d[, .(years = uniqueN(year), missing = paste(setdiff(years, year), collapse = ",")), by = .(city, gcm)]
  bad_years <- year_domain[years != expected_year_count | missing != ""]
  if (nrow(bad_years)) {
    print(bad_years)
    stop("One or more city/GCM tables do not contain exactly 2015–2099.", call. = FALSE)
  }
  if (nrow(d) != length(expected_cities) * length(gcm_names) * length(years) ||
      anyDuplicated(d, by = c("city", "gcm", "year")) ||
      !setequal(unique(d$year), expected_years$year)) {
    stop("Threshold tables do not cover every city, GCM, and year exactly once.", call. = FALSE)
  }
  cat(sprintf(
    "SSP%s: all city-GCM pairs have %d years; missing years: 0\n",
    ssp, length(years)
  ))

  d[, difference_from_static := moving_p97_5 - static_p97_5]
  d[, expected_clamped := moving_p97_5_anchored_raw < static_p97_5]
  if (anyNA(d[, .(static_p97_5, moving_p97_5, moving_p97_5_anchored_raw, clamped)]) ||
      any(d$difference_from_static < 0) ||
      any(d$clamped != d$expected_clamped)) {
    stop("Thresholds contain missing, below-static, or inconsistent clamp values for SSP", ssp, ".", call. = FALSE)
  }
  baseline_diff <- d[year == 2015L, max(abs(difference_from_static))]
  if (length(baseline_diff) != 1L || baseline_diff > 1e-8) {
    stop("2015 moving thresholds do not match the static threshold for SSP", ssp, ".", call. = FALSE)
  }
  cat(sprintf(
    "SSP%s: max |2015 moving − static| = %.12g °C; below-static values: 0; ",
    ssp, baseline_diff
  ))
  cat(sprintf("clamped city-GCM-years = %d\n", sum(d$clamped)))
  d
}

thresholds <- rbindlist(lapply(ssp_order, read_scenario))
thresholds <- merge(thresholds, city_metadata, by = "city", all.x = TRUE, sort = FALSE)
thresholds[, scenario := factor(ssp_labels[ssp], levels = ssp_labels[ssp_order])]

region_city_set <- c(
  list(Europe = expected_cities),
  setNames(lapply(region_order, function(x) city_metadata[region == x, city]), region_order)
)

summarise_geography <- function(cities, geography) {
  d <- thresholds[city %chin% cities]
  city_gcm_means <- d[, .(
    change = mean(difference_from_static),
    static_threshold = unique(static_p97_5),
    moving_threshold = mean(moving_p97_5)
  ), by = .(ssp, scenario, city, gcm, year)]
  if (d[, uniqueN(static_p97_5), by = .(city, gcm)][, any(V1 != 1L)]) {
    stop("Static thresholds are not unique within a city-GCM-year.", call. = FALSE)
  }

  version_a_city <- city_gcm_means[, .(
    change = mean(change),
    static_threshold = unique(static_threshold),
    moving_threshold = mean(moving_threshold)
  ), by = .(ssp, scenario, city, year)]
  version_a <- version_a_city[, .(
    mean = mean(change),
    minimum = min(change),
    maximum = max(change),
    n_units = uniqueN(city)
  ), by = .(ssp, scenario, year)]
  version_a[, `:=`(
    version = "A",
    geography = geography,
    clamped_city_gcm_years = d[
      , .(clamped_city_gcm_years = sum(clamped)), by = .(ssp, year)
    ][version_a, on = .(ssp, year), clamped_city_gcm_years]
  )]
  version_a[, clamped_city_gcm_share_percent :=
              100 * clamped_city_gcm_years / (length(cities) * length(gcm_names))]

  independent_mean <- aggregate(
    version_a_city$change,
    by = list(ssp = version_a_city$ssp, year = version_a_city$year),
    FUN = mean
  )
  data_table_mean <- version_a[scenario == ssp_labels[ssp], .(ssp, year, mean)]
  setnames(independent_mean, "x", "independent_mean")
  independent_check <- merge(data_table_mean, independent_mean, by = c("ssp", "year"))
  max_independent_diff <- max(abs(independent_check$mean - independent_check$independent_mean))
  if (max_independent_diff > 1e-12) {
    stop("Version A independent mean calculation differs by more than 1e-12 °C.", call. = FALSE)
  }
  cat(sprintf(
    "%s: Version A independent mean max difference = %.12g °C\n",
    geography, max_independent_diff
  ))

  version_b_model <- city_gcm_means[, .(
    change = mean(change)
  ), by = .(ssp, scenario, gcm, year)]
  version_b <- version_b_model[, .(
    mean = mean(change),
    minimum = min(change),
    maximum = max(change),
    n_units = uniqueN(gcm)
  ), by = .(ssp, scenario, year)]
  version_b[, `:=`(
    version = "B",
    geography = geography,
    clamped_city_gcm_years = d[
      , .(clamped_city_gcm_years = sum(clamped)), by = .(ssp, year)
    ][version_b, on = .(ssp, year), clamped_city_gcm_years]
  )]
  version_b[, clamped_city_gcm_share_percent :=
              100 * clamped_city_gcm_years / (length(cities) * length(gcm_names))]
  version_a_city[, `:=`(geography = geography, ssp = as.character(ssp))]

  list(
    version_a = version_a,
    version_b = version_b,
    city_data = version_a_city,
    model_data = version_b_model,
    max_independent_diff = max_independent_diff,
    clamps = d[, .(clamped_city_gcm_years = sum(clamped)), by = .(ssp, year)]
  )
}

make_change_plot <- function(summary, version, geography, y_range) {
  geography_title <- if (geography == "Europe") "Europe" else paste(geography, "Europe")
  title <- if (version == "A") {
    paste0(
      geography_title,
      ": how far the extreme-heat threshold has moved above its 2000–2014 level,\nby year"
    )
  } else {
    paste0(
      geography_title,
      ": how far the extreme-heat threshold has moved above its 2000–2014 level,\nby climate model"
    )
  }
  scenario_summary <- copy(summary)
  scenario_summary[, scenario := factor(scenario, levels = ssp_labels[ssp_order])]

  p <- ggplot(scenario_summary, aes(x = year)) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.45) +
    geom_ribbon(
      aes(ymin = minimum, ymax = maximum, fill = scenario, group = scenario),
      alpha = 0.18, colour = NA
    ) +
    geom_line(aes(y = maximum, colour = scenario, group = scenario), linewidth = 0.45) +
    geom_line(aes(y = minimum, colour = scenario, group = scenario), linewidth = 0.45) +
    geom_line(aes(y = mean, colour = scenario, group = scenario), linewidth = 1.05) +
    facet_grid(rows = vars(scenario), switch = "y") +
    scale_colour_manual(values = scenario_colours, guide = "none") +
    scale_fill_manual(values = scenario_colours, guide = "none") +
    scale_x_continuous(
      breaks = seq(2015, 2095, by = 10),
      limits = c(2015, 2099),
      expand = expansion(mult = c(0, 0.01)),
      labels = scales::label_number(big.mark = "", style_negative = "minus")
    ) +
    scale_y_continuous(
      labels = scales::label_number(style_negative = "minus"),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(ylim = y_range, clip = "off") +
    labs(
      title = title,
      subtitle = paste(
        if (version == "A") "Version A: across cities." else "Version B: across GCMs.",
        "Thick line: mean; ribbon and thin lines: minimum to maximum."
      ),
      x = "Year",
      y = "Moving 97.5th percentile minus static 97.5th percentile (°C)",
      caption = footnote_for(version)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.spacing.y = grid::unit(0.15, "lines"),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.title = element_text(face = "bold", size = 9.5, lineheight = 1),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 7, hjust = 0, margin = margin(t = 8)),
      plot.caption.position = "plot",
      plot.margin = margin(6, 10, 6, 6)
    )
  p
}

make_absolute_plot <- function(city_data) {
  absolute <- city_data[, .(
    moving_mean = mean(moving_threshold),
    static_mean = mean(static_threshold)
  ), by = year]
  ggplot(absolute, aes(x = year)) +
    geom_line(aes(y = moving_mean), colour = "#24558A", linewidth = 1.05) +
    geom_hline(
      aes(yintercept = static_mean),
      data = absolute[1L], colour = "grey35", linetype = "dashed", linewidth = 0.6
    ) +
    geom_text(
      data = absolute[year == max(year)],
      aes(x = year, y = moving_mean, label = "Moving threshold"),
      inherit.aes = FALSE, hjust = -0.05, vjust = 1.2, size = 2.8,
      colour = "#24558A"
    ) +
    geom_text(
      data = absolute[year == max(year)],
      aes(x = year, y = static_mean, label = "Static threshold (2000–2014)"),
      inherit.aes = FALSE, hjust = -0.05, vjust = -0.2, size = 2.8,
      colour = "grey35"
    ) +
    scale_x_continuous(
      breaks = seq(2015, 2095, by = 10),
      limits = c(2015, 2118),
      expand = expansion(mult = c(0, 0)),
      labels = scales::label_number(big.mark = "", style_negative = "minus")
    ) +
    scale_y_continuous(
      labels = scales::label_number(style_negative = "minus")
    ) +
    labs(
      x = "Year",
      y = "Mean-city extreme-heat threshold (°C)"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(margin = margin(r = 8))
    )
}

write_figure <- function(plot, stem, data) {
  fwrite(data, file.path(output_dir, paste0(stem, ".csv")))
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width = 10, height = 8, dpi = 300)
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = 10, height = 8, device = grDevices::cairo_pdf
  )
}

for (geography in names(region_city_set)) {
  cities <- region_city_set[[geography]]
  result <- summarise_geography(cities, geography)
  geography_title <- if (geography == "Europe") "Europe" else paste(geography, "Europe")
  label <- gsub(" ", "_", geography_title, fixed = TRUE)
  shared_limits <- range(c(
    result$version_a$minimum, result$version_a$maximum,
    result$version_b$minimum, result$version_b$maximum, 0
  ), finite = TRUE)
  shared_span <- diff(shared_limits)
  if (!is.finite(shared_span) || shared_span == 0) shared_span <- 1
  shared_y_range <- c(
    shared_limits[1] - shared_span * 0.08,
    shared_limits[2] + shared_span * 0.08
  )
  cat(sprintf(
    "%s shared Version A/B y-limits: %.6f to %.6f °C\n",
    geography, shared_y_range[1], shared_y_range[2]
  ))

  for (version in c("A", "B")) {
    summary <- if (version == "A") result$version_a else result$version_b
    plot <- make_change_plot(summary, version, geography, shared_y_range)
    if (version == "A" && geography == "Europe") {
      absolute <- result$city_data[, .(
        mean_city_absolute_threshold = mean(moving_threshold),
        mean_city_static_threshold = mean(static_threshold)
      ), by = year]
      data <- merge(copy(summary), absolute, by = "year", all.x = TRUE)
      plot <- ((plot + labs(caption = NULL)) / make_absolute_plot(result$city_data)) +
        plot_layout(heights = c(3.2, 1.1)) +
        plot_annotation(
          caption = footnote_for("A"),
          theme = theme(
            plot.caption = element_text(size = 7, hjust = 0, margin = margin(t = 8)),
            plot.caption.position = "plot"
          )
        )
    } else {
      data <- copy(summary)
      data[, `:=`(
        mean_city_absolute_threshold = NA_real_,
        mean_city_static_threshold = NA_real_
      )]
    }
    stem <- sprintf("extreme_temperature_%s_%s", version, label)
    write_figure(plot, stem, data)
    clamp_totals <- result$clamps[, .(total = sum(clamped_city_gcm_years)), by = ssp][order(ssp)]
    cat(sprintf(
      "%s Version %s clamp counts by SSP (2015–2099): %s\n",
      geography, version,
      paste(sprintf("%s=%d", ssp_labels[clamp_totals$ssp], clamp_totals$total), collapse = ", ")
    ))
  }

  if (geography == "Europe") {
    absolute_plot <- make_absolute_plot(result$city_data) +
      labs(
        title = "Europe: absolute extreme-heat threshold for the mean city",
        caption = footnote_for("A")
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.caption = element_text(size = 7, hjust = 0, margin = margin(t = 8)),
        plot.caption.position = "plot"
      )
    absolute_data <- result$city_data[, .(
      mean_city_absolute_threshold = mean(moving_threshold),
      mean_city_static_threshold = mean(static_threshold)
    ), by = year]
    write_figure(
      absolute_plot,
      "extreme_temperature_A_Europe_absolute",
      absolute_data
    )
  }
}
