#!/usr/bin/env Rscript

# Figure 1: annual absolute risk by age and temperature range.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

batch_root <- Sys.getenv("BATCH_ROOT", "results/europe/heatsplit_20261013")
figure_root <- Sys.getenv("FIGURE_ROOT", "results/figures_simon_20oct")
geo_level <- tolower(Sys.getenv("GEO_LEVEL", "all"))
geo_id_filter <- Sys.getenv("GEO_ID", "")
ssp_set <- as.integer(strsplit(Sys.getenv("SSP_SET", "1,2,3"), ",", fixed = TRUE)[[1]])
if (!length(ssp_set) || anyNA(ssp_set) || any(!ssp_set %in% 1:3))
  stop("SSP_SET must be a comma-separated subset of 1,2,3.", call. = FALSE)
if (!geo_level %chin% c("all", "europe", "region"))
  stop("GEO_LEVEL must be all, europe, or region.", call. = FALSE)
dir.create(file.path(figure_root, "data"), recursive = TRUE, showWarnings = FALSE)

ranges <- c("ExtrCold", "ModCold", "ModHeat", "NormExtrHeat", "ExtrHeat")
range_labels <- c(ExtrCold = "Extreme cold", ModCold = "Moderate cold",
                  ModHeat = "Moderate heat", NormExtrHeat = "Normal extreme heat",
                  ExtrHeat = "Extreme heat")
range_colours <- c(ExtrCold = "#24558A", ModCold = "#86BCE2",
                   ModHeat = "#F3C623", NormExtrHeat = "#F06B18", ExtrHeat = "#AD1A2B")
stack_order <- c("ModCold", "ExtrCold", "ModHeat", "NormExtrHeat", "ExtrHeat")
ssp_labels <- c(`1` = "SSP1–2.6", `2` = "SSP2–4.5", `3` = "SSP3–7.0")
period_labels <- c(`2055` = "2055–59", `2095` = "2095–99")
source_bands <- c("65--69", "70--74", "75--79", "80--84",
                  "85--89", "90--94", "95--99", "100+")
age_levels <- c("65–69", "70–74", "75–79", "80–84", "85–89", "90–94", "95+")
age_map <- data.table(
  source_band = source_bands,
  age_band = c(age_levels[1:6], "95+", "95+"),
  agegroup = c(rep("65-74", 2), rep("75-84", 2), rep("85+", 4))
)
starts <- c(2055L, 2095L)

needed <- c("data/city_results.csv", "data/wittgenstein_pop.csv", "data/wittgenstein_assr.csv")
if (any(!file.exists(needed))) stop("Missing input: ", paste(needed[!file.exists(needed)], collapse = ", "))

message("[Figure 1] Preparing exact-band demography...")
meta <- fread("data/city_results.csv")
meta <- unique(meta[agegroup == "75-84", .(
  city = URAU_CODE, city_name = LABEL, country = CNTR_CODE,
  country_name = cntr_name, region, city_pop = agepop, city_death = death
)])
if (nrow(meta) != 854L || anyDuplicated(meta$city)) stop("Expected 854 city calibration rows.")

pop <- fread("data/wittgenstein_pop.csv")[
  ssp %in% ssp_set & age %chin% source_bands & sex %chin% c("Female", "Male"),
  .(ssp = as.integer(ssp), country = CNTR_CODE, source_band = as.character(age), sex,
    year_start = as.integer(year), pop_k = as.numeric(pop))]
assr <- fread("data/wittgenstein_assr.csv")[
  ssp %in% ssp_set & age %chin% source_bands & sex %chin% c("Female", "Male"),
  .(ssp = as.integer(ssp), country = CNTR_CODE, source_band = as.character(age), sex,
    year_start = as.integer(sub("^([0-9]{4}).*$", "\\1", period)), assr = as.numeric(assr))]
assr[ssp == 1L & assr == 1, assr := 0.9995]
country_band <- merge(pop, assr, by = c("ssp", "country", "source_band", "sex", "year_start"))
country_band <- merge(country_band, age_map, by = "source_band", all.x = TRUE)
if (!nrow(country_band) || anyNA(country_band[, .(pop_k, assr, agegroup)]))
  stop("The demographic source tables did not merge cleanly.")

# The three snapshots are each held for five years in the 2000–2014 calibration.
calibration <- country_band[
  year_start %in% c(2000L, 2005L, 2010L) & agegroup == "75-84",
  .(country_pop = sum(pop_k) * 1000,
    country_death = sum(pop_k * (1 - assr) / 5) * 1000),
  by = .(ssp, country, year_start)
][, .(country_pop = mean(country_pop), country_death = mean(country_death), n = .N),
  by = .(ssp, country)]
if (any(calibration$n != 3L) || any(calibration$country_death <= 0))
  stop("The 2000–2014 calibration snapshots are incomplete.")
factors <- merge(meta, calibration, by = "country", allow.cartesian = TRUE)
factors[, `:=`(popfac = city_pop / country_pop, dfac = city_death / country_death)]
if (nrow(factors) != 854L * length(ssp_set) || any(!is.finite(factors$popfac + factors$dfac)))
  stop("City calibration factors are incomplete.")

future <- country_band[year_start %in% starts,
  .(country_pop = sum(pop_k) * 1000,
    country_deaths_5y = sum(pop_k * (1 - assr)) * 1000),
  by = .(ssp, country, year_start, source_band, age_band, agegroup)]
bands <- merge(factors[, .(city, city_name, country, country_name, region, ssp, popfac, dfac)],
               future, by = c("ssp", "country"), allow.cartesian = TRUE)
bands[, `:=`(population = country_pop * popfac, deaths_5y = country_deaths_5y * dfac)]
bands <- bands[, .(population = sum(population), deaths_5y = sum(deaths_5y)),
  by = .(city, city_name, country, country_name, region, ssp, year_start, age_band, agegroup)]
if (any(bands$population <= 0) || any(bands$deaths_5y < 0)) stop("Invalid band demography.")
weights <- bands[, .(age_band, death_weight = deaths_5y / sum(deaths_5y)),
                 by = .(city, ssp, year_start, agegroup)]
weight_sums <- weights[, sum(death_weight), by = .(city, ssp, year_start, agegroup)]$V1
if (max(abs(weight_sums - 1)) > 1e-12) stop("Band death weights do not sum to one.")

message("[Figure 1] Reading 19-GCM ensemble attribution outputs...")
read_city <- function(ssp, city) {
  path <- file.path(batch_root, paste0("ssp", ssp), city, "ENSEMBLE", "01_attribution_grouped.csv")
  if (!file.exists(path)) stop("Missing attribution file: ", path)
  x <- fread(path)[
    year %in% c(2055:2059, 2095:2099) & branch %chin% c("with_cc", "without_cc") & range %chin% ranges,
    .(an = sum(an)),
    by = .(city = geo_id, branch, agegroup, range,
           year_start = fifelse(year <= 2059L, 2055L, 2095L))]
  x[, ssp := as.integer(ssp)]
  x
}
attribution <- rbindlist(lapply(ssp_set, function(s) {
  root <- file.path(batch_root, paste0("ssp", s))
  cities <- sort(list.dirs(root, recursive = FALSE, full.names = FALSE))
  cities <- cities[file.exists(file.path(root, cities, "ENSEMBLE", "01_attribution_grouped.csv"))]
  if (length(cities) != 854L) stop("Expected 854 ensemble city files for SSP", s, ".")
  rbindlist(lapply(cities, function(city) read_city(s, city)))
}))
expected <- length(ssp_set) * 854L * 2L * 3L * 5L * 2L
if (nrow(attribution) != expected || any(!is.finite(attribution$an))) stop("Incomplete attribution domain.")

allocated <- merge(attribution, weights, by = c("city", "ssp", "year_start", "agegroup"),
                   allow.cartesian = TRUE)
allocated[, band_an := an * death_weight]
allocation_check <- merge(
  attribution[, .(source_an = sum(an)), by = .(ssp, year_start, city, branch, agegroup, range)],
  allocated[, .(allocated_an = sum(band_an)), by = .(ssp, year_start, city, branch, agegroup, range)],
  by = c("ssp", "year_start", "city", "branch", "agegroup", "range"))
allocation_check[, absolute_difference := abs(source_an - allocated_an)]
if (max(allocation_check$absolute_difference) > 1e-10) stop("Band allocation changed attribution totals.")

specs <- rbind(
  data.table(level = "europe", id = "EUROPE", label = "Europe", region = NA_character_),
  data.table(level = "region", id = c("Northern", "Western", "Eastern", "Southern"),
             label = paste(c("Northern", "Western", "Eastern", "Southern"), "Europe"),
             region = c("Northern", "Western", "Eastern", "Southern")))
if (geo_level == "europe") specs <- specs[level == "europe"]
if (geo_level == "region") {
  specs <- specs[level == "region" & (!nzchar(geo_id_filter) | id == geo_id_filter)]
  if (!nrow(specs)) stop("Unknown region GEO_ID: ", geo_id_filter)
}

pool_geo <- function(spec) {
  cities <- if (spec$level == "europe") meta$city else meta[region == spec$region, city]
  base <- bands[city %chin% cities, .(population = sum(population), deaths_5y = sum(deaths_5y)),
                by = .(ssp, year_start, age_band)]
  a <- allocated[city %chin% cities, .(band_an = sum(band_an)),
                 by = .(ssp, year_start, age_band, branch, range)]
  climate <- a[, .(climate_adjustment = sum(band_an[branch == "with_cc"]) -
                                          sum(band_an[branch == "without_cc"])),
               by = .(ssp, year_start, age_band)]
  d <- merge(a, merge(base, climate, by = c("ssp", "year_start", "age_band")),
             by = c("ssp", "year_start", "age_band"))
  d[, total_deaths_5y := fifelse(branch == "with_cc", deaths_5y + climate_adjustment, deaths_5y)]
  d[, annual_probability := 1 - exp(-total_deaths_5y / (5 * population))]
  d[, risk := band_an / total_deaths_5y * annual_probability * 1e5]
  d[, old_rate := band_an / (5 * population) * 1e5]
  d[, `:=`(geo_level = spec$level, geo_id = spec$id, geography = spec$label)]
  d
}
plot_data <- rbindlist(lapply(seq_len(nrow(specs)), function(i) pool_geo(specs[i])))
if (any(!is.finite(plot_data$risk)) || any(plot_data$risk < 0)) stop("Invalid Figure 1 risks.")

# Full-domain totals must reproduce the ensemble pooled master inputs.
master_checks <- rbindlist(lapply(seq_len(nrow(specs)), function(i) {
  spec <- specs[i]
  rbindlist(lapply(ssp_set, function(s) {
    path <- file.path(batch_root, "geographies", paste0("ssp", s), spec$level, spec$id,
                      "03_master_table.csv")
    x <- fread(path)[year %in% c(2055:2059, 2095:2099) & range %chin% ranges,
      .(master_an = sum(an)), by = .(ssp, branch, year_start = fifelse(year <= 2059L, 2055L, 2095L))]
    y <- plot_data[geo_id == spec$id & ssp == s, .(figure_an = sum(band_an)),
                   by = .(ssp, branch, year_start)]
    merge(x, y, by = c("ssp", "branch", "year_start"))[
      , `:=`(geo_id = spec$id, absolute_difference = abs(master_an - figure_an))]
  }))
}))
if (max(master_checks$absolute_difference) > 1e-8) stop("Pooled attribution totals do not match master tables.")

# Exact source population versus the pipeline PCLM population, SSP3 ages 80–84.
population_checks <- rbindlist(lapply(meta$city, function(city_id) {
  path <- file.path(batch_root, "ssp3", city_id, "demography", "00_demography_single_age.csv")
  pclm <- fread(path)[year %in% starts & age %between% c(80L, 84L),
                      .(pclm_population = sum(pop)), by = .(city = geo_id, year_start = year)]
  exact <- bands[city == city_id & ssp == 3L & age_band == "80–84",
                 .(city, year_start, exact_population = population)]
  merge(exact, pclm, by = c("city", "year_start"))
}))
population_checks[, relative_difference_percent :=
                    100 * (exact_population - pclm_population) / pclm_population]
population_summary <- population_checks[, .(
  maximum_absolute_difference_percent = max(abs(relative_difference_percent))), by = year_start]
print(population_summary)
if (any(population_summary$maximum_absolute_difference_percent > 0.1))
  stop("Exact-band and PCLM population differ by more than 0.1%.")

# Required nonlinear pooling example.
ex_demo <- bands[ssp == 3L & year_start == 2095L & age_band == "85–89"]
ex_an <- allocated[ssp == 3L & year_start == 2095L & age_band == "85–89",
                   .(band_an = sum(band_an)), by = .(city, branch, range)]
ex_wide <- dcast(ex_an, city + range ~ branch, value.var = "band_an")
ex_climate <- ex_wide[, .(climate_adjustment = sum(with_cc - without_cc)), by = city]
ex_full <- ex_an[branch == "with_cc", .(band_an = sum(band_an)), by = city]
ex_city <- Reduce(function(x, y) merge(x, y, by = "city"),
                  list(ex_demo[, .(city, population, deaths_5y)], ex_climate, ex_full))
ex_city[, total_deaths_5y := deaths_5y + climate_adjustment]
ex_city[, city_risk := band_an / total_deaths_5y *
          (1 - exp(-total_deaths_5y / (5 * population))) * 1e5]
pooling_example <- data.table(
  pooled_risk = plot_data[geo_id == "EUROPE" & ssp == 3L & year_start == 2095L &
                            age_band == "85–89" & branch == "with_cc", sum(risk)],
  population_weighted_city_risk = weighted.mean(ex_city$city_risk, ex_city$population),
  unweighted_mean_city_risk = mean(ex_city$city_risk))
pooling_example[, `:=`(
  pooled_minus_population_weighted = pooled_risk - population_weighted_city_risk,
  pooled_minus_unweighted = pooled_risk - unweighted_mean_city_risk)]
print(pooling_example)
if (abs(pooling_example$pooled_minus_population_weighted) < 1e-12)
  stop("The pooling check did not demonstrate a nonlinear difference.")

fwrite(allocation_check, file.path(figure_root, "data", "figure1_allocation_checks.csv"))
fwrite(master_checks, file.path(figure_root, "data", "figure1_master_total_checks.csv"))
fwrite(population_checks, file.path(figure_root, "data", "figure1_population_checks.csv"))
fwrite(pooling_example, file.path(figure_root, "data", "figure1_pooling_example.csv"))

footnote <- paste(
  "Absolute risk is the annual probability of a temperature-attributable death among people alive at the start of the year, by 5-year age group; deaths at all non-optimal temperatures are included, not only those attributable to climate change.",
  "Historical 2.5th and 97.5th percentiles (P2.5, P97.5): ERA5-Land, 2000–2014. Projected P97.5: the historical P97.5 plus the change in the calibrated GCM P97.5 over the 15 preceding years.",
  "Without climate change, the normal extreme band is empty or nearly so. 95+ includes 100+. For SSP1, Wittgenstein survival ratios stored as exactly 1 were set to 0.9995. Central estimates, mean of 19 GCMs; no uncertainty intervals.",
  sep = "\n")

make_plot <- function(d) {
  label <- unique(d$geography)
  d <- copy(d)
  d[, `:=`(
    age_band = factor(age_band, levels = age_levels),
    scenario = factor(ssp_labels[as.character(ssp)], levels = ssp_labels[as.character(ssp_set)]),
    period = factor(period_labels[as.character(year_start)], levels = period_labels),
    range_stack = factor(range, levels = stack_order),
    direction = fifelse(range %chin% c("ExtrCold", "ModCold"), "cold", "heat"),
    y = as.numeric(factor(age_band, levels = age_levels)) + fifelse(branch == "with_cc", 0.16, -0.16))]
  setorder(d, scenario, period, branch, age_band, direction, range_stack)
  d[, `:=`(stack_end = cumsum(risk), stack_start = cumsum(risk) - risk),
    by = .(scenario, period, branch, age_band, direction)]
  d[direction == "heat", `:=`(xmin = stack_start, xmax = stack_end)]
  d[direction == "cold", `:=`(xmin = -stack_end, xmax = -stack_start)]
  half_width <- max(abs(c(d$xmin, d$xmax))) * 1.10
  left_period <- levels(d$period)[1]
  bar_labels <- unique(d[period == left_period & age_band == "95+",
                         .(scenario, period, branch, y)])
  bar_labels[, `:=`(x = -0.98 * half_width,
                    text = fifelse(branch == "with_cc", "with climate change", "without climate change"))]
  sides <- unique(d[, .(scenario, period)])
  sides <- rbind(sides[, .(scenario, period, x = -0.55 * half_width, y = 7.65, text = "Cold")],
                 sides[, .(scenario, period, x = 0.55 * half_width, y = 7.65, text = "Heat")])
  p <- ggplot(d) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - 0.135, ymax = y + 0.135,
                  fill = range_stack, alpha = branch), colour = NA) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#333333") +
    geom_text(data = bar_labels, aes(x, y, label = text), inherit.aes = FALSE,
              hjust = 0, size = 2.15, colour = "#333333") +
    geom_text(data = sides, aes(x, y, label = text), inherit.aes = FALSE,
              fontface = "bold", size = 2.8, colour = "#333333") +
    facet_grid(scenario ~ period, scales = "fixed", switch = "y") +
    scale_x_continuous(limits = c(-half_width, half_width),
      labels = function(z) format(abs(z), big.mark = ",", trim = TRUE), expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = seq_along(age_levels), labels = age_levels,
                       limits = c(0.6, 7.9), expand = expansion(mult = 0)) +
    scale_fill_manual(values = range_colours, labels = range_labels,
                      breaks = ranges, drop = FALSE, name = NULL) +
    scale_alpha_manual(values = c(with_cc = 1, without_cc = 0.32),
      breaks = c("with_cc", "without_cc"),
      labels = c("With climate change", "Without climate change"), name = NULL) +
    guides(fill = guide_legend(order = 1, nrow = 1),
           alpha = guide_legend(order = 2, override.aes = list(fill = "#555555"))) +
    labs(title = paste0(label, ": absolute risk of temperature-attributable death by age and temperature range"),
      x = "Annual risk of temperature-attributable death per 100,000 people (cold to the left; heat to the right)",
      y = "Age group", caption = footnote) +
    theme_minimal(base_size = 10.5) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      strip.placement = "outside", strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      plot.title = element_text(face = "bold", size = 15),
      plot.caption = element_text(size = 7.8, hjust = 0),
      legend.position = "bottom", legend.box = "vertical",
      panel.spacing = unit(0.75, "lines"))
  list(plot = p, half_width = half_width)
}

message("[Figure 1] Rendering figures...")
for (i in seq_len(nrow(specs))) {
  spec <- specs[i]
  d <- plot_data[geo_id == spec$id]
  built <- make_plot(d)
  out <- file.path(figure_root, spec$level, spec$id)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(out, "figure1_absolute_risk")
  ggsave(paste0(stem, ".png"), built$plot, width = 14, height = 14, units = "in", dpi = 300)
  ggsave(paste0(stem, ".pdf"), built$plot, width = 14, height = 14, units = "in", device = cairo_pdf)
  d[, `:=`(scenario = ssp_labels[as.character(ssp)], period = period_labels[as.character(year_start)])]
  fwrite(d, paste0(stem, "_data.csv"))
  cat(sprintf("%s shared horizontal half-width: %.6f\n", spec$label, built$half_width))
}

readme <- c(
  "# Paper figures and tables", "", "## Figure 1", "",
  paste("Figure 1 uses exact Wittgenstein five-year age bands at each period start.",
        "City population and deaths use separate Masselot calibration factors.",
        "Each five-year death total is spread over 1,825 modelled days.",
        "Attributable deaths are averaged over 19 GCMs before cities are pooled and annual risk is calculated.",
        "The 65–69 and 70–74 bands use the 65–74 exposure-response function.",
        "The 75–79 and 80–84 bands use 75–84. The 85–89, 90–94 and 95+ bands use 85+."),
  "", paste("SSP1–2.6 sets survival ratios stored as exactly 1 to 0.9995 as a sensitivity.",
            "SSP2–4.5 includes 12 Swiss cities with within-age-band mortality-shape warnings; the maximum decline is 11.09%."),
  "", "## Extreme-temperature thresholds", "",
  paste("The early Europe Version A maximum is set by Benidorm (ES044C), Spain.",
        "Clamp counts and annual shares are reported in the threshold data CSV files and the run console output.",
        "The figure does not include a clamp strip."))
writeLines(readme, file.path(figure_root, "README.md"), useBytes = TRUE)
message("Figure 1 complete. Outputs: ", figure_root)
