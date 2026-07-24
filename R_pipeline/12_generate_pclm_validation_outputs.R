#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sf)
})

message("\n[12] Generating PCLM validation outputs...")

dir.create("results/validation/pclm_validation", recursive = TRUE, showWarnings = FALSE)

input_single <- "results/le_li_input/lloyd_fig_s1_baseline_city_year_single_age.csv"
if (!file.exists(input_single)) {
  stop("Missing Lloyd-style single-age input table. Run R_pipeline/11_build_lloyd_fig_s1_baseline.R first.", call. = FALSE)
}

range_cols <- c("AN_ExtrCold", "AN_ModCold", "AN_ModHeat", "AN_ExtrHeat")
range_names <- sub("^AN_", "", range_cols)
age_groups_65plus <- c("65-74", "75-84", "85+")
example_cities <- c("AT001C", "ES001C", "UK001C", "BG001C", "SE001C", "IT001C")
example_year <- 2010L
profile_years <- 2000:2019

single_age <- fread(input_single)
single_age[, agegroup := fifelse(age <= 74L, "65-74", fifelse(age <= 84L, "75-84", "85+"))]
checks_file <- "results/le_li_input/lloyd_fig_s1_baseline_checks.csv"
if (!file.exists(checks_file)) {
  stop("Missing Lloyd-style additivity checks. Run R_pipeline/11_build_lloyd_fig_s1_baseline.R first.", call. = FALSE)
}
checks_exact <- fread(checks_file)

city_meta <- unique(fread("data/city_results.csv")[, .(
  URAU_CODE, LABEL, CNTR_CODE, cntr_name, region, lon, lat
)])

grouped_baseline <- rbindlist(lapply(list.files("temp_results_baseline", pattern = "\\.rds$", full.names = TRUE), function(f) {
  city_id <- sub("\\.rds$", "", basename(f))
  d <- as.data.table(readRDS(f))[gcm == "ERA5" & sim == 0 & agegroup %in% age_groups_65plus,
                                .(year, agegroup, range, an)]
  d[, URAU_CODE := city_id]
  d
}), use.names = TRUE)

grouped_wide <- dcast(grouped_baseline, URAU_CODE + year + agegroup ~ range, value.var = "an", fill = 0)
setnames(grouped_wide, range_names, paste0("grouped_", range_names))

single_grouped <- single_age[, .(
  AN_ExtrCold = sum(AN_ExtrCold),
  AN_ModCold = sum(AN_ModCold),
  AN_ModHeat = sum(AN_ModHeat),
  AN_ExtrHeat = sum(AN_ExtrHeat)
), by = .(URAU_CODE, year, agegroup)]
setnames(single_grouped, range_cols, paste0("single_", range_names))

validation <- merge(grouped_wide, single_grouped, by = c("URAU_CODE", "year", "agegroup"), all = TRUE)
validation <- merge(validation, city_meta, by = "URAU_CODE", all.x = TRUE)

for (rng in range_names) {
  validation[, (paste0("error_", rng)) := get(paste0("single_", rng)) - get(paste0("grouped_", rng))]
}

fwrite(validation, "results/validation/pclm_validation/pclm_grouped_vs_single_validation.csv")

message("Creating Figure S1: additivity validation check...")

fig1_data <- validation[URAU_CODE %in% example_cities & year == example_year]
fig1_long <- rbindlist(lapply(range_names, function(rng) {
  data.table(
    URAU_CODE = fig1_data$URAU_CODE,
    LABEL = fig1_data$LABEL,
    year = fig1_data$year,
    agegroup = fig1_data$agegroup,
    range = rng,
    source = "Grouped input",
    an = fig1_data[[paste0("grouped_", rng)]]
  )
}), use.names = TRUE)
fig1_long <- rbind(fig1_long,
                   rbindlist(lapply(range_names, function(rng) {
                     data.table(
                       URAU_CODE = fig1_data$URAU_CODE,
                       LABEL = fig1_data$LABEL,
                       year = fig1_data$year,
                       agegroup = fig1_data$agegroup,
                       range = rng,
                       source = "Sum of single-year ANs",
                       an = fig1_data[[paste0("single_", rng)]]
                     )
                   }), use.names = TRUE))
fig1_long[, LABEL := factor(LABEL, levels = city_meta[match(example_cities, URAU_CODE)]$LABEL)]
fig1_long[, range := factor(range, levels = range_names)]
fig1_long[, agegroup := factor(agegroup, levels = age_groups_65plus)]

p_add <- ggplot(fig1_long, aes(x = agegroup, y = an, fill = source)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  facet_grid(LABEL ~ range, scales = "free_y") +
  scale_fill_manual(values = c("Grouped input" = "#4c78a8", "Sum of single-year ANs" = "#f58518")) +
  labs(
    title = sprintf("PCLM Additivity Check in %d", example_year),
    subtitle = "Grouped 65+ baseline ANs are recovered exactly after single-year disaggregation",
    x = NULL,
    y = "Attributable deaths",
    fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")

ggsave("results/validation/pclm_validation/figS1_pclm_additivity_check.pdf", p_add, width = 12, height = 10)

message("Creating Figure S2: single-year age profiles by city...")

fig2_data <- single_age[URAU_CODE %in% example_cities & year %in% profile_years,
                         .(AN_ExtrCold = mean(AN_ExtrCold),
                           AN_ModCold = mean(AN_ModCold),
                           AN_ModHeat = mean(AN_ModHeat),
                           AN_ExtrHeat = mean(AN_ExtrHeat)),
                         by = .(URAU_CODE, LABEL, age)]
fig2_long <- melt(fig2_data, id.vars = c("URAU_CODE", "LABEL", "age"), variable.name = "range", value.name = "an")
fig2_long[, range := factor(sub("^AN_", "", range), levels = range_names)]
fig2_long[, LABEL := factor(LABEL, levels = city_meta[match(example_cities, URAU_CODE)]$LABEL)]

p_city_profile <- ggplot(fig2_long, aes(x = age, y = an, color = range)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ LABEL, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(ExtrCold = "#2c7fb8", ModCold = "#7fcdbb", ModHeat = "#fdae61", ExtrHeat = "#d7191c")) +
  labs(
    title = "Single-Year Age Profiles by City",
    subtitle = "Mean annual attributable deaths, averaged across 2000-2019",
    x = "Age",
    y = "Attributable deaths",
    color = "Range"
  ) +
  theme_minimal(base_size = 10)

ggsave("results/validation/pclm_validation/figS2_single_year_age_profiles_by_city.pdf", p_city_profile, width = 12, height = 8)

message("Creating Figure S3: age profiles by region...")

fig3_data <- single_age[year %in% profile_years,
                        .(AN_ExtrCold = mean(AN_ExtrCold),
                          AN_ModCold = mean(AN_ModCold),
                          AN_ModHeat = mean(AN_ModHeat),
                          AN_ExtrHeat = mean(AN_ExtrHeat)),
                        by = .(region, age)]
fig3_long <- melt(fig3_data, id.vars = c("region", "age"), variable.name = "range", value.name = "an")
fig3_long[, range := factor(sub("^AN_", "", range), levels = range_names)]
fig3_long[, region := factor(region, levels = c("Northern", "Western", "Eastern", "Southern"))]

p_region_profile <- ggplot(fig3_long, aes(x = age, y = an, color = range)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ region, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(ExtrCold = "#2c7fb8", ModCold = "#7fcdbb", ModHeat = "#fdae61", ExtrHeat = "#d7191c")) +
  labs(
    title = "Regional Single-Year Age Profiles",
    subtitle = "Mean annual attributable deaths averaged across cities and years 2000-2019",
    x = "Age",
    y = "Attributable deaths",
    color = "Range"
  ) +
  theme_minimal(base_size = 10)

ggsave("results/validation/pclm_validation/figS3_age_profiles_by_region.pdf", p_region_profile, width = 10, height = 8)

message("Creating Figure S4: negative ModHeat geography...")

modheat_city_year <- grouped_baseline[range == "ModHeat" & agegroup %in% age_groups_65plus,
                                      .(an_65plus = sum(an)), by = .(URAU_CODE, year)]
neg_map_data <- modheat_city_year[, .(
  pct_negative_years = mean(an_65plus < 0) * 100,
  n_negative_years = sum(an_65plus < 0)
), by = URAU_CODE]
neg_map_data <- merge(neg_map_data, city_meta, by = "URAU_CODE", all.x = TRUE)

map_bbox <- st_bbox(c(xmin = -12, xmax = 35, ymin = 35, ymax = 72), crs = st_crs(4326))
europe_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
europe_map <- st_crop(europe_map, map_bbox)

p_neg_map <- ggplot() +
  geom_sf(data = europe_map, fill = "grey96", color = "grey75", linewidth = 0.2) +
  geom_point(data = neg_map_data,
             aes(x = lon, y = lat, color = pct_negative_years, size = n_negative_years),
             alpha = 0.9) +
  scale_color_gradientn(colours = c("#f7fbff", "#6baed6", "#08306b"), name = "% years\nModHeat < 0") +
  scale_size_continuous(name = "Negative\nyears", range = c(0.5, 4)) +
  coord_sf(xlim = c(-12, 35), ylim = c(35, 72), expand = FALSE) +
  labs(
    title = "Geographic Distribution of Negative Moderate-Heat ANs",
    subtitle = "City color shows the share of years from 1990-2019 where total 65+ ModHeat AN is negative",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave("results/validation/pclm_validation/figS4_negative_an_geographic_distribution.pdf", p_neg_map, width = 9, height = 6.5)

message("Creating Table S1: PCLM disaggregation summary...")

error_long <- rbindlist(lapply(range_names, function(rng) {
  data.table(
    range = rng,
    abs_error = abs(checks_exact[[paste0("diff_", rng)]])
  )
}), use.names = TRUE)

negative_grouped <- grouped_baseline[, .(any_negative_grouped = any(an < 0)), by = .(URAU_CODE, year, range)]

mean_age_by_range <- rbindlist(lapply(range_cols, function(col) {
  positive_idx <- single_age[[col]] > 0
  data.table(
    range = sub("^AN_", "", col),
    mean_age_positive_an = weighted.mean(single_age$age[positive_idx], w = single_age[[col]][positive_idx])
  )
}), use.names = TRUE)

table_s1 <- error_long[, .(
  TotalCitiesProcessed = uniqueN(checks_exact$URAU_CODE),
  MaxAdditivityError = max(abs_error, na.rm = TRUE),
  MeanAdditivityError = mean(abs_error, na.rm = TRUE),
  MedianAdditivityError = median(abs_error, na.rm = TRUE)
), by = range]

neg_summary <- negative_grouped[any_negative_grouped == TRUE, .N, by = range]
setnames(neg_summary, "N", "NegativeGroupedCityYears")
table_s1 <- merge(table_s1, neg_summary, by = "range", all.x = TRUE)
table_s1[is.na(NegativeGroupedCityYears), NegativeGroupedCityYears := 0L]
table_s1 <- merge(table_s1, mean_age_by_range, by = "range", all.x = TRUE)

fwrite(table_s1, "results/validation/pclm_validation/tableS1_pclm_disaggregation_summary.csv")

summary_lines <- c(
  "# PCLM Validation Outputs",
  "",
  "## Inputs",
  "- Single-age input: `results/le_li_input/lloyd_fig_s1_baseline_city_year_single_age.csv`",
  "- Grouped baseline comparison source: `temp_results_baseline/*.rds` filtered to `gcm == \"ERA5\"` and `sim == 0`",
  "- Exemplary cities: AT001C, ES001C, UK001C, BG001C, SE001C, IT001C",
  sprintf("- Example year for Figure S1: %d", example_year),
  "- Averaging window for age-profile figures: 2000-2019",
  "",
  "## Outputs",
  "- figS1_pclm_additivity_check.pdf",
  "- figS2_single_year_age_profiles_by_city.pdf",
  "- figS3_age_profiles_by_region.pdf",
  "- figS4_negative_an_geographic_distribution.pdf",
  "- tableS1_pclm_disaggregation_summary.csv",
  "- pclm_grouped_vs_single_validation.csv"
)
writeLines(summary_lines, "results/validation/pclm_validation/README.md")

message("✅ PCLM validation outputs generated in results/validation/pclm_validation/")