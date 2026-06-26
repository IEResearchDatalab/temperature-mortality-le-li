################################################################################
#
# Validation and reporting for Madrid demo run
#   Generates summary tables and diagnostic plots for each pipeline step.
#
################################################################################

library(data.table)
library(ggplot2)
library(scales)

if (!exists("path_out")) path_out <- "../results/production"
plot_dir <- "../results/validation_plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Pipeline validation and reporting\n")
cat("========================================\n\n")

#---------------------------
# Step 1: Prep validation
#---------------------------

cat("\n--- Step 1: Raw input validation ---\n")
coefs <- readRDS(file.path(path_out, "coefs_all.rds"))
city_res <- readRDS(file.path(path_out, "city_res.rds"))
pop_eu <- readRDS(file.path(path_out, "pop_eu_raw_filtered.rds"))
lt_eu <- readRDS(file.path(path_out, "lt_eu_raw_filtered.rds"))

cat(sprintf("  ERF coefficients: %d rows, %d cities\n", nrow(coefs), uniqueN(coefs$URAU_CODE)))
cat(sprintf("  City metadata: %d rows, %d cities, %d countries\n",
    nrow(city_res), uniqueN(city_res$URAU_CODE), uniqueN(city_res$CNTR_CODE)))
cat(sprintf("  Eurostat population: %d rows, %d countries\n", nrow(pop_eu), uniqueN(pop_eu$country)))
cat(sprintf("  Eurostat life tables: %d rows, %d countries\n", nrow(lt_eu), uniqueN(lt_eu$country)))

#---------------------------
# Step 2: AN computation
#---------------------------

cat("\n--- Step 2: Climate-only ANs ---\n")
ans <- fread(file.path(path_out, "ans_annual_all_cities.csv"))
cat(sprintf("  AN rows: %s\n", format(nrow(ans), big.mark = ",")))
cat(sprintf("  Cities: %d | SSPs: %s | Temp ranges: %s\n",
    uniqueN(ans$city_code), paste(unique(ans$ssp), collapse = ","),
    paste(unique(ans$temp_range), collapse = ",")))
cat(sprintf("  Years: %d-%d\n", min(ans$year), max(ans$year)))
cat(sprintf("  Attribution type: %s\n", unique(ans$attribution_type)))

# AN summary by temp range
an_summary <- ans[, .(total_AN = sum(an_est)), by = .(temp_range, year)]
setorder(an_summary, temp_range, year)

p <- ggplot(an_summary, aes(x = year, y = total_AN, color = temp_range)) +
  geom_line() + geom_point() +
  labs(title = "Madrid: Annual ANs by temperature range",
       x = "Year", y = "Attributable deaths") +
  theme_minimal()
ggsave(file.path(plot_dir, "02_ans_by_range.png"), p, width = 10, height = 6)

cat("  Plot saved: 02_ans_by_range.png\n")

#---------------------------
# Step 3: Mortality backbone
#---------------------------

cat("\n--- Step 3: Mortality backbone ---\n")
backbone <- readRDS(file.path(path_out, "mortality_backbone.rds"))
cat(sprintf("  Backbone rows: %s\n", format(nrow(backbone), big.mark = ",")))
cat(sprintf("  Countries: %d | Ages: %d-%d | Years: %d-%d\n",
    uniqueN(backbone$country), min(backbone$age), max(backbone$age),
    min(backbone$year), max(backbone$year)))

# Plot mx over time for selected ages (Spain)
spain_mx <- backbone[country == "ES" & age %in% c(65, 75, 85)]
p <- ggplot(spain_mx, aes(x = year, y = mx_proj, color = factor(age))) +
  geom_line() +
  labs(title = "Spain: Projected mortality rates (mx) by age",
       x = "Year", y = "mx", color = "Age") +
  theme_minimal()
ggsave(file.path(plot_dir, "03_backbone_spain.png"), p, width = 10, height = 6)
cat("  Plot saved: 03_backbone_spain.png\n")

#---------------------------
# Step 4: City mortality schedules
#---------------------------

cat("\n--- Step 4: City mortality schedules ---\n")
city_mx <- readRDS(file.path(path_out, "city_mortality_schedules.rds"))
cat(sprintf("  City schedule rows: %s\n", format(nrow(city_mx), big.mark = ",")))
cat(sprintf("  Cities: %d | Ages: %d-%d | Years: %d-%d\n",
    uniqueN(city_mx$city_code), min(city_mx$age), max(city_mx$age),
    min(city_mx$year), max(city_mx$year)))

madrid_mx <- city_mx[city_code == "ES001C" & age %in% c(65, 75, 85)]
p <- ggplot(madrid_mx, aes(x = year, y = mx, color = factor(age))) +
  geom_line() +
  labs(title = "Madrid: Projected all-cause mortality rates",
       x = "Year", y = "mx", color = "Age") +
  theme_minimal()
ggsave(file.path(plot_dir, "04_madrid_mortality.png"), p, width = 10, height = 6)
cat("  Plot saved: 04_madrid_mortality.png\n")

#---------------------------
# Step 5: AN disaggregation
#---------------------------

cat("\n--- Step 5: AN disaggregation ---\n")
single_an <- fread(file.path(path_out, "ans_single_age_all_cities.csv"))
cat(sprintf("  Single-age AN rows: %s\n", format(nrow(single_an), big.mark = ",")))
cat(sprintf("  Cities: %d | Ages: %d-%d\n",
    uniqueN(single_an$city_code), min(single_an$age), max(single_an$age)))

# Check totals preserved
grouped_an <- ans[, .(total = sum(an_est)), by = .(city_code, year, temp_range)]
single_check <- single_an[, .(total_single = sum(AN)), by = .(city_code, year, temp_range)]
comparison <- merge(grouped_an, single_check, by = c("city_code", "year", "temp_range"))
comparison[, diff := total - total_single]
cat(sprintf("  Disaggregation max diff: %.6f\n", max(abs(comparison$diff))))

# Age profile of ANs (Madrid, year 2050, SSP370)
madrid_an <- single_an[city_code == "ES001C" & year == 2050 & ssp == "370"]
if (nrow(madrid_an) > 0) {
  p <- ggplot(madrid_an, aes(x = age, y = AN, fill = temp_range)) +
    geom_bar(stat = "identity") +
    labs(title = "Madrid 2050 SSP3-7.0: Single-age ANs by temp range",
         x = "Age", y = "Attributable deaths") +
    theme_minimal()
  ggsave(file.path(plot_dir, "05_madrid_an_age_profile.png"), p, width = 10, height = 6)
  cat("  Plot saved: 05_madrid_an_age_profile.png\n")
}

#---------------------------
# Step 6: Bookkeeping
#---------------------------

cat("\n--- Step 6: Bookkeeping ---\n")
analysis <- fread(file.path(path_out, "analysis_dataset_all_cities.csv"))
cat(sprintf("  Analysis dataset rows: %s\n", format(nrow(analysis), big.mark = ",")))
cat(sprintf("  Cities: %d | GCMs: %d | SSPs: %s\n",
    uniqueN(analysis$city_code), uniqueN(analysis$gcm),
    paste(unique(analysis$ssp), collapse = ",")))

# Check identity
analysis[, check := deaths - (residual_deaths + extreme_cold + moderate_cold + moderate_heat + extreme_heat)]
max_err <- max(abs(analysis$check))
cat(sprintf("  Identity check max error: %.10f\n", max_err))
cat(sprintf("  Negative residuals: %d\n", nrow(analysis[residual_deaths < 0])))

# Madrid summary
madrid_analysis <- analysis[city_code == "ES001C"]
madrid_summary <- madrid_analysis[,
  .(deaths = sum(deaths),
    residual = sum(residual_deaths),
    an_total = sum(extreme_cold + moderate_cold + moderate_heat + extreme_heat)),
  by = .(year, gcm, ssp)]

p <- ggplot(madrid_summary, aes(x = year)) +
  geom_line(aes(y = deaths, color = "All-cause")) +
  geom_line(aes(y = residual, color = "Non-temperature")) +
  geom_line(aes(y = an_total, color = "Temperature-attributable")) +
  facet_grid(gcm ~ ssp, scales = "free_y") +
  labs(title = "Madrid: Death decomposition", x = "Year", y = "Deaths") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(plot_dir, "06_madrid_decomposition.png"), p, width = 14, height = 10)
cat("  Plot saved: 06_madrid_decomposition.png\n")

#---------------------------
# Step 7: Life tables
#---------------------------

cat("\n--- Step 7: Life tables ---\n")
lt <- fread(file.path(path_out, "lifespan_inequality_all_cities.csv"))
cat(sprintf("  Life-table rows: %s\n", format(nrow(lt), big.mark = ",")))
cat(sprintf("  e65 range: [%.2f, %.2f]\n", min(lt$e65, na.rm = TRUE), max(lt$e65, na.rm = TRUE)))
cat(sprintf("  SD range: [%.2f, %.2f]\n", min(lt$sd, na.rm = TRUE), max(lt$sd, na.rm = TRUE)))

madrid_lt <- lt[city_code == "ES001C"]
p1 <- ggplot(madrid_lt, aes(x = year, y = e65, color = gcm)) +
  geom_line() + facet_wrap(~ssp) +
  labs(title = "Madrid: e65 over time", x = "Year", y = "Remaining LE at 65") +
  theme_minimal()
ggsave(file.path(plot_dir, "07_madrid_e65.png"), p1, width = 12, height = 8)

p2 <- ggplot(madrid_lt, aes(x = year, y = sd, color = gcm)) +
  geom_line() + facet_wrap(~ssp) +
  labs(title = "Madrid: Lifespan inequality (SD) above 65", x = "Year", y = "SD") +
  theme_minimal()
ggsave(file.path(plot_dir, "07_madrid_sd.png"), p2, width = 12, height = 8)
cat("  Plots saved: 07_madrid_e65.png, 07_madrid_sd.png\n")

#---------------------------
# Step 8: Decomposition
#---------------------------

cat("\n--- Step 8: Decomposition ---\n")
decomp_annual <- fread(file.path(path_out, "decomposition_annual_all_cities.csv"))
decomp_period <- fread(file.path(path_out, "decomposition_period_all_cities.csv"))
cat(sprintf("  Annual decomposition rows: %s\n", format(nrow(decomp_annual), big.mark = ",")))
cat(sprintf("  Period decomposition rows: %s\n", format(nrow(decomp_period), big.mark = ",")))

if (nrow(decomp_period) > 0) {
  madrid_decomp <- decomp_period[city_code == "ES001C"]
  
  # Sum of contributions vs total change
  check_decomp <- madrid_decomp[, .(sum_decomp = sum(delta_e65)), by = .(gcm, ssp, period)]
  cat(sprintf("  Decomposition check (sum delta_e65 by period):\n"))
  print(check_decomp)

  p <- ggplot(madrid_decomp, aes(x = age, y = delta_e65, fill = cause)) +
    geom_bar(stat = "identity") +
    facet_grid(period ~ gcm * ssp) +
    labs(title = "Madrid: Decomposition of e65 change by age and cause",
         x = "Age", y = "Contribution to e65 change") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(plot_dir, "08_madrid_decomp_e65.png"), p, width = 14, height = 10)
  cat("  Plot saved: 08_madrid_decomp_e65.png\n")
}

#---------------------------
# Step 9: Aggregation
#---------------------------

cat("\n--- Step 9: Aggregation ---\n")
country_analysis <- fread(file.path(path_out, "country_analysis.csv"))
europe_analysis <- fread(file.path(path_out, "europe_analysis.csv"))
country_lt <- fread(file.path(path_out, "country_lifespan_inequality.csv"))
europe_lt <- fread(file.path(path_out, "europe_lifespan_inequality.csv"))

cat(sprintf("  Country analysis: %s rows\n", format(nrow(country_analysis), big.mark = ",")))
cat(sprintf("  Europe analysis: %s rows\n", format(nrow(europe_analysis), big.mark = ",")))
cat(sprintf("  Country life tables: %s rows\n", format(nrow(country_lt), big.mark = ",")))
cat(sprintf("  Europe life tables: %s rows\n", format(nrow(europe_lt), big.mark = ",")))

#---------------------------
# Step 10: Absolute risk
#---------------------------

cat("\n--- Step 10: Absolute risk ---\n")
risk_summary <- fread(file.path(path_out, "absolute_risk_summary.csv"))
if (nrow(risk_summary) > 0) {
  cat(sprintf("  Risk summary rows: %s\n", format(nrow(risk_summary), big.mark = ",")))
  
  madrid_risk <- risk_summary[city_code == "ES001C"]
  if (nrow(madrid_risk) > 0) {
    cat(sprintf("  Madrid pct_change range: [%.2f, %.2f]\n",
        min(madrid_risk$pct_change, na.rm = TRUE),
        max(madrid_risk$pct_change, na.rm = TRUE)))
  }
}

#---------------------------
# Final summary
#---------------------------

cat("\n========================================\n")
cat("Validation complete\n")
cat("========================================\n")
cat(sprintf("Plots saved to: %s\n", normalizePath(plot_dir)))
cat(sprintf("Results in: %s\n", normalizePath(path_out)))