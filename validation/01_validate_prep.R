################################################################################
#
# Validation: Step 01 — Raw input preparation
#   Diagnose the loaded ERF coefficients, city metadata, and Eurostat data.
#
# Plots:
#   01a_erf_curves.png   — RR vs temperature for each age group (Madrid)
#   01b_city_map.png     — City locations (lat/lon)
#   01c_baseline_age.png — Baseline deaths and population by age (Madrid)
#   01d_hist_mx.png      — Historical mx trends for Spain (selected ages)
#
################################################################################

library(ggplot2); library(scales); library(data.table)
library(dlnm)

validation_dir <- "../results/validation_plots"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Validation Step 01: Raw input preparation\n")
cat("========================================\n\n")

#---------------------------
# Load data from Step 01 outputs
#---------------------------

coefs <- readRDS(file.path(path_out, "coefs_all.rds"))
city_meta <- readRDS(file.path(path_out, "city_res.rds"))
lt_hist <- readRDS(file.path(path_out, "lt_eu_raw_filtered.rds"))
pop_proj <- readRDS(file.path(path_out, "pop_eu_raw_filtered.rds"))

cat(sprintf("  ERF coefficients: %d rows x %d cities\n",
    nrow(coefs), uniqueN(coefs$URAU_CODE)))
cat(sprintf("  City metadata: %d rows x %d cities x %d countries\n",
    nrow(city_meta), uniqueN(city_meta$URAU_CODE), uniqueN(city_meta$CNTR_CODE)))
cat(sprintf("  Life table (hist): %d rows x %d countries\n",
    nrow(lt_hist), uniqueN(lt_hist$country)))
cat(sprintf("  Population (proj): %d rows x %d countries\n",
    nrow(pop_proj), uniqueN(pop_proj$country)))

#---------------------------
# 01a: ERF curves (RR vs temperature by age group)
#---------------------------

# Build basis from a representative temperature range (e.g., Madrid's)
madrid_coefs <- coefs[URAU_CODE == "ES001C"]
temp_seq <- seq(-10, 40, by = 0.5)
argvar <- list(fun = varfun, degree = vardegree,
  knots = quantile(temp_seq, varper / 100), Bound = range(temp_seq))
basis <- do.call(onebasis, c(list(x = temp_seq), argvar))

rr_data <- rbindlist(lapply(unique(madrid_coefs$agegroup), function(ag) {
  coef_row <- madrid_coefs[agegroup == ag]
  coef_vec <- as.numeric(coef_row[, .(b1, b2, b3, b4, b5)])
  log_rr <- as.numeric(basis %*% coef_vec)
  # Find MMT within 25th-99th percentile of typical temps
  search_idx <- which(temp_seq >= quantile(temp_seq, 0.25) &
    temp_seq <= quantile(temp_seq, 0.99))
  mmt <- temp_seq[search_idx][which.min(log_rr[search_idx])]
  cenvec <- do.call(onebasis, c(list(x = mmt), argvar))
  rr <- as.numeric(exp(log_rr - drop(cenvec %*% coef_vec)))
  data.table(temp = temp_seq, agegroup = ag, RR = rr, MMT = mmt)
}))

p <- ggplot(rr_data, aes(x = temp, y = RR, color = agegroup)) +
  geom_line(linewidth = 1) +
  geom_vline(data = unique(rr_data[, .(agegroup, MMT)]),
    aes(xintercept = MMT, color = agegroup), linetype = "dashed", alpha = 0.5) +
  labs(title = "Madrid: Exposure-response curves by age group (Masselot ERFs)",
       subtitle = paste0("B-spline basis, degree=", vardegree,
         ", knots at ", paste(varper, collapse = "/"), " percentiles"),
       caption = "Dashed vertical lines = MMT per age group. RR centered at MMT.",
       x = "Temperature (°C)", y = "Relative Risk") +
  scale_y_log10() + theme_minimal() + theme(legend.position = "bottom") +
  geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.4)

ggsave(file.path(validation_dir, "01a_erf_curves.png"), p, width = 10, height = 6)
cat(sprintf("  01a: ERF curves saved (%s)\n", "01a_erf_curves.png"))

#---------------------------
# 01b: City map
#---------------------------

p <- ggplot(city_meta, aes(x = lon, y = lat, color = CNTR_CODE)) +
  geom_point(alpha = 0.6, size = 1) +
  labs(title = sprintf("City locations (%d cities across Europe)",
       uniqueN(city_meta$URAU_CODE)),
       caption = "Each point = one city. Color = country code.",
       x = "Longitude", y = "Latitude") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file.path(validation_dir, "01b_city_map.png"), p, width = 10, height = 7)
cat(sprintf("  01b: City map saved (%s)\n", "01b_city_map.png"))

#---------------------------
# 01c: Baseline age distribution (Madrid)
#---------------------------

madrid_demo <- city_meta[URAU_CODE == "ES001C"]
p <- ggplot(madrid_demo, aes(x = agegroup)) +
  geom_bar(aes(y = agepop, fill = "Population"), stat = "identity", alpha = 0.7) +
  geom_bar(aes(y = death, fill = "Deaths"), stat = "identity", alpha = 0.7) +
  labs(title = "Madrid: Baseline population and deaths by age group",
       caption = "Bars show total population (blue) and annual deaths (red) in each age group.",
       x = "Age group", y = "Count") +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("Population" = "#3498db", "Deaths" = "#e74c3c")) +
  theme_minimal() + theme(legend.position = "bottom", legend.title = element_blank())

ggsave(file.path(validation_dir, "01c_baseline_age.png"), p, width = 8, height = 5)
cat(sprintf("  01c: Baseline age distribution saved (%s)\n", "01c_baseline_age.png"))

#---------------------------
# 01d: Historical mx trends (Spain)
#---------------------------

spain_mx <- lt_hist[country == "ES" & age %in% c(30, 50, 65, 75, 85)]
p <- ggplot(spain_mx, aes(x = year, y = mx, color = factor(age))) +
  geom_line(linewidth = 1) + geom_point(size = 0.5) +
  scale_y_log10(labels = comma) +
  labs(title = "Spain: Historical mortality rates by age (Eurostat demo_mlifetable)",
       subtitle = paste0("Years ", min(spain_mx$year), "-", max(spain_mx$year)),
       caption = "Log scale. Eurostat combined sex (T) death rates.",
       x = "Year", y = "Mortality rate (mx)", color = "Age") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "01d_hist_mx.png"), p, width = 10, height = 6)
cat(sprintf("  01d: Historical mx saved (%s)\n", "01d_hist_mx.png"))

#---------------------------
# Summary
#---------------------------
cat(sprintf("\n  Age groups in ERF: %s\n", paste(unique(coefs$agegroup), collapse = ", ")))
cat(sprintf("  Countries with data: %s\n", paste(sort(unique(lt_hist$country)), collapse = ", ")))
cat(sprintf("  Life-table years: %d-%d\n", min(lt_hist$year), max(lt_hist$year)))
cat(sprintf("  Projection years: %d-%d\n", min(pop_proj$year), max(pop_proj$year)))
cat("\n--- Step 01 validation complete ---\n")