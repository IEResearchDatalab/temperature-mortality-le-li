################################################################################
#
# Validation: Step 07 — Life tables (e65 and LI)
#   Diagnose remaining life expectancy at 65 and lifespan inequality.
#
# Plots:
#   07a_e65_trajectory.png    — e65 over time by SSP and GCM
#   07b_sd_trajectory.png     — Lifespan inequality (SD) over time
#   07c_mx_vs_age.png         — Full mx profile at selected years
#   07d_life_table_columns.png— Key life-table columns for one year
#
################################################################################

library(ggplot2); library(scales); library(data.table)

validation_dir <- "../results/validation_plots"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Validation Step 07: Life tables (e65 and LI)\n")
cat("========================================\n\n")

lt_file <- file.path(path_out, "lifespan_inequality_all_cities.csv")
if (!file.exists(lt_file)) {
  cat("  WARNING: lifespan_inequality_all_cities.csv not found. Run 07 first.\n")
  quit(save = "no")
}

lt <- fread(lt_file)
cat(sprintf("  Life-table rows: %s\n", format(nrow(lt), big.mark = ",")))
cat(sprintf("  Cities: %d | GCMs: %d | SSPs: %s\n",
    uniqueN(lt$city_code), uniqueN(lt$gcm),
    paste(unique(lt$ssp), collapse = ", ")))
cat(sprintf("  e65 range: [%.2f, %.2f]\n", min(lt$e65, na.rm = TRUE),
    max(lt$e65, na.rm = TRUE)))
cat(sprintf("  SD range: [%.2f, %.2f]\n", min(lt$sd, na.rm = TRUE),
    max(lt$sd, na.rm = TRUE)))

#---------------------------
# 07a: e65 trajectory (Madrid)
#---------------------------

madrid_lt <- lt[city_code == "ES001C"]

p <- ggplot(madrid_lt, aes(x = year, y = e65, color = gcm)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ssp, labeller = labeller(ssp = ssp_labels), ncol = 3) +
  labs(title = "Madrid: Remaining life expectancy at 65 (e65)",
       subtitle = "Period life tables from full all-cause mortality schedules.",
       caption = "e65 = remaining years expected to live at age 65.
Higher e65 = improving longevity. Each line = one GCM.
SSP panels show different emissions scenarios.",
       x = "Year", y = "Remaining life expectancy at 65") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "07a_e65_trajectory.png"), p, width = 12, height = 7)
cat(sprintf("  07a: e65 trajectory saved\n"))

#---------------------------
# 07b: Lifespan inequality (SD)
#---------------------------

p <- ggplot(madrid_lt, aes(x = year, y = sd, color = gcm)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ssp, labeller = labeller(ssp = ssp_labels), ncol = 3) +
  labs(title = "Madrid: Lifespan inequality above 65 (SD of age at death)",
       subtitle = "SD = sqrt(Σ dx(x + ax − mean)² / lx_65). Higher = more unequal.",
       caption = "SD of attained age at death above 65.
Decreasing SD = compression of mortality (deaths concentrate at fewer ages).
Increasing SD = expansion of mortality.",
       x = "Year", y = "SD of age at death above 65") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "07b_sd_trajectory.png"), p, width = 12, height = 7)
cat(sprintf("  07b: SD trajectory saved\n"))

#---------------------------
# 07c: mx profile at selected years (Madrid, one GCM, SSP 2)
#---------------------------

analysis <- fread(file.path(path_out, "analysis_dataset_all_cities.csv"))
one_gcm <- unique(analysis$gcm)[1]
madrid_mx <- analysis[city_code == "ES001C" & gcm == one_gcm & ssp == "2" &
  year %in% c(2020, 2050, 2090)]

madrid_mx[, mx_total := deaths / population]

p <- ggplot(madrid_mx, aes(x = age, y = mx_total, color = factor(year))) +
  geom_line(linewidth = 1) +
  scale_y_log10(labels = comma) +
  labs(title = sprintf("Madrid SSP3-7.0: All-cause mx at selected years (%s)", one_gcm),
       subtitle = "Log scale. Age 20-100. Full mortality schedule including temperature effects.",
       caption = "mx = deaths/population. Each line = one year.
Shift downward over time = mortality improvement.
Kink at 85+ reflects open-age aggregation in Eurostat inputs.",
       x = "Age", y = "Mortality rate (mx)", color = "Year") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "07c_mx_vs_age.png"), p, width = 10, height = 6)
cat(sprintf("  07c: mx vs age saved\n"))

#---------------------------
# 07d: Life-table columns (selected year)
#---------------------------

# Reconstruct life table columns for Madrid 2020
source("../R/period_lifetable.R")
madrid_2020 <- analysis[city_code == "ES001C" & gcm == one_gcm & ssp == "2" & year == 2020]
mx_2020 <- madrid_2020$deaths / madrid_2020$population
mx_2020[is.nan(mx_2020) | is.infinite(mx_2020)] <- 1e-10
lt_2020 <- lifetable(mx_2020, madrid_2020$age)

lt_plot <- lt_2020[age >= 65]
lt_long <- melt(lt_plot, id.vars = "age",
  measure.vars = c("mx", "qx", "lx", "dx", "ex"),
  variable.name = "column", value.name = "value")

# Filter out mx for better scale
lt_long <- lt_long[column != "mx"]

p <- ggplot(lt_long, aes(x = age, y = value, color = column)) +
  geom_line(linewidth = 1) +
  facet_wrap(~column, scales = "free_y", ncol = 3) +
  labs(title = sprintf("Madrid 2020: Life-table columns (age 65+, %s)", one_gcm),
       subtitle = "Key life-table functions from period life table.",
       caption = "qx = death probability. lx = survivors (of 100k). dx = deaths.
ex = remaining life expectancy.
Built from all-cause mx including temperature effects.",
       x = "Age", y = "Value") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file.path(validation_dir, "07d_life_table_columns.png"), p, width = 12, height = 8)
cat(sprintf("  07d: Life-table columns saved\n"))

cat("\n--- Step 07 validation complete ---\n")