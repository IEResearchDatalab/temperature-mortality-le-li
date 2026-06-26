################################################################################
#
# Validation: Step 04 — City mortality schedules
#   Diagnose the translation from national backbone to city schedules.
#
# Plots:
#   04a_city_vs_national.png — Madrid mx vs Spain national mx
#   04b_city_deaths.png      — Madrid projected deaths over time
#   04c_age_pyramid.png      — Madrid population age structure (selected years)
#
################################################################################

library(ggplot2); library(scales); library(data.table); library(ungroup)

validation_dir <- "../results/validation_plots"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Validation Step 04: City mortality schedules\n")
cat("========================================\n\n")

sched_file <- file.path(path_out, "city_mortality_schedules.rds")
if (!file.exists(sched_file)) {
  cat("  WARNING: city_mortality_schedules.rds not found. Run 04 first.\n")
  quit(save = "no")
}

city_mx <- readRDS(sched_file)
backbone <- readRDS(file.path(path_out, "mortality_backbone.rds"))

cat(sprintf("  City schedule rows: %s\n", format(nrow(city_mx), big.mark = ",")))
cat(sprintf("  Cities: %s\n", paste(unique(city_mx$city_code), collapse = ", ")))
cat(sprintf("  Countries: %s\n", paste(unique(city_mx$country), collapse = ", ")))
cat(sprintf("  Years: %d-%d\n", min(city_mx$year), max(city_mx$year)))

#---------------------------
# 04a: Madrid vs Spain national mx
#---------------------------

madrid <- city_mx[city_code == "ES001C" & age %in% c(65, 75, 85)]
spain_nat <- backbone[country == "ES" & age %in% c(65, 75, 85)]

p <- ggplot() +
  geom_line(data = spain_nat, aes(x = year, y = mx_proj, color = factor(age),
    linetype = "National"), linewidth = 1) +
  geom_line(data = madrid, aes(x = year, y = mx, color = factor(age),
    linetype = "Madrid"), linewidth = 1) +
  scale_y_log10(labels = comma) +
  labs(title = "Madrid vs Spain: Projected mortality rates",
       subtitle = "Madrid inherits national improvement ratios applied to city baseline.",
       caption = "City method: city_baseline_mx(age) × (national_projected / national_baseline).
National = solid, Madrid = dashed. Log scale.",
       x = "Year", y = "Mortality rate (mx)", color = "Age") +
  scale_linetype_manual(values = c("National" = "solid", "Madrid" = "dashed")) +
  guides(linetype = guide_legend(title = "Source")) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "04a_city_vs_national.png"), p, width = 10, height = 6)
cat(sprintf("  04a: City vs national saved\n"))

#---------------------------
# 04b: Madrid projected deaths over time
#---------------------------

madrid_deaths <- city_mx[city_code == "ES001C",
  .(total_deaths = sum(deaths)), by = year]

madrid_by_age <- city_mx[city_code == "ES001C" & age %% 10 == 5,
  .(deaths = sum(deaths)), by = .(year, age)]

p <- ggplot(madrid_by_age[age %in% c(65, 75, 85, 95)],
    aes(x = year, y = deaths, color = factor(age))) +
  geom_line(linewidth = 1) +
  labs(title = "Madrid: Projected all-cause deaths by age",
       subtitle = "Derived from city mortality schedules (mx × population).",
       caption = "Deaths = mx × population for each single age, summed to selected ages.
Shows the combined effect of projected mx improvement and population aging.",
       x = "Year", y = "Annual deaths", color = "Age") +
  scale_y_continuous(labels = comma) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "04b_city_deaths.png"), p, width = 10, height = 6)
cat(sprintf("  04b: City deaths saved\n"))

#---------------------------
# 04c: Population age structure (Madrid)
#---------------------------

madrid_pop <- city_mx[city_code == "ES001C" & year %in% c(2020, 2050, 2080),
  .(total_pop = sum(population)), by = .(year, age)]

p <- ggplot(madrid_pop, aes(x = age, y = total_pop, fill = factor(year))) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  labs(title = "Madrid: Projected population age structure",
       subtitle = "Disaggregated from city baseline using PCLM. Fixed across scenarios.",
       caption = "Population is projected by applying city baseline age structure.
Changes over time reflect the mortality improvement ratios applied to population counts.",
       x = "Age", y = "Population", fill = "Year") +
  scale_y_continuous(labels = comma) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(validation_dir, "04c_age_pyramid.png"), p, width = 10, height = 6)
cat(sprintf("  04c: Age pyramid saved\n"))

cat("\n--- Step 04 validation complete ---\n")