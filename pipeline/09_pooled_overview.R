#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# Pipeline Part 09: Cross-geography pooled summary
#   Combines the independently constructed European, regional and country
#   pooled life tables. Southern Europe is deliberately shown in orange, not
#   blue, following Simon Lloyd's presentation guidance.
#
################################################################################

source("pipeline/00_pkg_params.R")

root <- Sys.getenv("POOLED_ROOT", file.path("results/europe/pooled", paste0("ssp", ssp_name)))
if (!dir.exists(root)) stop("Missing pooled-results root: ", root, call. = FALSE)

meta <- unique(fread("data/city_results.csv")[, .(
  city = URAU_CODE, country = CNTR_CODE, country_name = cntr_name, region
)])
geo_def <- rbind(
  data.table(geo_level = "europe", geo_id = "EUROPE", label = "Europe"),
  unique(meta[, .(geo_level = "region", geo_id = region, label = paste(region, "Europe"))]),
  unique(meta[, .(geo_level = "country", geo_id = country, label = country_name)])
)
geo_def[, level_order := match(geo_level, c("europe", "region", "country"))]
setorder(geo_def, level_order, geo_id)
geo_def[, level_order := NULL]

read_one <- function(level, id, label) {
  d <- file.path(root, level, id)
  lev_file <- file.path(d, "07_le_li_levels.csv")
  bet_file <- file.path(d, "07_between_branch_decomposition.csv")
  chk_file <- file.path(d, "checks", "07_pooled_checks.csv")
  if (!all(file.exists(c(lev_file, bet_file, chk_file)))) {
    stop("Incomplete pooled output: ", d, call. = FALSE)
  }

  lev <- fread(lev_file)
  bet <- fread(bet_file)
  chk <- fread(chk_file)
  lev[, period_start := 2020L + ((year - 2020L) %/% 5L) * 5L]
  pl <- lev[, .(LE65 = mean(LE65), LI65 = mean(LI65)), by = .(branch, period_start)]
  first_p <- min(pl$period_start)
  last_p <- max(pl$period_start)
  gains <- pl[, .(
    gain_LE = LE65[period_start == last_p] - LE65[period_start == first_p],
    change_LI = LI65[period_start == last_p] - LI65[period_start == first_p]
  ), by = branch]
  last <- bet[as.integer(substr(period, 1, 4)) == last_p,
    .(dLE = sum(le_contribution), dLI = sum(li_contribution)), by = cause]
  total <- last[, .(cause = "total", dLE = sum(dLE), dLI = sum(dLI))]
  last <- rbind(last, total)
  n_cities <- as.integer(sub("/.*", "", chk[check_name == "city_coverage", value]))
  gain_without <- gains[branch == "without_cc", gain_LE]

  overall <- data.table(
    geo_level = level, geo_id = id, label = label, n_cities = n_cities,
    first_period = sprintf("%d-%d", first_p, first_p + 4L),
    last_period = sprintf("%d-%d", last_p, last_p + 4L),
    gain_LE_with_cc = gains[branch == "with_cc", gain_LE],
    gain_LE_without_cc = gain_without,
    change_LI_with_cc = gains[branch == "with_cc", change_LI],
    change_LI_without_cc = gains[branch == "without_cc", change_LI],
    endcentury_dLE = last[cause == "total", dLE],
    endcentury_dLE_pct_without_gain = 100 * last[cause == "total", dLE] / gain_without,
    endcentury_dLI = last[cause == "total", dLI]
  )
  last[, `:=`(
    geo_level = level, geo_id = id, label = label, n_cities = n_cities,
    period = sprintf("%d-%d", last_p, last_p + 4L),
    dLE_pct_without_gain = 100 * dLE / gain_without
  )]
  list(overall = overall, cause = last)
}

ans <- Map(read_one, geo_def$geo_level, geo_def$geo_id, geo_def$label)
summary_dt <- rbindlist(lapply(ans, `[[`, "overall"), fill = TRUE)
cause_dt <- rbindlist(lapply(ans, `[[`, "cause"), fill = TRUE)
fwrite(summary_dt, file.path(root, "09_geography_summary.csv"))
fwrite(cause_dt, file.path(root, "09_geography_cause_summary.csv"))

region_colors <- c(
  "Eastern Europe" = "#E6AB02",
  "Northern Europe" = "#1B9E77",
  "Southern Europe" = "#D95F02",
  "Western Europe" = "#7570B3"
)

reg <- summary_dt[geo_level == "region"]
reg_long <- melt(reg,
  id.vars = c("label", "n_cities"),
  measure.vars = c("endcentury_dLE", "endcentury_dLE_pct_without_gain"),
  variable.name = "measure", value.name = "value"
)
reg_long[, measure := factor(measure,
  levels = c("endcentury_dLE", "endcentury_dLE_pct_without_gain"),
  labels = c("Change in LE65 (years)", "% of without-CC LE65 gain"))]
p_reg <- ggplot(reg_long, aes(label, value, fill = label)) +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_col(width = 0.7) +
  facet_wrap(~measure, scales = "free_y") +
  scale_fill_manual(values = region_colors, guide = "none") +
  labs(
    title = "Regional climate-change effect on remaining life expectancy at 65",
    subtitle = sprintf("%s, 2095-2099, pooled life tables", ssplabs[ssp_name]),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(root, "09_region_endcentury.png"), p_reg, width = 10, height = 5.5, dpi = 160)

cty <- summary_dt[geo_level == "country"]
cty <- merge(cty, unique(meta[, .(geo_id = country, region)]), by = "geo_id", all.x = TRUE)
cty[, label := factor(label, levels = label[order(endcentury_dLE)])]
p_country <- ggplot(cty, aes(endcentury_dLE, label, fill = paste(region, "Europe"))) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_col(width = 0.75) +
  scale_fill_manual(values = region_colors, name = NULL) +
  labs(
    title = "Country climate-change effect on remaining life expectancy at 65",
    subtitle = sprintf("%s, 2095-2099, pooled life tables", ssplabs[ssp_name]),
    x = "With minus without climate change (years)", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(file.path(root, "09_country_endcentury_ranking.png"), p_country, width = 10, height = 8, dpi = 160)

print(summary_dt[geo_level %in% c("europe", "region")], digits = 5)
message("Saved cross-geography summary to ", root)
