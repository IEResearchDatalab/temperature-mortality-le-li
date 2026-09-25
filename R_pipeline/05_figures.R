#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# R Pipeline Part 05: Summary figures and headline numbers
#   Built only from Part 04 outputs. The figures Simon asked for:
#     Fig 1  LE65 and LI65+ trajectories, with vs without CC, plus the gap
#            (Delta LE65 and Delta LI65+ on a dual axis)      [10 Sep §3.3, §8]
#     Fig 2  Contributions of each temperature range to the change in LE65 and
#            LI65+ within each branch, by 5-year age band and ~20-year block
#            (changes between consecutive 5-year-period means)
#            (Lloyd et al. 2024 Figs 3-4 layout)              [26 Aug, 7 Sep]
#     Fig 3  Age profile of the climate-change effect (with - without CC) by
#            single age and temperature range, ~2050 and ~2090 [10 Sep §8]
#   05_summary.csv: headline numbers, including the climate-change loss as a
#   share of the LE65 gain without CC (relative change, Simon 8 Sep).
#
################################################################################

source("R_pipeline/00_pkg_params.R")

message("\n[05] Building summary figures...")

time_blocks <- c(2020, 2040, 2060, 2080, 2100)
age_band_breaks <- c(seq(65, 100, 5), Inf)
snapshot_periods <- list("2045-2054" = c(2045, 2054), "2085-2094" = c(2085, 2094))

levels_dt <- fread(file.path(out_dir, "04_le_li_levels.csv"))
between <- fread(file.path(out_dir, "04_between_branch_decomposition.csv"))
between[, period_start := as.integer(substr(period, 1, 4))]

#----- Fig 1: trajectories and the with - without CC gap

lev_long <- melt(levels_dt, id.vars = c("branch", "year"), measure.vars = c("LE65", "LI65"), variable.name = "measure")
lev_long[, measure := factor(measure, levels = c("LE65", "LI65"),
  labels = c("Remaining life expectancy at 65 (years)", "Lifespan inequality 65+ (SD, years)"))]
p1a <- ggplot(lev_long, aes(year, value, colour = branch)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~measure, scales = "free_y") +
  scale_colour_manual(values = branch_colors, labels = branch_labels, name = NULL) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

gap <- between[, .(dLE = sum(le_contribution), dLI = sum(li_contribution)), by = .(period_start)]
gap[, mid := period_start + 2]
li_scale <- max(abs(gap$dLE)) / max(abs(gap$dLI), 1e-12)
p1b <- ggplot(gap, aes(mid)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(aes(y = dLE, colour = "Delta LE65 (left axis)"), linewidth = 0.8) +
  geom_point(aes(y = dLE, colour = "Delta LE65 (left axis)")) +
  geom_line(aes(y = dLI * li_scale, colour = "Delta LI65+ (right axis)"), linewidth = 0.8, linetype = 2) +
  geom_point(aes(y = dLI * li_scale, colour = "Delta LI65+ (right axis)")) +
  scale_y_continuous(name = "Delta LE65, with - without CC (years)",
    sec.axis = sec_axis(~ . / li_scale, name = "Delta LI65+, with - without CC (SD)")) +
  scale_colour_manual(values = c("Delta LE65 (left axis)" = "black", "Delta LI65+ (right axis)" = "#7b3294"), name = NULL) +
  labs(x = "5-year period (midpoint)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
p1 <- (p1a / p1b) + plot_annotation(title = sprintf("%s: LE65 and LI65+, with vs without climate change", city_name), subtitle = scenario_label)
ggsave(file.path(fig_dir, "05_fig1_trajectories.png"), p1, width = 11, height = 9, dpi = 160)

#----- Fig 2: within-branch contributions by age band and time block
# Uses the decomposition between consecutive 5-year-period means (Part 04), so
# block sums are changes between period means and are not driven by the
# weather of single endpoint years.

wp <- fread(file.path(out_dir, "04_within_branch_period_decomposition.csv"))
wp[, `:=`(p_from = as.integer(substr(period_from, 1, 4)), p_to = as.integer(substr(period_to, 1, 4)))]
wp <- wp[cause %in% range_levels]
wp[, block_id := findInterval(p_to, time_blocks)]
block_lab <- wp[, .(lab = sprintf("%s to %s", period_from[which.min(p_from)], period_to[which.max(p_to)])), by = block_id]
wp <- merge(wp, block_lab, by = "block_id")
wp[, age_band := cut(age, age_band_breaks, right = FALSE,
  labels = c(paste(head(age_band_breaks, -2), head(age_band_breaks, -2) + 4, sep = "-"), "100+"))]
blk <- rbind(
  wp[, .(contribution = sum(le_contribution), measure = "LE65 (years)"), by = .(branch, block = lab, age_band, cause)],
  wp[, .(contribution = sum(li_contribution), measure = "LI65+ (SD)"), by = .(branch, block = lab, age_band, cause)]
)
blk[, cause := factor(cause, levels = range_levels)]
blk[, block := factor(block, levels = block_lab[order(block_id)]$lab)]
p2 <- ggplot(blk, aes(age_band, contribution, fill = cause)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, colour = "grey40") +
  facet_grid(measure ~ branch + block, scales = "free_y", labeller = labeller(branch = branch_labels)) +
  scale_fill_manual(values = range_colors, labels = range_labels, name = NULL) +
  labs(title = sprintf("%s: contribution of temperature-related mortality to the change in LE65 and LI65+", city_name),
       subtitle = paste(scenario_label, "- change between 5-year-period means, summed over ~20-year blocks; all other causes not shown"),
       x = "Age group", y = "Contribution to change") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")
ggsave(file.path(fig_dir, "05_fig2_within_branch_by_age_block.png"), p2, width = 16, height = 7, dpi = 160)

#----- Fig 3: age profile of the climate-change effect

snap <- rbindlist(lapply(names(snapshot_periods), function(n) {
  pr <- snapshot_periods[[n]]
  between[period_start >= pr[1] & period_start + 4 <= pr[2] & cause %in% range_levels,
    .(LE = mean(le_contribution), LI = mean(li_contribution)), by = .(age, cause)][, snapshot := n]
}))
snap <- melt(snap, id.vars = c("age", "cause", "snapshot"), variable.name = "measure")
snap[, `:=`(cause = factor(cause, levels = range_levels),
  measure = factor(measure, levels = c("LE", "LI"), labels = c("Delta LE65 (years)", "Delta LI65+ (SD)")))]
p3 <- ggplot(snap, aes(age, value, colour = cause, linetype = snapshot)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~measure, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = range_colors, labels = range_labels, name = NULL) +
  labs(title = sprintf("%s: age profile of the climate-change effect (with - without CC)", city_name),
       subtitle = paste(scenario_label, "- mean of 5-year-period decompositions"),
       x = "Age", y = "Contribution by single year of age", linetype = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "05_fig3_age_profile_cc_effect.png"), p3, width = 10, height = 8, dpi = 160)

#----- Headline numbers

first_y <- min(levels_dt$year); last_y <- max(levels_dt$year)
gain <- levels_dt[, .(gain_LE = LE65[year == last_y] - LE65[year == first_y],
                      change_LI = LI65[year == last_y] - LI65[year == first_y]), by = branch]
last_p <- max(between$period_start)
cc_last <- between[period_start == last_p, .(dLE = sum(le_contribution), dLI = sum(li_contribution)), by = cause]
gain_wo <- gain[branch == "without_cc", gain_LE]
summary_dt <- rbind(
  data.table(item = sprintf("LE65 gain %d-%d, %s", first_y, last_y, gain$branch), value = gain$gain_LE, unit = "years"),
  data.table(item = sprintf("LI65+ change %d-%d, %s", first_y, last_y, gain$branch), value = gain$change_LI, unit = "SD"),
  data.table(item = sprintf("CC effect on LE65, %d-%d, %s", last_p, last_p + 4, c(cc_last$cause, "total")),
             value = c(cc_last$dLE, sum(cc_last$dLE)), unit = "years"),
  data.table(item = sprintf("CC effect on LE65 as %% of the without-CC gain, %d-%d, %s", last_p, last_p + 4, c(cc_last$cause, "total")),
             value = 100 * c(cc_last$dLE, sum(cc_last$dLE)) / gain_wo, unit = "%"),
  data.table(item = sprintf("CC effect on LI65+, %d-%d, %s", last_p, last_p + 4, c(cc_last$cause, "total")),
             value = c(cc_last$dLI, sum(cc_last$dLI)), unit = "SD")
)
fwrite(summary_dt, file.path(out_dir, "05_summary.csv"))
print(summary_dt, digits = 4)
message("Saved figures to ", fig_dir, " and summary to ", file.path(out_dir, "05_summary.csv"))
