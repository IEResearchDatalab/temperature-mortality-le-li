#!/usr/bin/env Rscript

################################################################################
#
# Temperature-related mortality and its impact on life expectancy and
# lifespan inequality at older ages in European cities
#
# Pipeline Part 08: Figures and headline results for pooled geographies
#   Follows the presentation decisions from Simon Lloyd on 25 Sep 2026:
#     - smooth LE65 / LI65+ trajectories with a roughly 10-year moving average;
#     - show with- versus without-CC trajectories and their gap;
#     - retain age and temperature-range decompositions;
#     - show cumulative contributions to LE65 change since 2020-2024.
#
################################################################################

source("pipeline/00_pkg_params.R")

levels_file <- file.path(out_dir, "07_le_li_levels.csv")
between_file <- file.path(out_dir, "07_between_branch_decomposition.csv")
within_file <- file.path(out_dir, "07_within_branch_period_decomposition.csv")
for (f in c(levels_file, between_file, within_file)) {
  if (!file.exists(f)) stop("Missing Part 07 output: ", f, call. = FALSE)
}

levels_dt <- fread(levels_file)
between <- fread(between_file)
within <- fread(within_file)
geo_label <- unique(levels_dt$label)
if (length(geo_label) != 1L) stop("Part 07 output must contain one geography.", call. = FALSE)

message(sprintf("\n[08] Building pooled figures for %s...", geo_label))

time_blocks <- c(2020, 2040, 2060, 2080, 2100)
age_band_breaks <- c(seq(65, 100, 5), Inf)
snapshot_periods <- list("2045-2054" = c(2045, 2054), "2085-2094" = c(2085, 2094))

between[, period_start := as.integer(substr(period, 1, 4))]
within[, `:=`(
  p_from = as.integer(substr(period_from, 1, 4)),
  p_to = as.integer(substr(period_to, 1, 4))
)]

#----- Fig 1: smoothed trajectories and the with - without CC gap

lev_long <- melt(levels_dt, id.vars = c("branch", "year"), measure.vars = c("LE65", "LI65"),
  variable.name = "measure", value.name = "raw_value")
setorder(lev_long, branch, measure, year)
lev_long[, value := frollmean(raw_value, n = 10L, align = "center"), by = .(branch, measure)]
lev_long[, measure := factor(measure, levels = c("LE65", "LI65"),
  labels = c("Remaining life expectancy at 65 (years)", "Lifespan inequality 65+ (SD, years)"))]
fwrite(lev_long, file.path(out_dir, "08_smoothed_levels.csv"))

p1a <- ggplot(lev_long[!is.na(value)], aes(year, value, colour = branch)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~measure, scales = "free_y") +
  scale_colour_manual(values = branch_colors, labels = branch_labels, name = NULL) +
  labs(x = NULL, y = NULL, caption = "Centred 10-year moving average; incomplete edge windows omitted") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

gap <- between[, .(dLE = sum(le_contribution), dLI = sum(li_contribution)), by = period_start]
gap[, mid := period_start + 2]
li_scale <- max(abs(gap$dLE)) / max(abs(gap$dLI), 1e-12)
p1b <- ggplot(gap, aes(mid)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(aes(y = dLE, colour = "Delta LE65 (left axis)"), linewidth = 0.9) +
  geom_point(aes(y = dLE, colour = "Delta LE65 (left axis)")) +
  geom_line(aes(y = dLI * li_scale, colour = "Delta LI65+ (right axis)"), linewidth = 0.9, linetype = 2) +
  geom_point(aes(y = dLI * li_scale, colour = "Delta LI65+ (right axis)")) +
  scale_y_continuous(name = "Delta LE65, with - without CC (years)",
    sec.axis = sec_axis(~ . / li_scale, name = "Delta LI65+, with - without CC (SD)")) +
  scale_colour_manual(values = c("Delta LE65 (left axis)" = "black", "Delta LI65+ (right axis)" = "#7b3294"), name = NULL) +
  labs(x = "5-year period (midpoint)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p1 <- (p1a / p1b) + plot_annotation(
  title = sprintf("%s: LE65 and LI65+, with vs without climate change", geo_label),
  subtitle = sprintf("%s, pooled life table, central ERF estimates", ssplabs[ssp_name])
)
ggsave(file.path(fig_dir, "08_fig1_smoothed_trajectories.png"), p1, width = 11, height = 9, dpi = 160)

#----- Fig 2: within-branch contributions by age band and time block

wp <- within[cause %in% range_levels]
wp[, block_id := findInterval(p_to, time_blocks)]
block_lab <- wp[, .(lab = sprintf("%s to %s", period_from[which.min(p_from)], period_to[which.max(p_to)])), by = block_id]
wp <- merge(wp, block_lab, by = "block_id")
wp[, age_band := cut(age, age_band_breaks, right = FALSE,
  labels = c(paste(head(age_band_breaks, -2), head(age_band_breaks, -2) + 4, sep = "-"), "100+"))]
blk <- rbind(
  wp[, .(contribution = sum(le_contribution), measure = "LE65 (years)"), by = .(branch, block = lab, age_band, cause)],
  wp[, .(contribution = sum(li_contribution), measure = "LI65+ (SD)"), by = .(branch, block = lab, age_band, cause)]
)
blk[, `:=`(
  cause = factor(cause, levels = range_levels),
  block = factor(block, levels = block_lab[order(block_id)]$lab)
)]
p2 <- ggplot(blk, aes(age_band, contribution, fill = cause)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, colour = "grey40") +
  facet_grid(measure ~ branch + block, scales = "free_y", labeller = labeller(branch = branch_labels)) +
  scale_fill_manual(values = range_colors, labels = range_labels, name = NULL) +
  labs(
    title = sprintf("%s: temperature-related contributions to changes in LE65 and LI65+", geo_label),
    subtitle = sprintf("%s, pooled life table; changes between 5-year-period means", ssplabs[ssp_name]),
    x = "Age group", y = "Contribution to change"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")
ggsave(file.path(fig_dir, "08_fig2_within_branch_by_age_block.png"), p2, width = 16, height = 7, dpi = 160)

#----- Fig 3: age profile of the climate-change effect

snap <- rbindlist(lapply(names(snapshot_periods), function(n) {
  pr <- snapshot_periods[[n]]
  between[period_start >= pr[1] & period_start + 4 <= pr[2] & cause %in% range_levels,
    .(LE = mean(le_contribution), LI = mean(li_contribution)), by = .(age, cause)][, snapshot := n]
}))
snap <- melt(snap, id.vars = c("age", "cause", "snapshot"), variable.name = "measure")
snap[, `:=`(
  cause = factor(cause, levels = range_levels),
  measure = factor(measure, levels = c("LE", "LI"), labels = c("Delta LE65 (years)", "Delta LI65+ (SD)"))
)]
p3 <- ggplot(snap, aes(age, value, colour = cause, linetype = snapshot)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~measure, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = range_colors, labels = range_labels, name = NULL) +
  labs(
    title = sprintf("%s: age profile of the climate-change effect", geo_label),
    subtitle = sprintf("%s, pooled life table; with minus without climate change", ssplabs[ssp_name]),
    x = "Age", y = "Contribution by single year of age", linetype = NULL
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "08_fig3_age_profile_cc_effect.png"), p3, width = 10, height = 8, dpi = 160)

#----- Fig 4: cumulative climate-change contribution to the LE65 gain

inc <- within[cause %in% range_levels,
  .(contribution = sum(le_contribution)), by = .(branch, period_to, p_to, cause)]
inc <- dcast(inc, period_to + p_to + cause ~ branch, value.var = "contribution")
inc[, increment := with_cc - without_cc]
setorder(inc, cause, p_to)
inc[, cumulative_contribution := cumsum(increment), by = cause]
base <- data.table(period_to = "2020-2024", p_to = 2020L, cause = range_levels,
  with_cc = 0, without_cc = 0, increment = 0, cumulative_contribution = 0)
cum <- rbind(base, inc, use.names = TRUE)
total <- cum[, .(
  with_cc = sum(with_cc), without_cc = sum(without_cc), increment = sum(increment),
  cumulative_contribution = sum(cumulative_contribution)
), by = .(period_to, p_to)][, cause := "total"]
cum_plot <- rbind(cum, total, use.names = TRUE)
cum_plot[, cause := factor(cause, levels = c(range_levels, "total"))]
fwrite(cum_plot, file.path(out_dir, "08_cumulative_cc_le_contributions.csv"))

p4 <- ggplot(cum_plot, aes(p_to + 2, cumulative_contribution, colour = cause)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.4) +
  scale_colour_manual(
    values = c(range_colors, total = "black"),
    labels = c(range_labels, total = "Total"), name = NULL
  ) +
  labs(
    title = sprintf("%s: cumulative climate-change contribution to the LE65 gain", geo_label),
    subtitle = sprintf("%s, pooled life table; with minus without climate change, cumulative since 2020-2024", ssplabs[ssp_name]),
    x = "5-year period (midpoint)", y = "Cumulative contribution to the LE65 gain (years)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "08_fig4_cumulative_le_contributions.png"), p4, width = 11, height = 6, dpi = 160)

#----- Headline numbers use the same 5-year mean mortality schedules as the
# decomposition, not means of annual LE/LI values (the nonlinear life-table
# transformation makes those quantities slightly different).

levels_dt[, period_start := 2020L + ((year - 2020L) %/% 5L) * 5L]
first_p <- min(levels_dt$period_start); last_p <- max(levels_dt$period_start)
gain <- within[, .(
  gain_LE = sum(le_contribution),
  change_LI = sum(li_contribution)
), by = branch]
cc_last <- between[period_start == last_p, .(dLE = sum(le_contribution), dLI = sum(li_contribution)), by = cause]
cc_first <- between[period_start == first_p, .(dLE_first = sum(le_contribution), dLI_first = sum(li_contribution)), by = cause]
cc_change <- merge(cc_last, cc_first, by = "cause", sort = FALSE)
cc_change[, `:=`(dLE_change = dLE - dLE_first, dLI_change = dLI - dLI_first)]
cc_change[, cause_order := match(cause, cause_levels)]
setorder(cc_change, cause_order)
cc_change[, cause_order := NULL]
gain_wo <- gain[branch == "without_cc", gain_LE]
summary_dt <- rbind(
  data.table(item = sprintf("LE65 gain %d-%d to %d-%d, %s", first_p, first_p + 4, last_p, last_p + 4, gain$branch),
    value = gain$gain_LE, unit = "years"),
  data.table(item = sprintf("LI65+ change %d-%d to %d-%d, %s", first_p, first_p + 4, last_p, last_p + 4, gain$branch),
    value = gain$change_LI, unit = "SD"),
  data.table(item = sprintf("CC effect on LE65, %d-%d, %s", last_p, last_p + 4, c(cc_last$cause, "total")),
    value = c(cc_last$dLE, sum(cc_last$dLE)), unit = "years"),
  data.table(item = sprintf("CC effect on LE65 as %% of the without-CC gain, %d-%d, %s", last_p, last_p + 4, c(cc_last$cause, "total")),
    value = 100 * c(cc_last$dLE, sum(cc_last$dLE)) / gain_wo, unit = "%"),
  data.table(item = sprintf("CC effect on LI65+, %d-%d, %s", last_p, last_p + 4, c(cc_last$cause, "total")),
    value = c(cc_last$dLI, sum(cc_last$dLI)), unit = "SD"),
  data.table(item = sprintf("CC-induced change in LE65 gain, %d-%d to %d-%d, %s", first_p, first_p + 4, last_p, last_p + 4,
      c(cc_change$cause, "total")),
    value = c(cc_change$dLE_change, sum(cc_change$dLE_change)), unit = "years"),
  data.table(item = sprintf("CC-induced change in LE65 gain as %% of the without-CC gain, %d-%d to %d-%d, %s",
      first_p, first_p + 4, last_p, last_p + 4, c(cc_change$cause, "total")),
    value = 100 * c(cc_change$dLE_change, sum(cc_change$dLE_change)) / gain_wo, unit = "%"),
  data.table(item = sprintf("CC-induced change in LI65+, %d-%d to %d-%d, %s", first_p, first_p + 4, last_p, last_p + 4,
      c(cc_change$cause, "total")),
    value = c(cc_change$dLI_change, sum(cc_change$dLI_change)), unit = "SD")
)
fwrite(summary_dt, file.path(out_dir, "08_summary.csv"))
print(summary_dt, digits = 5)
message("Saved pooled figures to ", fig_dir)

