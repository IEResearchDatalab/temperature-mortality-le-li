#!/usr/bin/env Rscript

################################################################################
# Diagnostic requested by Simon Lloyd (25 Sep 2026): separate changes in the
# climate-change LE65 gap into (a) the evolving baseline mortality schedule and
# (b) the evolving temperature-attributable fractions. Uses an exact two-factor
# Shapley decomposition between consecutive five-year periods.
################################################################################

source("pipeline/00_pkg_params.R")

master_file <- Sys.getenv("MASTER_FILE",
  file.path("results/europe", paste0("ssp", ssp_name), city_id, "ENSEMBLE", "03_master_table.csv"))
diag_dir <- Sys.getenv("DIAG_DIR", file.path("results/europe/diagnostics", city_id))
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(master_file)) stop("Missing city master table: ", master_file, call. = FALSE)

x <- fread(master_file)
x[, period_start := 2020L + ((year - 2020L) %/% 5L) * 5L]
x[, mx := deaths_component / pop]
pm <- x[, .(mx = mean(mx)), by = .(branch, period_start, age, range)]

base <- pm[branch == "without_cc", .(baseline_mx = sum(mx)), by = .(period_start, age)]
wide <- dcast(pm, period_start + age + range ~ branch, value.var = "mx")
wide[, gap_mx := with_cc - without_cc]
gap <- wide[, .(gap_mx = sum(gap_mx)), by = .(period_start, age)]
sched <- merge(base, gap, by = c("period_start", "age"))
sched[, gap_fraction := gap_mx / baseline_mx]
if (any(!is.finite(unlist(sched[, .(baseline_mx, gap_mx, gap_fraction)]))) || any(sched$baseline_mx <= 0)) {
  stop("Invalid baseline or climate-gap mortality schedule.", call. = FALSE)
}

le65 <- function(mx) {
  px <- exp(-mx * nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(mx == 0, lx * nx, dx / mx)
  sum(Lx)
}
climate_gap <- function(baseline_mx, gap_fraction) {
  le65(baseline_mx * (1 + gap_fraction)) - le65(baseline_mx)
}
get_sched <- function(p, col) sched[period_start == p][order(age), get(col)]

periods <- sort(unique(sched$period_start))
levels <- rbindlist(lapply(periods, function(p) {
  b <- get_sched(p, "baseline_mx"); g <- get_sched(p, "gap_fraction")
  data.table(period_start = p, climate_gap_LE65 = climate_gap(b, g),
    baseline_LE65 = le65(b), with_cc_LE65 = le65(b * (1 + g)))
}))

steps <- rbindlist(lapply(seq_len(length(periods) - 1L), function(i) {
  p0 <- periods[i]; p1 <- periods[i + 1L]
  b0 <- get_sched(p0, "baseline_mx"); b1 <- get_sched(p1, "baseline_mx")
  g0 <- get_sched(p0, "gap_fraction"); g1 <- get_sched(p1, "gap_fraction")
  f00 <- climate_gap(b0, g0); f10 <- climate_gap(b1, g0)
  f01 <- climate_gap(b0, g1); f11 <- climate_gap(b1, g1)
  baseline_component <- 0.5 * ((f10 - f00) + (f11 - f01))
  temperature_component <- 0.5 * ((f01 - f00) + (f11 - f10))
  data.table(
    period_from = p0, period_to = p1,
    climate_gap_from = f00, climate_gap_to = f11,
    gap_change = f11 - f00,
    baseline_mortality_component = baseline_component,
    temperature_fraction_component = temperature_component,
    closure_error = baseline_component + temperature_component - (f11 - f00)
  )
}))
steps[, `:=`(
  cumulative_baseline_component = cumsum(baseline_mortality_component),
  cumulative_temperature_component = cumsum(temperature_fraction_component),
  cumulative_gap_change = cumsum(gap_change)
)]

if (max(abs(steps$closure_error)) > 1e-12 ||
    abs(tail(steps$cumulative_gap_change, 1) - (tail(levels$climate_gap_LE65, 1) - levels$climate_gap_LE65[1])) > 1e-12) {
  stop("Two-factor diagnostic failed closure.", call. = FALSE)
}

fwrite(levels, file.path(diag_dir, "demography_temperature_gap_levels.csv"))
fwrite(steps, file.path(diag_dir, "demography_temperature_gap_decomposition.csv"))
fwrite(data.table(
  check_name = c("step_closure", "whole_period_closure"), status = "PASS",
  value = c(sprintf("max %.3e", max(abs(steps$closure_error))),
    sprintf("error %.3e", tail(steps$cumulative_gap_change, 1) - (tail(levels$climate_gap_LE65, 1) - levels$climate_gap_LE65[1])))
), file.path(diag_dir, "demography_temperature_gap_checks.csv"))

plot_dt <- melt(steps,
  id.vars = c("period_to"),
  measure.vars = c("cumulative_baseline_component", "cumulative_temperature_component", "cumulative_gap_change"),
  variable.name = "component", value.name = "value"
)
plot_dt[, component := factor(component,
  levels = c("cumulative_baseline_component", "cumulative_temperature_component", "cumulative_gap_change"),
  labels = c("Baseline mortality schedule", "Temperature-attributable fractions", "Total change in climate LE65 gap"))]
p <- ggplot(plot_dt, aes(period_to + 2, value, colour = component)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c("#7570B3", "#D95F02", "black"), name = NULL) +
  labs(
    title = sprintf("%s: why the climate-change LE65 gap evolves", city_name),
    subtitle = sprintf("%s; exact two-factor Shapley decomposition since 2020-2024", ssplabs[ssp_name]),
    x = "5-year period (midpoint)", y = "Cumulative contribution (years)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(diag_dir, "demography_temperature_gap_decomposition.png"), p, width = 10, height = 6, dpi = 160)

print(levels, digits = 5)
print(steps[, .(
  total_baseline_component = sum(baseline_mortality_component),
  total_temperature_component = sum(temperature_fraction_component),
  total_gap_change = sum(gap_change),
  max_abs_closure_error = max(abs(closure_error))
)], digits = 6)
message("Saved demography/temperature gap diagnostic to ", diag_dir)
