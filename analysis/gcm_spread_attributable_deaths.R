suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript gcm_spread_attributable_deaths.R <repo_root> <ssp_root> <output_dir>")
}

repo_root <- normalizePath(args[[1]], mustWork = TRUE)
ssp_root <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta <- unique(
  fread(file.path(repo_root, "data", "city_results.csv"))[
    , .(city = URAU_CODE, region)
  ]
)
cities <- sort(meta$city)

source(file.path(repo_root, "pipeline", "00_pkg_params.R"), local = TRUE)
gcms <- gcmlist

expected <- CJ(city = cities, gcm = gcms)
expected[, file := file.path(ssp_root, city, gcm, "01_attribution_grouped.csv")]
missing <- expected[!file.exists(file)]
if (nrow(missing)) {
  fwrite(missing, file.path(output_dir, "missing_attribution_files.csv"))
  stop(sprintf("Missing %d city-GCM attribution files", nrow(missing)))
}
input_inventory <- copy(expected)
input_inventory[, `:=`(
  size = file.info(file)$size,
  mtime = as.character(file.info(file)$mtime)
)]
fwrite(input_inventory, file.path(output_dir, "input_inventory.csv"))

read_gcm <- function(gcm_name) {
  files <- expected[gcm == gcm_name]
  pieces <- lapply(seq_len(nrow(files)), function(i) {
    x <- fread(files$file[[i]], select = c("year", "branch", "agegroup", "range", "an", "geo_id"))
    x[, city := files$city[[i]]]
    x
  })
  x <- rbindlist(pieces, use.names = TRUE)
  x <- meta[x, on = "city"]
  x[, gcm := gcm_name]
  europe <- x[, .(an = sum(an)), by = .(gcm, year, branch, range)]
  europe[, geography := "Europe"]
  setcolorder(europe, c("gcm", "geography", "year", "branch", "range", "an"))
  regions <- x[, .(an = sum(an)), by = .(gcm, geography = paste(region, "Europe"), year, branch, range)]
  rbind(europe, regions)
}

message(sprintf("Reading %d files across %d GCMs", nrow(expected), length(gcms)))
ncores <- as.integer(Sys.getenv("NCORES", as.character(min(length(gcms), parallel::detectCores()))))
per_gcm_annual <- rbindlist(
  parallel::mclapply(gcms, read_gcm, mc.cores = ncores),
  use.names = TRUE
)
setorder(per_gcm_annual, geography, gcm, year, branch, range)

period_map <- function(year) {
  fifelse(year >= 2020 & year <= 2024, "2020-2024",
    fifelse(year >= 2095 & year <= 2099, "2095-2099", NA_character_)
  )
}
per_gcm_annual[, period := period_map(year)]
period <- per_gcm_annual[!is.na(period), .(annual_deaths = mean(an)),
  by = .(geography, gcm, period, branch, range)
]
wide <- dcast(period, geography + gcm + period + range ~ branch, value.var = "annual_deaths")
wide[, additional_deaths := with_cc - without_cc]
setorder(wide, geography, gcm, period, range)

net <- wide[, .(additional_deaths = sum(additional_deaths)),
  by = .(geography, gcm, period)
]
net[, range := "Net"]
gcm_period_cause <- rbind(
  wide[, .(geography, gcm, period, range, additional_deaths)],
  net,
  use.names = TRUE
)

spread <- gcm_period_cause[, .(
    n_gcms = .N,
    mean = mean(additional_deaths),
    median = median(additional_deaths),
    p02_5 = as.numeric(quantile(additional_deaths, 0.025, names = FALSE)),
    p05 = as.numeric(quantile(additional_deaths, 0.05, names = FALSE)),
    p95 = as.numeric(quantile(additional_deaths, 0.95, names = FALSE)),
    p97_5 = as.numeric(quantile(additional_deaths, 0.975, names = FALSE)),
    minimum = min(additional_deaths),
    maximum = max(additional_deaths),
    n_positive = sum(additional_deaths > 0),
    n_negative = sum(additional_deaths < 0),
    n_zero = sum(additional_deaths == 0),
    relative_sd = sd(additional_deaths) / abs(mean(additional_deaths))
  ),
  by = .(geography, period, range)
]
setorder(spread, period, geography, range)

fwrite(gcm_period_cause, file.path(output_dir, "gcm_period_cause_values.csv"))
fwrite(spread, file.path(output_dir, "gcm_spread_summary.csv"))

end_net <- gcm_period_cause[period == "2095-2099" & range == "Net"]
geo_order <- c("Europe", "Southern Europe", "Eastern Europe", "Western Europe", "Northern Europe")
end_net[, geography := factor(geography, levels = rev(geo_order))]

p1 <- ggplot(end_net, aes(x = additional_deaths, y = geography)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.5) +
  geom_boxplot(width = 0.52, outlier.shape = NA, fill = "#DCE6F1", colour = "#17365D") +
  geom_point(position = position_jitter(height = 0.10, width = 0), alpha = 0.65, size = 1.8, colour = "#17365D") +
  scale_x_continuous(labels = scales::label_comma()) +
  labs(
    title = "Climate-model spread in additional temperature-attributable deaths",
    subtitle = "SSP3-7.0, 2095-2099 annual mean, with climate change minus without climate change",
    x = "Additional deaths per year among ages 65+",
    y = NULL,
    caption = "Each point is one of 19 GCMs. Central ERF, 854 cities. Spread is descriptive, not a confidence interval."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0, colour = "grey35")
  )
ggsave(file.path(output_dir, "01_gcm_spread_net_deaths.png"), p1, width = 10, height = 5.8, dpi = 180)

cause_order <- c("ExtrHeat", "ModHeat", "ModCold", "ExtrCold")
cause_labels <- c(
  ExtrHeat = "Extreme heat", ModHeat = "Moderate heat",
  ModCold = "Moderate cold", ExtrCold = "Extreme cold"
)
eu_causes <- gcm_period_cause[
  geography == "Europe" & period == "2095-2099" & range %in% cause_order
]
eu_causes[, range := factor(range, levels = rev(cause_order), labels = rev(cause_labels[cause_order]))]

p2 <- ggplot(eu_causes, aes(x = additional_deaths, y = range)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.5) +
  geom_boxplot(width = 0.52, outlier.shape = NA, fill = "#FCE4D6", colour = "#A13D2D") +
  geom_point(position = position_jitter(height = 0.10, width = 0), alpha = 0.65, size = 1.8, colour = "#A13D2D") +
  scale_x_continuous(labels = scales::label_comma()) +
  labs(
    title = "Extreme heat dominates the SSP3 mortality signal across climate models",
    subtitle = "Europe, 2095-2099 annual mean, with climate change minus without climate change",
    x = "Change in attributable deaths per year among ages 65+",
    y = NULL,
    caption = "Each point is one of 19 GCMs. Negative cold values represent avoided cold-attributable deaths."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0, colour = "grey35")
  )
ggsave(file.path(output_dir, "02_gcm_spread_europe_causes.png"), p2, width = 10, height = 5.4, dpi = 180)

checks <- data.table(
  check = c("expected_city_gcm_files", "unique_gcms", "unique_cities", "finite_values"),
  value = c(nrow(expected), uniqueN(expected$gcm), uniqueN(expected$city), sum(is.finite(gcm_period_cause$additional_deaths))),
  expected = c(length(cities) * length(gcms), length(gcms), length(cities), nrow(gcm_period_cause))
)
checks[, status := fifelse(value == expected, "PASS", "FAIL")]
fwrite(checks, file.path(output_dir, "validation_checks.csv"))
if (any(checks$status == "FAIL")) stop("Validation failed")

git_value <- function(...) {
  value <- tryCatch(
    system2("git", c("-C", repo_root, ...), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  paste(value, collapse = " ")
}
manifest <- c(
  sprintf("generated_at=%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  sprintf("server=%s", Sys.info()[["nodename"]]),
  sprintf("repository=%s", repo_root),
  sprintf("branch=%s", git_value("branch", "--show-current")),
  sprintf("commit=%s", git_value("rev-parse", "HEAD")),
  sprintf("git_status=%s", ifelse(nzchar(git_value("status", "--porcelain")), "dirty", "clean")),
  "script=analysis/gcm_spread_attributable_deaths.R",
  sprintf("input_pattern=%s", file.path(ssp_root, "<CITY>", "<GCM>", "01_attribution_grouped.csv")),
  sprintf("output_directory=%s", normalizePath(output_dir, mustWork = TRUE)),
  sprintf("command=Rscript analysis/gcm_spread_attributable_deaths.R %s %s %s", repo_root, ssp_root, output_dir),
  sprintf("cities=%d", length(cities)),
  sprintf("gcms=%d", length(gcms)),
  "validation=PASS"
)
writeLines(manifest, file.path(output_dir, "run_manifest.txt"))

message("Completed GCM-spread analysis: ", output_dir)
