#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ungroup)
  library(parallel)
})

source("R_pipeline/functions/pclm_utils.R")

message("\n[11] Building Lloyd-style single-age baseline AN table...")

dir.create("results/le_li_input", recursive = TRUE, showWarnings = FALSE)

city_filter <- trimws(Sys.getenv("CITY_FILTER", unset = ""))
year_filter <- trimws(Sys.getenv("YEAR_FILTER", unset = ""))
range_levels <- c("ExtrCold", "ModCold", "ModHeat", "ExtrHeat")
age_groups_65plus <- c("65-74", "75-84", "85+")

city_meta <- fread("data/city_results.csv")
city_meta <- unique(city_meta[agegroup %in% age_groups_65plus, .(
  URAU_CODE,
  LABEL,
  CNTR_CODE,
  cntr_name,
  region,
  lon,
  lat,
  agegroup,
  age_start = fifelse(agegroup == "65-74", 65L, fifelse(agegroup == "75-84", 75L, 85L)),
  age_width = fifelse(agegroup == "85+", 16L, 10L),
  pop = agepop,
  death_baseline = death
)])

baseline_files <- list.files("temp_results_baseline", pattern = "\\.rds$", full.names = TRUE)
if (nzchar(city_filter)) {
  wanted <- trimws(strsplit(city_filter, ",", fixed = TRUE)[[1]])
  baseline_files <- baseline_files[sub("\\.rds$", "", basename(baseline_files)) %in% wanted]
}
if (!length(baseline_files)) {
  stop("No matching baseline files found in temp_results_baseline/", call. = FALSE)
}

process_city <- function(file_path) {
  city_id <- sub("\\.rds$", "", basename(file_path))
  d <- as.data.table(readRDS(file_path))[gcm == "ERA5" & sim == 0 & agegroup %in% age_groups_65plus]
  if (nzchar(year_filter)) {
    wanted_years <- as.integer(trimws(strsplit(year_filter, ",", fixed = TRUE)[[1]]))
    d <- d[year %in% wanted_years]
  }
  if (!nrow(d)) return(NULL)

  meta_city <- city_meta[URAU_CODE == city_id][order(age_start)]
  if (nrow(meta_city) != 3) stop(sprintf("Metadata for %s does not contain the three expected 65+ age groups.", city_id))

  x <- meta_city$age_start
  nlast <- tail(meta_city$age_width, 1)
  ages_single <- 65:(65 + pclm_expected_length(x, nlast) - 1)

  pop_single <- pclm_disaggregate_nonnegative(x = x, y = meta_city$pop, nlast = nlast)
  death_single <- pclm_disaggregate_nonnegative(x = x, y = meta_city$death_baseline, nlast = nlast)

  annual_grouped <- d[, .(an = sum(an)), by = .(year, agegroup, range)]
  annual_grouped <- dcast(annual_grouped, year + agegroup ~ range, value.var = "an", fill = 0)
  annual_grouped <- merge(annual_grouped, meta_city[, .(agegroup, age_start)], by = "agegroup")
  setorder(annual_grouped, year, age_start)

  years <- sort(unique(annual_grouped$year))
  rows <- vector("list", length(years))
  checks <- vector("list", length(years))

  for (i in seq_along(years)) {
    yr <- years[[i]]
    grp <- annual_grouped[year == yr][order(age_start)]
    single <- data.table(
      URAU_CODE = city_id,
      LABEL = meta_city$LABEL[1],
      CNTR_CODE = meta_city$CNTR_CODE[1],
      cntr_name = meta_city$cntr_name[1],
      region = meta_city$region[1],
      lon = meta_city$lon[1],
      lat = meta_city$lat[1],
      year = yr,
      age = ages_single,
      age_label = fifelse(ages_single == 100L, "100+", as.character(ages_single)),
      pop = pop_single,
      death_baseline = death_single,
      AN_ExtrCold = pclm_disaggregate_signed(x = x, y = grp$ExtrCold, nlast = nlast),
      AN_ModCold = pclm_disaggregate_signed(x = x, y = grp$ModCold, nlast = nlast),
      AN_ModHeat = pclm_disaggregate_signed(x = x, y = grp$ModHeat, nlast = nlast),
      AN_ExtrHeat = pclm_disaggregate_signed(x = x, y = grp$ExtrHeat, nlast = nlast)
    )
    single[, AN_total := AN_ExtrCold + AN_ModCold + AN_ModHeat + AN_ExtrHeat]
    rows[[i]] <- single

    checks[[i]] <- data.table(
      URAU_CODE = city_id,
      year = yr,
      diff_ExtrCold = sum(single$AN_ExtrCold) - sum(grp$ExtrCold),
      diff_ModCold = sum(single$AN_ModCold) - sum(grp$ModCold),
      diff_ModHeat = sum(single$AN_ModHeat) - sum(grp$ModHeat),
      diff_ExtrHeat = sum(single$AN_ExtrHeat) - sum(grp$ExtrHeat),
      any_na = anyNA(single[, .(AN_ExtrCold, AN_ModCold, AN_ModHeat, AN_ExtrHeat, pop, death_baseline)])
    )
  }

  list(data = rbindlist(rows), checks = rbindlist(checks))
}

num_cores <- min(8L, max(1L, detectCores() - 1L))
results <- mclapply(baseline_files, process_city, mc.cores = num_cores)
results <- Filter(Negate(is.null), results)

lloyd_s1 <- rbindlist(lapply(results, `[[`, "data"), use.names = TRUE)
checks <- rbindlist(lapply(results, `[[`, "checks"), use.names = TRUE)

output_file <- "results/le_li_input/lloyd_fig_s1_baseline_city_year_single_age.csv"
check_file <- "results/le_li_input/lloyd_fig_s1_baseline_checks.csv"

fwrite(lloyd_s1, output_file)
fwrite(checks, check_file)

message("Saved Lloyd-style baseline table to ", output_file)
message("Saved additivity checks to ", check_file)
message(sprintf("Rows: %d | Cities: %d | Years: %d | Any NA rows in checks: %s",
                nrow(lloyd_s1), uniqueN(lloyd_s1$URAU_CODE), uniqueN(lloyd_s1$year), any(checks$any_na)))