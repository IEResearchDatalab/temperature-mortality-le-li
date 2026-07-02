################################################################################
# 
# Contrasting future heat and cold-related mortality under climate change, 
# demographic and adaptation scenarios in 854 European cities
#
# R Code Part 5: Generate Tables (Masselot 2023 & Lloyd 2024 Reproductions)
#
# Adapted from Pierre Masselot & Antonio Gasparrini
#
################################################################################

library(data.table)

source("R_pipeline/01_initialize.R")

# Create output directories
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

message("\n[5/5] Generating Summary Tables...")

city_meta <- fread("data/city_results.csv")
city_meta <- city_meta[, .(URAU_CODE, cntr_name, agegroup, pop = agepop, deaths_baseline = death)]

rds_files <- list.files("temp_results", pattern = "\\.rds$", full.names = TRUE)

file_info <- data.table(
  path = rds_files, 
  URAU_CODE = gsub(".rds", "", basename(rds_files))
)
file_info <- merge(file_info, unique(city_meta[, .(URAU_CODE, cntr_name)]), by = "URAU_CODE")

message("Processing ", nrow(file_info), " cities in batches...")

process_city_agg <- function(f, city_meta_full) {
  city_id <- gsub(".rds", "", basename(f))
  d <- try(readRDS(f), silent = TRUE)
  if(inherits(d, "try-error")) return(NULL)
  
  meta_v <- city_meta_full[URAU_CODE == city_id]
  if(nrow(meta_v) == 0) return(NULL)
  cn <- meta_v$cntr_name[1]
  
  d[, decade := (year %/% 10) * 10]
  d_agg <- d[, .(an = sum(an)), by = .(decade, range, ssp, gcm, agegroup, sim)]
  
  d_agg[, pop_baseline := meta_v$pop[match(agegroup, meta_v$agegroup)]]
  d_agg[, deaths_baseline := meta_v$deaths_baseline[match(agegroup, meta_v$agegroup)]]
  d_agg[, country := cn]
  
  return(d_agg)
}

# Aggregate by chunks to file
dir.create("temp_table_batches", showWarnings = FALSE)
batch_size <- 50
n_batches <- ceiling(nrow(file_info) / batch_size)

for(b in 1:n_batches) {
  start_idx <- (b-1) * batch_size + 1
  end_idx <- min(b * batch_size, nrow(file_info))
  
  message("  Batch ", b, "/", n_batches, " (Cities ", start_idx, "-", end_idx, ")...")
  
  chunk_list <- list()
  for(i in start_idx:end_idx) {
    res <- process_city_agg(file_info$path[i], city_meta)
    if(!is.null(res)) chunk_list[[length(chunk_list) + 1]] <- res
  }
  
  chunk_dt <- rbindlist(chunk_list)
  # Collapse to country level within the chunk to save intermediate space
  chunk_agg <- chunk_dt[, .(
    an = sum(an),
    pop_baseline = sum(pop_baseline),
    deaths_baseline = sum(deaths_baseline)
  ), by = .(country, decade, range, ssp, gcm, agegroup, sim)]
  
  saveRDS(chunk_agg, paste0("temp_table_batches/batch_", b, ".rds"))
  rm(chunk_list, chunk_dt, chunk_agg); gc(verbose = FALSE)
}

message("Merging country-level batches...")
batch_files <- list.files("temp_table_batches", full.names = TRUE)
batch_list <- lapply(batch_files, readRDS)
country_results <- rbindlist(batch_list)[, .(
  an = sum(an),
  pop_baseline = sum(pop_baseline),
  deaths_baseline = sum(deaths_baseline)
), by = .(country, decade, range, ssp, gcm, agegroup, sim)]

rm(batch_list); gc()

# Save granular data for LE/LI script
message("Saving input for LE/LI decomposition...")
le_li_data <- country_results[, .(
  an = mean(an),
  pop = mean(pop_baseline),
  death_baseline = mean(deaths_baseline)
), by = .(country, decade, range, ssp, gcm, agegroup)]

setnames(le_li_data, "country", "cntr_name")
dir.create("results/le_li_input", recursive = TRUE, showWarnings = FALSE)
fwrite(le_li_data, "results/le_li_input/le_li_input_ans.csv")

message("Generating Table S1/S6...")
t1_raw <- country_results[, .(
  an_total = sum(an),
  pop_total = sum(pop_baseline),
  deaths_total = sum(deaths_baseline)
), by = .(country, decade, ssp, range, gcm, sim)]

t1_raw[, type := ifelse(grepl("Cold", range), "Cold", "Heat")]
t1_combined <- t1_raw[, .(
  an = sum(an_total),
  pop = mean(pop_total),
  deaths_base = mean(deaths_total)
), by = .(country, decade, ssp, type, gcm, sim)]

table_s6 <- t1_combined[, .(
  AN = mean(an),
  AN_low = quantile(an, 0.025),
  AN_hi = quantile(an, 0.975),
  AF = mean(an / (deaths_base + 1e-6)) * 100,
  AF_low = quantile(an / (deaths_base + 1e-6), 0.025) * 100,
  AF_hi = quantile(an / (deaths_base + 1e-6), 0.975) * 100,
  Rate = mean(an / (pop + 1e-6)) * 100000,
  Rate_low = quantile(an / (pop + 1e-6), 0.025) * 100000,
  Rate_hi = quantile(an / (pop + 1e-6), 0.975) * 100000
), by = .(country, decade, ssp, type)]

fwrite(table_s6, "tables/Table_Masselot2023_S6_Country.csv")
fwrite(table_s6, "results/tables/Table_Masselot2023_S6_Country.csv")

table_s1 <- t1_raw[, .(
  an_total = sum(an_total)
), by = .(decade, ssp, range, gcm, sim)]

table_s1_summary <- table_s1[, .(
  AN = mean(an_total),
  low = quantile(an_total, 0.025),
  hi = quantile(an_total, 0.975)
), by = .(decade, ssp, range)]

fwrite(table_s1_summary, "results/tables/Table_Lloyd2024_S1_Ranges.csv")

message("Table generation complete. Files saved in results/
