################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 05_LifeTables.R
#   Compute period life tables and lifespan inequality
#   For all cities, GCMs, SSPs
#
################################################################################

if (length(ls()) == 0) source("04_AnalysisDataset.R")

source("../R/period_lifetable.R")
source("../R/utils.R")

#---------------------------
# Load analysis dataset
#---------------------------

analysis_all <- fread(file.path(path_out, "analysis_dataset_all_cities.csv"))

#---------------------------
# Compute LE and LI per city-GCM-SSP-year
#---------------------------

compute_lifetable_stats <- function(mx_vec, age_vec, ax = 0.5) {
  keep <- age_vec >= 65
  mx_vec <- mx_vec[keep]
  age_vec <- age_vec[keep]
  lt <- lifetable(mx = mx_vec, age = age_vec, ax = ax)
  lt <- lt[order(age)]

  e65 <- lt[age == 65, ex]

  lx_start <- lt[1, lx]
  mean_age <- sum(lt$dx * (lt$age + lt$ax)) / lx_start
  sd <- sqrt(sum(lt$dx * (lt$age + lt$ax - mean_age)^2) / lx_start)

  list(e65 = e65, sd = sd, lt = lt)
}

all_stats <- list()
combos <- unique(analysis_all[, .(city_code, gcm, ssp, year)])
n_combos <- nrow(combos)

cat(sprintf("Life tables for %d city-GCM-SSP-year combos...\n", n_combos))

for (i in seq_len(n_combos)) {
  if (i %% 1000 == 0)
    cat(sprintf("  %d / %d\n", i, n_combos))

  cc <- combos$city_code[i]
  gc <- combos$gcm[i]
  ss <- combos$ssp[i]
  yr <- combos$year[i]

  yr_data <- analysis_all[city_code == cc & gcm == gc &
    ssp == ss & year == yr]

  if (nrow(yr_data) == 0) next

  mx_base <- yr_data$deaths_baseline / yr_data$population
  mx_base[is.nan(mx_base) | is.infinite(mx_base)] <- 1e-10

  mx_clim <- yr_data$deaths / yr_data$population
  mx_clim[is.nan(mx_clim) | is.infinite(mx_clim)] <- 1e-10

  s_base <- compute_lifetable_stats(mx_base, yr_data$age)
  s_clim <- compute_lifetable_stats(mx_clim, yr_data$age)

  all_stats[[length(all_stats) + 1]] <- data.table(
    city_code = cc, gcm = gc, ssp = ss, year = yr,
    e65_base = s_base$e65, e65_clim = s_clim$e65,
    sd_base = s_base$sd, sd_clim = s_clim$sd)
}

stats_all <- rbindlist(all_stats)
fwrite(stats_all, file.path(path_out, "lifespan_inequality_all_cities.csv"))

cat(sprintf("Life tables complete: %s rows\n",
  format(nrow(stats_all), big.mark = ",")))