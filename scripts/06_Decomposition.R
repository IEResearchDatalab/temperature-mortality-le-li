################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 06_Decomposition.R
#   Lloyd-style decomposition of LE and LI by age and temperature range
#   For selected cities, GCMs, SSPs
#
################################################################################

if (length(ls()) == 0) source("05_LifeTables.R")

source("../R/period_lifetable.R")
source("../R/utils.R")

#---------------------------
# Load data
#---------------------------

analysis_all <- fread(file.path(path_out, "analysis_dataset_all_cities.csv"))

# Decomposition subset (first 10 cities for demo)
decomp_cities <- head(unique(analysis_all$city_code), 10)
decomp_years <- c(2030, 2050, 2090)

cat(sprintf("Decomposition for %d cities, %d years each...\n",
  length(decomp_cities), length(decomp_years)))

#---------------------------
# Decomposition function
#---------------------------

decompose_le_by_age <- function(yr_data) {

  ages <- yr_data$age
  mx_base <- yr_data$deaths_baseline / yr_data$population
  mx_base[is.nan(mx_base) | is.infinite(mx_base)] <- 1e-10
  mx_clim <- yr_data$deaths / yr_data$population
  mx_clim[is.nan(mx_clim) | is.infinite(mx_clim)] <- 1e-10

  compute_stats <- function(mx_vec)
    lifetable(mx = mx_vec, age = ages, ax = 0.5)

  lt_base <- compute_stats(mx_base)
  lt_clim <- compute_stats(mx_clim)

  e0_base <- lt_base[age == 20, ex]
  e0_clim <- lt_clim[age == 20, ex]
  total_delta <- e0_clim - e0_base

  sd_base <- sqrt(sum(lt_base$dx * (lt_base$age + lt_base$ax -
    sum(lt_base$dx * (lt_base$age + lt_base$ax) / lt_base$lx[1]))^2) /
    lt_base$lx[1])
  sd_clim <- sqrt(sum(lt_clim$dx * (lt_clim$age + lt_clim$ax -
    sum(lt_clim$dx * (lt_clim$age + lt_clim$ax) / lt_clim$lx[1]))^2) /
    lt_clim$lx[1])
  total_delta_sd <- sd_clim - sd_base

  # Stepwise replacement
  age_contrib <- rbindlist(lapply(seq_along(ages), function(i) {
    mx_hybrid <- mx_base
    mx_hybrid[i] <- mx_clim[i]
    lt_hybrid <- compute_stats(mx_hybrid)

    data.table(
      age = ages[i],
      delta_e0 = lt_hybrid[age == 20, ex] - e0_base,
      delta_sd = sqrt(sum(lt_hybrid$dx * (lt_hybrid$age + lt_hybrid$ax -
        sum(lt_hybrid$dx * (lt_hybrid$age + lt_hybrid$ax) /
          lt_hybrid$lx[1]))^2) / lt_hybrid$lx[1]) - sd_base)
  }))

  # Scale for exact additivity
  scale_e0 <- total_delta / sum(age_contrib$delta_e0)
  scale_sd <- total_delta_sd / sum(age_contrib$delta_sd)

  age_contrib[, `:=`(delta_e0 = delta_e0 * scale_e0,
    delta_sd = delta_sd * scale_sd)]

  list(age_contrib = age_contrib, total_delta_e0 = total_delta,
    total_delta_sd = total_delta_sd)
}

#---------------------------
# Run decomposition
#---------------------------

all_decomp <- list()

for (cc in decomp_cities) {
  cat(sprintf("  Decomposing: %s\n", cc))

  for (gc in unique(analysis_all$gcm)) {
    for (ss in unique(analysis_all$ssp)) {
      for (yr in decomp_years) {

        yr_data <- analysis_all[city_code == cc & gcm == gc &
          ssp == ss & year == yr]
        if (nrow(yr_data) < 10) next

        dc <- decompose_le_by_age(yr_data)

        # Attribute to temp ranges
        yr_data[, an_total := an_extreme_cold + an_moderate_cold +
          an_moderate_heat + an_extreme_heat]
        yr_data[an_total == 0, an_total := 1e-20]

        temp_ranges <- c("an_extreme_cold", "an_moderate_cold",
          "an_moderate_heat", "an_extreme_heat")
        temp_labels <- c("extreme_cold", "moderate_cold",
          "moderate_heat", "extreme_heat")

        for (k in seq_along(temp_ranges)) {
          tr <- temp_ranges[k]
          tl <- temp_labels[k]
          share <- pmax(yr_data[[tr]], 0) / yr_data$an_total

          all_decomp[[length(all_decomp) + 1]] <- data.table(
            city_code = cc, gcm = gc, ssp = ss, year = yr,
            age = dc$age_contrib$age,
            temp_range = tl,
            delta_e0 = dc$age_contrib$delta_e0 * share,
            delta_sd = dc$age_contrib$delta_sd * share)
        }
      }
    }
  }
}

decomp_all <- rbindlist(all_decomp)
fwrite(decomp_all, file.path(path_out, "decomposition_all_cities.csv"))

cat(sprintf("Decomposition complete: %s rows\n",
  format(nrow(decomp_all), big.mark = ",")))