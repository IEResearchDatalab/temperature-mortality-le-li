################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 07_LifeTables.R
#   Build period life tables from age 65 onward.
#   Calculate remaining life expectancy at 65 (e65) and lifespan
#   inequality (LI) above age 65.
#
# LI definition: SD of attained age at death above 65.
#   SD = sqrt( Σ dx * (x + ax - A)² / lx_65 )
#   where A = Σ dx * (x + ax) / lx_65
#
# All life-table calculations are explicit in this script (not imported
# from a custom package helper).
#
# Input: analysis_dataset_all_cities.csv
#   Contains: city_code, gcm, ssp, year, age, population, mx, deaths,
#             extreme_cold, moderate_cold, moderate_heat, extreme_heat,
#             residual_deaths
#
# Output: lifespan_inequality_all_cities.csv
#   Contains: city_code, gcm, ssp, year, e65, sd (LI above 65)
#
################################################################################

if (!exists("path_out")) source("00_Packages_Parameters.R")

#---------------------------
# Load analysis dataset
#---------------------------

analysis_all <- fread(file.path(path_out, "analysis_dataset_all_cities.csv"))

#---------------------------
# Life-table function (explicit, ages 65+ only)
#---------------------------

build_lifetable_65plus <- function(mx_vec, age_vec, ax = 0.5) {
  keep <- age_vec >= 65
  mx <- mx_vec[keep]
  age <- age_vec[keep]
  ord <- order(age)

  mx <- mx[ord]
  age <- age[ord]
  n <- c(diff(age), Inf)

  # qx from mx (Chiang method)
  qx <- (n * mx) / (1 + n * (1 - ax) * mx)
  qx[length(qx)] <- 1

  # lx, dx
  lx <- numeric(length(qx))
  lx[1] <- 100000
  for (i in 2:length(qx)) {
    lx[i] <- lx[i - 1] * (1 - qx[i - 1])
  }
  dx <- lx * qx

  # Lx, Tx, ex
  Lx <- numeric(length(qx))
  for (i in 1:(length(qx) - 1)) {
    Lx[i] <- n[i] * (lx[i] - (1 - ax) * dx[i])
  }
  Lx[length(qx)] <- if (mx[length(qx)] > 0) lx[length(qx)] / mx[length(qx)] else 0

  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  ex[is.nan(ex)] <- 0

  # Remaining life expectancy at 65
  e65 <- ex[age == 65]

  # Lifespan inequality: SD of attained age at death above 65
  lx_start <- lx[1]
  mean_age <- sum(dx * (age + ax)) / lx_start
  sd <- sqrt(sum(dx * (age + ax - mean_age)^2) / lx_start)

  list(e65 = e65, sd = sd, lt = data.table(age = age, n = n, mx = mx,
    qx = qx, lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex, ax = ax))
}

#---------------------------
# Compute per city-GCM-SSP-year
#---------------------------

combos <- unique(analysis_all[, .(city_code, gcm, ssp, year)])
n_combos <- nrow(combos)

cat(sprintf("Life tables for %d city-GCM-SSP-year combos...\n", n_combos))

all_stats <- list()

for (i in seq_len(n_combos)) {
  if (i %% 1000 == 0)
    cat(sprintf("  %d / %d\n", i, n_combos))

  cc <- combos$city_code[i]
  gc <- combos$gcm[i]
  ss <- combos$ssp[i]
  yr <- combos$year[i]

  yr_data <- analysis_all[city_code == cc & gcm == gc & ssp == ss & year == yr]
  if (nrow(yr_data) == 0) next

  setorder(yr_data, age)
  ages <- yr_data$age

  # All-cause mortality (includes temperature component)
  mx_all <- yr_data$deaths / yr_data$population
  mx_all[is.nan(mx_all) | is.infinite(mx_all)] <- 1e-10

  result <- build_lifetable_65plus(mx_all, ages)

  all_stats[[length(all_stats) + 1]] <- data.table(
    city_code = cc, gcm = gc, ssp = ss, year = yr,
    e65 = result$e65,
    sd = result$sd)
}

stats_all <- rbindlist(all_stats)

#---------------------------
# Save
#---------------------------

fwrite(stats_all, file.path(path_out, "lifespan_inequality_all_cities.csv"))
cat(sprintf("Life tables complete: %s rows\n",
  format(nrow(stats_all), big.mark = ",")))