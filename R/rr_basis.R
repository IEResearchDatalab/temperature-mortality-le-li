################################################################################
#
# Excess mortality attributed to heat and cold:
#   a health impact assessment study in 854 cities in Europe
#
# The Lancet Planetary Health, 2023
# https://doi.org/10.1016/S2542-5196(23)00023-2
#
# RR basis functions
#
################################################################################

library(dlnm)
library(splines)

#---------------------------
# Build basis parameters
#---------------------------

build_basis_params <- function(hist_temps, varfun = "bs", vardegree = 2,
  varper = c(10, 75, 90)) {

  varknots <- quantile(hist_temps, varper / 100, na.rm = TRUE)
  varbound <- range(hist_temps, na.rm = TRUE)

  argvar <- list(fun = varfun, degree = vardegree,
    knots = varknots, Bound = varbound)

  list(argvar = argvar, varknots = varknots, varbound = varbound)
}

#---------------------------
# Compute RR curves for each age group
#---------------------------

compute_rr_curves <- function(coefs_city, agelabs, argvar, varbound,
  temp_step = 0.5) {

  temp_seq <- seq(varbound[1], varbound[2], by = temp_step)
  n_temp <- length(temp_seq)

  basis <- do.call(onebasis, c(list(x = temp_seq), argvar))

  rr_matrix <- matrix(NA, nrow = n_temp, ncol = length(agelabs))
  mmt_vec <- numeric(length(agelabs))

  for (i in seq_along(agelabs)) {
    age <- agelabs[i]
    coef_row <- coefs_city[agegroup == age]
    coefs <- as.numeric(coef_row[, .(b1, b2, b3, b4, b5)])

    log_rr <- basis %*% coefs

    ind <- temp_seq >= quantile(temp_seq, 0.25) &
           temp_seq <= quantile(temp_seq, 0.99)
    mmt <- temp_seq[ind][which.min(log_rr[ind])]
    mmt_vec[i] <- mmt

    cenvec <- do.call(onebasis, c(list(x = mmt), argvar))
    log_rr_centered <- log_rr - drop(cenvec %*% coefs)

    rr <- exp(log_rr_centered)
    rr_matrix[, i] <- as.vector(rr)
  }

  list(temp_seq = temp_seq, rr_matrix = rr_matrix, mmt_vec = mmt_vec,
    basis = basis)
}

#---------------------------
# Compute daily AN (Masselot-style)
#---------------------------

compute_daily_an <- function(temps, coefs, mmt, annual_deaths,
  argvar, n_days_year) {

  bvar <- do.call(onebasis, c(list(x = temps), argvar))
  cenvec <- do.call(onebasis, c(list(x = rep(mmt, length(temps))),
    argvar))
  bvarcen <- scale(bvar, center = cenvec, scale = FALSE)

  af_day <- 1 - exp(-bvarcen %*% coefs)
  af_day[af_day < 0] <- 0  # Gasparrini & Leone (2014): range approach
  daily_deaths <- annual_deaths / n_days_year
  an_day <- af_day * daily_deaths

  list(af = af_day, an = an_day)
}

#---------------------------
# Classify temperatures into ranges
#---------------------------

classify_temp_range <- function(temp, p025, p975, mmt) {
  fifelse(temp < p025, "extreme_cold",
    fifelse(temp < mmt, "moderate_cold",
      fifelse(temp <= p975, "moderate_heat", "extreme_heat")))
}

#---------------------------
# Interpolate RR to single-year ages
#---------------------------

interpolate_rr_to_single_age <- function(rr_matrix, mmt_vec,
  age_midpoints, age_range = 20:100) {

  n_temp <- nrow(rr_matrix)
  rr_single_age <- matrix(NA, nrow = n_temp, ncol = length(age_range))
  colnames(rr_single_age) <- age_range

  for (t_idx in seq_len(n_temp)) {
    rr_at_temp <- rr_matrix[t_idx, ]
    rr_single_age[t_idx, ] <- approx(x = age_midpoints, y = rr_at_temp,
      xout = age_range, rule = 2)$y
  }

  mmt_single_age <- approx(x = age_midpoints, y = mmt_vec,
    xout = age_range, rule = 2)$y

  list(rr_single_age = rr_single_age, mmt_single_age = mmt_single_age)
}

#---------------------------
# Compute average RR by age
#---------------------------

compute_avg_rr_by_age <- function(temps, temp_seq, rr_single_age,
  mmt_single_age, age_range, component = "total",
  doys = NULL, sw_matrix = NULL) {

  temps <- temps[!is.na(temps)]
  if (length(temps) == 0) return(rep(NA_real_, length(age_range)))

  temp_indices <- vapply(temps, function(t) which.min(abs(temp_seq - t)),
    integer(1))

  rr_vals <- rr_single_age[temp_indices, , drop = FALSE]

  if (component != "total") {
    for (j in seq_along(age_range)) {
      mmt <- mmt_single_age[j]
      if (component == "heat") {
        rr_vals[temps <= mmt, j] <- 1
      } else if (component == "cold") {
        rr_vals[temps > mmt, j] <- 1
      }
    }
  }

  if (!is.null(sw_matrix) && !is.null(doys)) {
    avg_rr <- numeric(length(age_range))
    for (j in seq_along(age_range)) {
      w <- sw_matrix[as.character(age_range[j]), doys]
      avg_rr[j] <- weighted.mean(rr_vals[, j], w)
    }
  } else {
    avg_rr <- colMeans(rr_vals)
  }

  names(avg_rr) <- as.character(age_range)
  avg_rr
}