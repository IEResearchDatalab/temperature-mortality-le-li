################################################################################
#
# Simulation and confidence intervals
#
################################################################################

library(MASS)

#---------------------------
# Simulate AN with CIs (Masselot-style)
#---------------------------

simulate_an <- function(temps, coefs, vcov_mat, mmt, annual_deaths,
  argvar, nsim = 1000, n_days_year = 365) {

  bvar <- do.call(onebasis, c(list(x = temps), argvar))
  cenvec <- do.call(onebasis,
    c(list(x = rep(mmt, length(temps))), argvar))
  bvarcen <- scale(bvar, center = cenvec, scale = FALSE)

  daily_deaths <- annual_deaths / n_days_year

  # Point estimate
  af_est <- 1 - exp(-bvarcen %*% coefs)
  an_est <- sum(af_est * daily_deaths)

  # Simulated coefficients
  set.seed(13042025)
  coef_sim <- mvrnorm(nsim, coefs, vcov_mat)

  # Simulated ANs
  an_sim <- apply(coef_sim, 1, function(cf) {
    af <- 1 - exp(-bvarcen %*% cf)
    af[af < 0] <- 0
    sum(af * daily_deaths)
  })

  # CIs
  ci <- quantile(an_sim, c(.025, .975), na.rm = TRUE)

  list(est = an_est, low = ci[1], hi = ci[2], sim = an_sim)
}

#---------------------------
# Simulate by temperature range
#---------------------------

simulate_an_by_range <- function(temps, coefs, vcov_mat, mmt, p025, p975,
  annual_deaths, argvar, nsim = 1000, n_days_year = 365) {

  bvar <- do.call(onebasis, c(list(x = temps), argvar))
  cenvec <- do.call(onebasis,
    c(list(x = rep(mmt, length(temps))), argvar))
  bvarcen <- scale(bvar, center = cenvec, scale = FALSE)

  daily_deaths <- annual_deaths / n_days_year

  # Temperature range indicators
  is_extreme_cold <- temps < p025
  is_moderate_cold <- temps >= p025 & temps < mmt
  is_moderate_heat <- temps >= mmt & temps <= p975
  is_extreme_heat <- temps > p975

  ranges <- list(
    extreme_cold = is_extreme_cold,
    moderate_cold = is_moderate_cold,
    moderate_heat = is_moderate_heat,
    extreme_heat = is_extreme_heat)

  set.seed(13042025)
  coef_sim <- mvrnorm(nsim, coefs, vcov_mat)

  res <- lapply(names(ranges), function(rname) {
    idx <- ranges[[rname]]
    if (sum(idx) == 0)
      return(data.table(range = rname, est = 0, low = 0, hi = 0))

    af_est <- 1 - exp(-bvarcen[idx, , drop = FALSE] %*% coefs)
    an_est <- sum(af_est * daily_deaths[idx])

    an_sim <- apply(coef_sim, 1, function(cf) {
      af <- 1 - exp(-bvarcen[idx, , drop = FALSE] %*% cf)
      sum(af * daily_deaths[idx])
    })

    ci <- quantile(an_sim, c(.025, .975), na.rm = TRUE)
    data.table(range = rname, est = an_est, low = ci[1], hi = ci[2])
  })

  rbindlist(res)
}