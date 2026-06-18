################################################################################
#
# Excess mortality attributed to heat and cold: 
#   a health impact assessment study in 854 cities in Europe
#
# The Lancet Planetary Health, 2023
# https://doi.org/10.1016/S2542-5196(23)00023-2
#
# Impact computation and aggregation functions
#
################################################################################

#--------------------------------
# Compute impact measures and aggregate results
#--------------------------------

impact <- function(res, measures = c("an", "af", "rate", "cuman"), 
  perlen = 1, ensonly = TRUE, warming_win = NULL) {

  # Compute attributable fraction
  if ("af" %in% measures)
    res[, af := an / death]

  # Compute death rates
  if ("rate" %in% measures)
    res[, rate := an / pop]

  # Compute cumulative AN
  if ("cuman" %in% measures)
    res[order(year), cuman := cumsum(an), 
      by = c("gcm", "range", "sc", "agegroup", "res")]

  # Fill NAs
  setnafill(res, fill = 0, cols = measures)

  # Results by period
  res[, period := floor(year / perlen) * perlen]

  estres <- res[res == "est", lapply(.SD, mean, na.rm = TRUE), 
    by = c("period", "sc", "range", "agegroup"), .SDcols = measures]

  cires <- res[res != "est", 
    as.list(unlist(lapply(.SD, fquantile, c(.025, .975), na.rm = TRUE))),
    by = c("period", "sc", "range", "agegroup"), .SDcols = measures]

  if (!ensonly) {
    estresgcm <- res[res == "est", lapply(.SD, mean, na.rm = TRUE), 
      by = c("period", "sc", "range", "agegroup", "gcm"), 
      .SDcols = measures]
    estres <- rbind(estres[, gcm := "ens"], estresgcm)
    ciresgcm <- res[res != "est", 
      as.list(unlist(lapply(.SD, fquantile, c(.025, .975), na.rm = TRUE))),
      by = c("period", "sc", "range", "agegroup", "gcm"), 
      .SDcols = measures]
    cires <- rbind(cires[, gcm := "ens"], ciresgcm)
    rm(estresgcm, ciresgcm)
  }

  names(cires) <- gsub("\\.97.*\\%", "_high", names(cires))
  names(cires) <- gsub("\\.2.*\\%", "_low", names(cires))

  periodres <- merge(estres, cires)
  rm(estres, cires); gc()

  out <- list(period = periodres)

  # Results by warming level
  if (!is.null(warming_win)) {
    res <- merge(res, warming_win, by = c("gcm", "year"), 
      allow.cartesian = TRUE)

    estres <- res[res == "est", lapply(.SD, mean, na.rm = TRUE), 
      by = c("level", "sc", "range", "agegroup"), .SDcols = measures]

    cires <- res[res != "est", 
      as.list(unlist(lapply(.SD, fquantile, c(.025, .975), na.rm = TRUE))),
      by = c("level", "sc", "range", "agegroup"), .SDcols = measures]

    if (!ensonly) {
      estresgcm <- res[res == "est", lapply(.SD, mean, na.rm = TRUE), 
        by = c("level", "sc", "range", "agegroup", "gcm"), 
        .SDcols = measures]
      estres <- rbind(estres[, gcm := "ens"], estresgcm)
      ciresgcm <- res[res != "est", 
        as.list(unlist(lapply(.SD, fquantile, c(.025, .975), na.rm = TRUE))),
        by = c("level", "sc", "range", "agegroup", "gcm"), 
        .SDcols = measures]
      cires <- rbind(cires[, gcm := "ens"], ciresgcm)
      rm(estresgcm, ciresgcm)
    }

    names(cires) <- gsub("\\.97.*\\%", "_high", names(cires))
    names(cires) <- gsub("\\.2.*\\%", "_low", names(cires))

    levelres <- merge(estres, cires)
    rm(estres, cires); gc()

    out$level <- levelres
  }

  rm(res)
  out
}

#--------------------------------
# Aggregate results
#--------------------------------

impact_aggregate <- function(res, agg, vars = "an", by = key(res)) {
  for (vi in names(agg)) {
    restot <- res[, lapply(.SD, sum), 
      by = setdiff(by, vi), .SDcols = vars]
    res <- rbind(res, restot[, (vi) := agg[vi]])
  }
  res
}

#--------------------------------
# Compute impact measures
#--------------------------------

impact_measures <- function(res, vars = "an", 
  measures = c("af", "rate", "cuman"), by = key(res)) {

  if ("af" %in% measures)
    res[, sprintf("af_%s", vars) := lapply(.SD, "/", death), 
      .SDcols = vars]

  if ("rate" %in% measures)
    res[, sprintf("rate_%s", vars) := lapply(.SD, "/", pop), 
      .SDcols = vars]

  if ("cuman" %in% measures)
    res[, sprintf("cuman_%s", vars) := lapply(.SD, cumsum), 
      .SDcols = vars, by = setdiff(by, "year")]

  newvars <- outer(measures, vars, paste, sep = "_") |> c()
  setnafill(res, fill = 0, cols = newvars)

  res
}

#--------------------------------
# Summarise with confidence intervals
#--------------------------------

impact_summarise <- function(res, vars = "an", by = key(res), 
  prob = .95) {

  estres <- res[res == "est", lapply(.SD, mean, na.rm = TRUE), 
    by = by, .SDcols = vars]

  alpha <- 1 - prob
  lims <- c(alpha / 2, 1 - (alpha / 2))
  cires <- res[res != "est", 
    as.list(unlist(lapply(.SD, fquantile, lims, na.rm = TRUE))),
    by = by, .SDcols = vars]

  setnames(estres, vars, sprintf("%s_est", vars))
  names(cires) <- gsub(sprintf("\\.%s\\%%", lims[1] * 100), "_low", 
    names(cires))
  names(cires) <- gsub(sprintf("\\.%s\\%%", lims[2] * 100), "_high", 
    names(cires))

  merge(estres, cires)
}