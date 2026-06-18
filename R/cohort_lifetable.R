################################################################################
#
# Cohort life table with climate-adjusted mortality
#
################################################################################

build_cohort_lifetable <- function(mort_proj_dt, mult_dt, ssp_val,
  adapt_lab, rcp_labels, cohort_start_age, cohort_years,
  radix = 100000) {

  rcp_lab <- rcp_labels[ssp_val]

  cohort_age <- cohort_start_age:(cohort_start_age +
    length(cohort_years) - 1)
  cohort_years_vec <- cohort_years

  lt <- data.table::data.table(
    age = cohort_age,
    year = cohort_years_vec[1:length(cohort_age)],
    rcp = rcp_lab,
    adaptation = adapt_lab)

  lt <- merge(lt,
    mort_proj_dt[, .(year, age, qx_base = qx, mx_base = mx, ax)],
    by = c("year", "age"), all.x = TRUE)

  if (any(is.na(lt$qx_base))) {
    lt[is.na(qx_base),
      qx_base := mort_proj_dt[year == max(year) & age == .BY$age, qx],
      by = age]
    lt[is.na(mx_base),
      mx_base := mort_proj_dt[year == max(year) & age == .BY$age, mx],
      by = age]
    lt[is.na(ax), ax := 0.5]
  }

  lt <- merge(lt,
    mult_dt[ssp == ssp_val & adaptation == adapt_lab,
      .(year, age, multiplier)],
    by = c("year", "age"), all.x = TRUE)
  lt[is.na(multiplier), multiplier := 1]

  lt[, mx_clim := mx_base * multiplier]
  lt[, qx_clim := mx_to_qx(mx_clim, ax)]

  lt[qx_base > 1, qx_base := 1]
  lt[qx_clim > 1, qx_clim := 1]

  lt <- lt[order(age)]

  lt[, lx_base := radix]
  for (i in 2:nrow(lt))
    lt$lx_base[i] <- lt$lx_base[i - 1] * (1 - lt$qx_base[i - 1])
  lt[, dx_base := lx_base * qx_base]

  lt[, lx_clim := radix]
  for (i in 2:nrow(lt))
    lt$lx_clim[i] <- lt$lx_clim[i - 1] * (1 - lt$qx_clim[i - 1])
  lt[, dx_clim := lx_clim * qx_clim]

  lt
}