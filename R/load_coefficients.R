################################################################################
#
# Excess mortality attributed to heat and cold: 
#   a health impact assessment study in 854 cities in Europe
#
# The Lancet Planetary Health, 2023
# https://doi.org/10.1016/S2542-5196(23)00023-2
#
# Coefficient loading and interpolation
#
################################################################################

#---------------------------
# Load city coefficients
#---------------------------

load_city_coefficients <- function(city_code,
  coefs_path = "data/coefs.csv") {

  coefs_all <- data.table::fread(coefs_path)
  coefs_city <- coefs_all[URAU_CODE == city_code]

  if (nrow(coefs_city) == 0)
    stop(sprintf("No coefficients found for city_code = '%s'",
      city_code))

  list(all = coefs_all, city = coefs_city)
}

#---------------------------
# Interpolate to single-year ages
#---------------------------

interpolate_coefs_to_single_age <- function(coefs_city, coefs_all,
  agelabs, age_midpoints, age_range = 20:100,
  city_code = NULL) {

  coef_cols <- names(coefs_city)[grepl("^b[0-9]+$", names(coefs_city))]
  if (length(coef_cols) == 0) stop("No coefficient columns b1..bk found.")

  coefs_city_ord <- coefs_city[match(agelabs, agegroup)]
  if (any(is.na(coefs_city_ord$agegroup)))
    stop("Some agelabs not found in coefs_city$agegroup.")

  coef_interp_dt <- data.table::data.table(
    agegroup = as.character(age_range))
  for (cc in coef_cols) {
    y <- coefs_city_ord[[cc]]
    coef_interp_dt[[cc]] <- stats::approx(x = age_midpoints, y = y,
      xout = age_range, rule = 2)$y
  }

  template <- coefs_city_ord[1]
  orig_cols <- names(coefs_all)
  coefs_single_age <- template[rep(1, length(age_range))]

  coefs_single_age[, agegroup := as.character(age_range)]
  for (cc in coef_cols)
    coefs_single_age[[cc]] <- coef_interp_dt[[cc]]

  if (!is.null(city_code) && "URAU_CODE" %in% names(coefs_single_age))
    coefs_single_age[, URAU_CODE := city_code]

  missing_cols <- setdiff(orig_cols, names(coefs_single_age))
  if (length(missing_cols) > 0)
    stop(sprintf("Interpolated coefs missing columns: %s",
      paste(missing_cols, collapse = ", ")))

  coefs_single_age[, ..orig_cols]
}