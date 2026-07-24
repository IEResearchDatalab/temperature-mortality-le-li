pclm_expected_length <- function(x, nlast) {
  sum(diff(x)) + nlast
}

pclm_disaggregate_nonnegative <- function(x, y, nlast, zero_guard = 1e-8) {
  stopifnot(length(x) == length(y))

  target_length <- pclm_expected_length(x, nlast)
  if (all(is.na(y))) {
    stop("Input contains only NA values.")
  }
  if (any(y < 0, na.rm = TRUE)) {
    stop("Nonnegative disaggregation received negative values.")
  }
  if (sum(y, na.rm = TRUE) == 0) {
    return(rep(0, target_length))
  }

  y_work <- as.numeric(y)
  y_work[y_work == 0] <- zero_guard
  fit <- suppressWarnings(ungroup::pclm(x = x, y = y_work, nlast = nlast)$fitted)
  fit <- as.numeric(fit)

  fit <- fit * (sum(y, na.rm = TRUE) / sum(fit, na.rm = TRUE))
  fit
}

pclm_disaggregate_signed <- function(x, y, nlast, zero_guard = 1e-8) {
  stopifnot(length(x) == length(y))

  y <- as.numeric(y)
  pos <- pmax(y, 0)
  neg <- pmax(-y, 0)

  fit_pos <- pclm_disaggregate_nonnegative(x = x, y = pos, nlast = nlast, zero_guard = zero_guard)
  fit_neg <- pclm_disaggregate_nonnegative(x = x, y = neg, nlast = nlast, zero_guard = zero_guard)
  fit <- fit_pos - fit_neg

  # Preserve signed totals exactly up to floating-point tolerance.
  correction <- sum(y, na.rm = TRUE) - sum(fit, na.rm = TRUE)
  if (abs(correction) > 0) {
    fit[length(fit)] <- fit[length(fit)] + correction
  }
  fit
}