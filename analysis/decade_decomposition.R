#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DemoDecomp)
})

decade_cause_levels <- c(
  "ExtrCold", "ModCold", "ModHeat", "NormExtrHeat", "ExtrHeat", "rest"
)
decade_age_levels <- 65:100
decade_nx <- c(rep(1, 100 - 65), Inf)

decade_life_expectancy <- function(mx) {
  dim(mx) <- c(length(decade_age_levels), length(mx) / length(decade_age_levels))
  total_mx <- rowSums(mx)
  px <- exp(-total_mx * decade_nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(total_mx == 0, lx * decade_nx, dx / total_mx)
  Tx <- rev(cumsum(rev(Lx)))
  Tx[1] / lx[1]
}

decade_lifespan_inequality <- function(mx) {
  dim(mx) <- c(length(decade_age_levels), length(mx) / length(decade_age_levels))
  total_mx <- rowSums(mx)
  px <- exp(-total_mx * decade_nx)
  lx <- head(cumprod(c(1, px)), -1)
  dx <- c(-diff(lx), tail(lx, 1))
  Lx <- ifelse(total_mx == 0, lx * decade_nx, dx / total_mx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  x_conditional <- decade_age_levels - 65
  sqrt(sum(dx * (x_conditional + 0.5 - ex[1])^2))
}

decade_grid <- function(values) {
  data.table(
    age = rep(decade_age_levels, length(decade_cause_levels)),
    cause = rep(decade_cause_levels, each = length(decade_age_levels)),
    contribution = as.numeric(values)
  )
}

decade_cause_vector <- function(schedule) {
  ordered <- schedule[order(match(cause, decade_cause_levels), age)]
  if (!identical(as.integer(unique(ordered$age)), decade_age_levels) ||
      !identical(unique(ordered$cause), decade_cause_levels)) {
    stop("Decade mortality schedule has an incomplete age/cause grid.", call. = FALSE)
  }
  ordered$mx_cause
}

compute_decade_decomposition <- function(
  batch_root,
  geo_level,
  geo_id,
  ssp_set = 1:3,
  output_dir,
  n_horiuchi = 400L,
  branch = "with_cc"
) {
  if (!length(ssp_set) || anyNA(ssp_set) || any(!ssp_set %in% 1:3)) {
    stop("ssp_set must be a subset of 1, 2, 3.", call. = FALSE)
  }
  if (length(n_horiuchi) != 1L || is.na(n_horiuchi) || n_horiuchi < 1L) {
    stop("n_horiuchi must be a positive integer.", call. = FALSE)
  }

  master_file <- function(ssp) file.path(
    batch_root, "geographies", paste0("ssp", ssp), geo_level, geo_id,
    "03_master_table.csv"
  )
  missing <- vapply(ssp_set, function(ssp) !file.exists(master_file(ssp)), logical(1))
  if (any(missing)) {
    stop("Missing annual pooled master table(s): ",
         paste(master_file(ssp_set[missing]), collapse = ", "), call. = FALSE)
  }

  branch_name <- branch
  all_results <- list()
  all_checks <- list()
  for (ssp in ssp_set) {
    master <- fread(master_file(ssp))
    master <- master[gcm == "ENSEMBLE" & branch == branch_name &
                       age %in% decade_age_levels]
    if (!nrow(master) || any(!is.finite(master$pop)) ||
        any(!is.finite(master$deaths_component))) {
      stop("Annual master table has no valid with-CC ensemble rows for SSP",
           ssp, ".", call. = FALSE)
    }

    annual <- master[, .(
      deaths = sum(deaths_component),
      person_years = sum(pop)
    ), by = .(year, age, cause = range)]
    annual[, decade := (year %/% 10L) * 10L]
    decade_rates <- annual[, .(
      mx_cause = sum(deaths) / sum(person_years)
    ), by = .(decade, age, cause)]

    expected <- length(seq(2020L, 2090L, by = 10L)) *
      length(decade_age_levels) * length(decade_cause_levels)
    if (nrow(decade_rates) != expected ||
        anyDuplicated(decade_rates, by = c("decade", "age", "cause")) ||
        any(!decade_rates$cause %in% decade_cause_levels)) {
      stop("Annual data do not form a complete decade/age/cause grid for SSP",
           ssp, ".", call. = FALSE)
    }

    decade_rates[, mx_total := sum(mx_cause), by = .(decade, age)]
    decade_levels <- decade_rates[cause == decade_cause_levels[1], .(
      LE65 = decade_life_expectancy(mx_total[order(age)]),
      LI65 = decade_lifespan_inequality(mx_total[order(age)])
    ), by = decade]
    transitions <- list()
    checks <- list()
    decades <- sort(unique(decade_rates$decade))

    for (i in seq_len(length(decades) - 1L)) {
      earlier <- decades[i]
      later <- decades[i + 1L]
      mx_earlier <- decade_cause_vector(decade_rates[decade == earlier])
      mx_later <- decade_cause_vector(decade_rates[decade == later])
      le_contrib <- as.vector(horiuchi(
        func = decade_life_expectancy,
        pars1 = mx_earlier,
        pars2 = mx_later,
        N = as.integer(n_horiuchi)
      ))
      li_contrib <- as.vector(horiuchi(
        func = decade_lifespan_inequality,
        pars1 = mx_earlier,
        pars2 = mx_later,
        N = as.integer(n_horiuchi)
      ))
      contribution <- decade_grid(le_contrib)
      le_change <- decade_levels[decade == later, LE65] -
        decade_levels[decade == earlier, LE65]
      li_change <- decade_levels[decade == later, LI65] -
        decade_levels[decade == earlier, LI65]
      closure_error <- sum(le_contrib) - le_change
      li_closure_error <- sum(li_contrib) - li_change
      checks[[length(checks) + 1L]] <- data.table(
        ssp = ssp,
        decade_from = earlier,
        decade_to = later,
        le_change = le_change,
        li_change = li_change,
        le_closure_error = closure_error,
        li_closure_error = li_closure_error,
        status = if (abs(closure_error) < 1e-6 &&
                     abs(li_closure_error) < 1e-6) "PASS" else "FAIL"
      )
      contribution[, `:=`(
        ssp = ssp,
        decade_from = earlier,
        decade_to = later,
        le_contribution_months = contribution * 12,
        li_contribution_days = li_contrib * 365
      )]
      transitions[[length(transitions) + 1L]] <- contribution
    }

    all_results[[as.character(ssp)]] <- rbindlist(transitions)
    all_checks[[as.character(ssp)]] <- rbindlist(checks)
  }

  result <- rbindlist(all_results)
  checks <- rbindlist(all_checks)
  if (any(checks$status != "PASS")) {
    failed <- checks[status != "PASS"]
    stop(
      "Decade Horiuchi LE closure failed: ",
      paste(sprintf(
        "SSP%s %d-%d error %.6g",
        failed$ssp, failed$decade_from, failed$decade_to,
          pmax(abs(failed$le_closure_error), abs(failed$li_closure_error))
      ), collapse = "; "),
      call. = FALSE
    )
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- if (identical(branch, "with_cc")) "" else paste0("_", branch)
  stem <- file.path(output_dir, paste0("decade_decomposition_", geo_level, "_", geo_id, suffix))
  fwrite(result, paste0(stem, ".csv"))
  fwrite(checks, paste0(stem, "_checks.csv"))
  list(data = result, checks = checks)
}
