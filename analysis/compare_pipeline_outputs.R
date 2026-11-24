#!/usr/bin/env Rscript

################################################################################
# Compare saved Part 01, Part 03 and Part 04 outputs between two batch roots,
# including city ENSEMBLE outputs and pooled geography decompositions.
# Byte-identical files are reported separately; other files are joined by
# their declared data keys and checked at a relative tolerance.
# Requires: candidate and reference batch roots with matching SSP output paths.
# Produces: a concise comparison report on stdout; no files are written.
# How to run:
#   CANDIDATE_ROOT=results/europe/heatsplit_20261124_rebuild \
#   REFERENCE_ROOT=results/europe/heatsplit_20261013 SSP=3 \
#   Rscript analysis/compare_pipeline_outputs.R
################################################################################

suppressPackageStartupMessages(library(data.table))

candidate_root <- Sys.getenv(
  "CANDIDATE_ROOT", "results/europe/heatsplit_20261124_rebuild"
)
reference_root <- Sys.getenv(
  "REFERENCE_ROOT", "results/europe/heatsplit_20261013"
)
ssp <- Sys.getenv("SSP", "3")
tolerance <- as.numeric(Sys.getenv("RELATIVE_TOLERANCE", "1e-9"))
if (!grepl("^[1-3]$", ssp)) stop("SSP must be 1, 2, or 3.", call. = FALSE)
if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
  stop("RELATIVE_TOLERANCE must be a finite nonnegative number.", call. = FALSE)
}

candidate_root <- normalizePath(candidate_root, mustWork = TRUE)
reference_root <- normalizePath(reference_root, mustWork = TRUE)
candidate_ssp <- file.path(candidate_root, paste0("ssp", ssp))
reference_ssp <- file.path(reference_root, paste0("ssp", ssp))
if (!dir.exists(candidate_ssp) || !dir.exists(reference_ssp)) {
  stop("Both batch roots must contain the selected SSP directory.", call. = FALSE)
}

contracts <- list(
  "01_attribution_grouped.csv" = list(
    keys = c(
      "year", "branch", "agegroup", "range", "geo_id", "label", "ssp",
      "gcm", "days_in_year", "annualization_rule"
    ),
    values = "an"
  ),
  "03_master_table.csv" = list(
    keys = c("geo_id", "label", "ssp", "gcm", "branch", "year", "age", "range"),
    values = c("pop", "death", "an", "temp_deaths", "rest", "deaths_component")
  ),
  "04_le_li_levels.csv" = list(
    keys = c("geo_id", "label", "gcm", "branch", "year"),
    values = c("LE65", "LI65")
  ),
  "04_between_branch_decomposition.csv" = list(
    keys = c("geo_id", "label", "gcm", "period", "age", "cause"),
    values = c("le_contribution", "li_contribution")
  ),
  "04_within_branch_period_decomposition.csv" = list(
    keys = c(
      "geo_id", "label", "gcm", "branch", "period_from", "period_to",
      "age", "cause"
    ),
    values = c("le_contribution", "li_contribution")
  )
)

candidate_files <- list.files(
  candidate_ssp, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE
)
candidate_files <- candidate_files[basename(candidate_files) %in% names(contracts)]
if (!length(candidate_files)) {
  stop("No supported Part 01, Part 03, or Part 04 files found.", call. = FALSE)
}

compare_pair <- function(candidate, reference, contract) {
  candidate_size <- file.info(candidate)$size
  reference_size <- file.info(reference)$size
  if (identical(candidate_size, reference_size) &&
      identical(
        unname(tools::md5sum(candidate)),
        unname(tools::md5sum(reference))
      )) {
    return(list(status = "byte-identical", max_abs = 0, max_rel = 0,
                cause = ""))
  }

  candidate_data <- fread(candidate)
  reference_data <- fread(reference)
  if (!setequal(names(candidate_data), names(reference_data)) ||
      !all(c(contract$keys, contract$values) %in% names(candidate_data))) {
    return(list(status = "beyond-tolerance", max_abs = NA_real_,
                max_rel = NA_real_, cause = "schema mismatch"))
  }
  setcolorder(candidate_data, names(reference_data))
  if (anyDuplicated(candidate_data, by = contract$keys) ||
      anyDuplicated(reference_data, by = contract$keys)) {
    return(list(status = "beyond-tolerance", max_abs = NA_real_,
                max_rel = NA_real_, cause = "duplicate data keys"))
  }

  joined <- merge(
    reference_data, candidate_data, by = contract$keys, all = TRUE,
    sort = FALSE, suffixes = c("_reference", "_candidate")
  )
  if (nrow(joined) != nrow(reference_data) ||
      nrow(joined) != nrow(candidate_data) ||
      anyNA(joined[, ..contract$keys])) {
    return(list(status = "beyond-tolerance", max_abs = NA_real_,
                max_rel = NA_real_, cause = "row/key mismatch"))
  }

  max_abs <- 0
  max_rel <- 0
  changed <- character()
  non_numeric_values <- FALSE
  missing_value_mismatch <- FALSE
  for (column in contract$values) {
    reference_values <- joined[[paste0(column, "_reference")]]
    candidate_values <- joined[[paste0(column, "_candidate")]]
    if (!is.numeric(reference_values) || !is.numeric(candidate_values)) {
      non_numeric_values <- TRUE
      next
    }
    both_na <- is.na(reference_values) & is.na(candidate_values)
    if (any(xor(is.na(reference_values), is.na(candidate_values)))) {
      changed <- c(changed, column)
      missing_value_mismatch <- TRUE
      next
    }
    active <- which(!both_na)
    if (!length(active)) next
    finite <- is.finite(reference_values[active]) & is.finite(candidate_values[active])
    if (any(!finite)) {
      nonfinite <- active[!finite]
      if (!all(reference_values[nonfinite] == candidate_values[nonfinite])) {
        missing_value_mismatch <- TRUE
      }
    }
    active <- active[finite]
    if (!length(active)) next
    abs_diff <- abs(candidate_values[active] - reference_values[active])
    denom <- pmax(abs(candidate_values[active]), abs(reference_values[active]))
    rel_diff <- ifelse(denom == 0, 0, abs_diff / denom)
    max_abs <- max(max_abs, abs_diff, na.rm = TRUE)
    max_rel <- max(max_rel, rel_diff, na.rm = TRUE)
    if (any(abs_diff > 0)) changed <- c(changed, column)
  }

  if (non_numeric_values || missing_value_mismatch) {
    status <- "beyond-tolerance"
  } else if (max_rel <= tolerance) {
    status <- "within-tolerance"
  } else {
    status <- "beyond-tolerance"
  }
  list(
    status = status, max_abs = max_abs, max_rel = max_rel,
    cause = paste(unique(changed), collapse = ",")
  )
}

details <- rbindlist(lapply(candidate_files, function(candidate) {
  relative_path <- substring(candidate, nchar(candidate_ssp) + 2L)
  reference <- file.path(reference_ssp, relative_path)
  file_type <- basename(candidate)
  if (!file.exists(reference)) {
    return(data.table(
      file_type, relative_path, status = "missing-reference",
      max_abs = NA_real_, max_rel = NA_real_, cause = "reference file absent"
    ))
  }
  result <- compare_pair(candidate, reference, contracts[[file_type]])
  data.table(
    file_type, relative_path, status = result$status,
    max_abs = result$max_abs, max_rel = result$max_rel, cause = result$cause
  )
}))

summary <- details[, .(
  files = .N,
  byte_identical = sum(status == "byte-identical"),
  within_tolerance = sum(status == "within-tolerance"),
  beyond_tolerance = sum(status %chin% c("beyond-tolerance", "missing-reference")),
  max_abs_difference = if (all(is.na(max_abs))) NA_real_ else max(max_abs, na.rm = TRUE),
  max_relative_difference = if (all(is.na(max_rel))) NA_real_ else max(max_rel, na.rm = TRUE),
  difference_columns = paste(
    sort(unique(unlist(strsplit(cause[nzchar(cause)], ",", fixed = TRUE)))),
    collapse = ","
  )
), by = file_type]
print(summary)
cat(
  "Candidate root:", candidate_root, "\n",
  "Reference root:", reference_root, "\n",
  "SSP:", ssp, " Relative tolerance:", format(tolerance, scientific = TRUE), "\n",
  "Files:", nrow(details),
  " Byte-identical:", sum(details$status == "byte-identical"),
  " Within tolerance:", sum(details$status == "within-tolerance"),
  " Beyond/missing:", sum(details$status %chin% c("beyond-tolerance", "missing-reference")),
  "\n"
)

exceptions <- details[status %chin% c("beyond-tolerance", "missing-reference")]
if (nrow(exceptions)) {
  print(exceptions[order(-max_rel)][seq_len(min(nrow(exceptions), 20L))])
  quit(status = 1L)
}
