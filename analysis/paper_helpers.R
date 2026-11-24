################################################################################
# Shared age-band and period labels for Figures 2 and 3.
# Requires: data.table when age_band() is called.
# Produces: age_band() and period_start() functions in the calling session.
# How to run: Rscript -e 'source("analysis/paper_helpers.R")'
# Notes: the functions label existing data; they do not alter the decomposition.
################################################################################

age_band <- function(age) {
  data.table::fifelse(age >= 100L, "100+", sprintf("%d-%d", (age %/% 5L) * 5L, (age %/% 5L) * 5L + 4L))
}

period_start <- function(x) as.integer(substr(x, 1L, 4L))
