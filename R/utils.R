################################################################################
#
# Utility functions
#
################################################################################

`%+%` <- function(a, b) paste0(a, b)

Rep <- function(x, n) paste(rep(x, n), collapse = "")

get_ax <- function(age) {
  ifelse(age == 0, 0.1, 0.5)
}

mx_to_qx <- function(mx, ax = 0.5) {
  mx / (1 + (1 - ax) * mx)
}

cat_header <- function(title, width = 72) {
  cat("\n" %+% Rep("=", width) %+% "\n")
  cat(title %+% "\n")
  cat(Rep("=", width) %+% "\n\n")
}

cat_step <- function(step_num, description) {
  cat(sprintf("\nStep %d: %s...\n", step_num, description))
}