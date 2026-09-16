#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_path <- if (length(file_arg)) normalizePath(file_arg) else normalizePath("tests/run_tests.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
tmp_out <- file.path(tempdir(), paste0("phase0-test-", format(Sys.time(), "%Y%m%d-%H%M%S")))
dir.create(tmp_out, recursive = TRUE, showWarnings = FALSE)
setwd(repo_root)

cmd <- "Rscript"
args <- c("validation/run_validation.R")
status <- system2(cmd, args, env = c(paste0("VALIDATION_OUT_ROOT=", tmp_out)), stdout = TRUE, stderr = TRUE)

required <- c(
  file.path(tmp_out, "tables", "gate_dashboard.csv"),
  file.path(tmp_out, "tables", "clamp_comparison.csv"),
  file.path(tmp_out, "tables", "masselot_reproduction.csv"),
  file.path(tmp_out, "figures", "gate_dashboard.png"),
  file.path(tmp_out, "session_info.txt"),
  file.path(tmp_out, "output_manifest.csv")
)

missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing test outputs: ", paste(missing, collapse = ", "))

cat("run_tests.R smoke test passed; output dir:", tmp_out, "\n")
