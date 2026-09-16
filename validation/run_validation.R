#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_path <- if (length(file_arg)) normalizePath(file_arg) else normalizePath("validation/run_validation.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
out_root <- Sys.getenv("VALIDATION_OUT_ROOT", unset = file.path(tempdir(), paste0("phase0-", format(Sys.time(), "%Y%m%d-%H%M%S"))))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(VALIDATION_OUT_ROOT = out_root)
setwd(repo_root)

source(file.path(repo_root, "validation", "phase0", "build_validation_pack.R"), local = TRUE)
