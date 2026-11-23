# Output file naming shared by the paper-figure scripts.
# Requires: no earlier outputs; source this file from the repository root.
# Produces: figure_file(), which creates FIGURE_ROOT/{figures,figures_pdf,data}.
# How to run: Rscript -e 'source("analysis/output_names.R")'
# Notes: figure content and the SSP caveats are supplied by the calling script.
#
# Convention: figureNN_<name>_<geography>.<ext>
#   figure02_le65_contributions_europe.png
#   figure03_li65_contributions_bars_data_western.csv
# Files are routed by type under FIGURE_ROOT:
#   .png -> figures/   .pdf -> figures_pdf/   anything else (csv, txt) -> data/
# The geography id (lower case) goes just before the extension.
figure_file <- function(figure_root, filename, geo_id) {
  ext <- tools::file_ext(filename)
  base <- tools::file_path_sans_ext(filename)
  dir_name <- switch(tolower(ext), png = "figures", pdf = "figures_pdf", "data")
  dir <- file.path(figure_root, dir_name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(base, "_", tolower(geo_id), ".", ext))
}
