#!/usr/bin/env bash
################################################################################
# Build the paper figures, tables and supplementary summaries in a new root.
# Requires: completed Parts 00-04 for all SSPs and 35 pooled geographies;
# per-city/GCM attribution and threshold sensitivity runs; extreme-heat-day
# runs; ERF uncertainty runs for SSP2 and SSP3; and the files checked below.
# Produces: figures, tables, checks, supplementary summaries and logs below
# FIGURE_ROOT. No source batch or sensitivity outputs are modified.
# How to run: BATCH_ROOT=results/europe/heatsplit_20261013 \
#   FIGURE_ROOT=results/paper_output_new GEO_SCOPE=main \
#   analysis/run_paper_outputs.sh
# Notes: SSP1 replacement and SSP1/SSP2 shape-check bypasses are exploratory;
#   see the caveats and prerequisite run commands in README.md.
################################################################################
set -euo pipefail

case "${1:-}" in
  "") check_only=0 ;;
  --check-inputs) check_only=1 ;;
  *) echo "ERROR: usage: $0 [--check-inputs]" >&2; exit 2 ;;
esac
[ "$#" -le 1 ] || { echo "ERROR: usage: $0 [--check-inputs]" >&2; exit 2; }

export BATCH_ROOT=${BATCH_ROOT:-results/europe/heatsplit_20261013}
export FIGURE_ROOT=${FIGURE_ROOT:?Set FIGURE_ROOT to a new output directory}
export GEO_SCOPE=${GEO_SCOPE:-main}
export SSP_SET=1,2,3
export N_HORIUCHI=400
export TABLE1_CORES=${TABLE1_CORES:-4}
export THRESHOLD_SUMMARY_CORES=${THRESHOLD_SUMMARY_CORES:-4}
export EXTREME_HEAT_DAYS_SUMMARY_CORES=${EXTREME_HEAT_DAYS_SUMMARY_CORES:-4}
export ATTRIBUTION_ROOT=${ATTRIBUTION_ROOT:-$BATCH_ROOT}
export THRESHOLD_ROOT=${THRESHOLD_ROOT:-results/europe/threshold_sensitivities_20261104}
export EXTREME_HEAT_DAYS_ROOT=${EXTREME_HEAT_DAYS_ROOT:-results/europe/extreme_heat_days_20261117/canonical_p975_p025}
export UNCERTAINTY_ROOT=${UNCERTAINTY_ROOT:-results/europe/uncertainty_20261110}

fail() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [ -s "$1" ] || fail "missing or empty input: $1"; }

case "$GEO_SCOPE" in main|all) ;; *) fail "GEO_SCOPE must be main or all";; esac
for value in "$TABLE1_CORES" "$THRESHOLD_SUMMARY_CORES" \
             "$EXTREME_HEAT_DAYS_SUMMARY_CORES"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || fail "core counts must be positive integers"
done
[ ! -e "$FIGURE_ROOT" ] && [ ! -L "$FIGURE_ROOT" ] ||
  fail "FIGURE_ROOT already exists; choose a new output directory: $FIGURE_ROOT"
output_path=$(realpath -m "$FIGURE_ROOT")
for root in "$BATCH_ROOT" "$ATTRIBUTION_ROOT" "$THRESHOLD_ROOT" \
            "$EXTREME_HEAT_DAYS_ROOT" "$UNCERTAINTY_ROOT"; do
  input_path=$(realpath -m "$root")
  case "$output_path/" in "$input_path/"*) fail "FIGURE_ROOT is inside input root $root";; esac
done

for f in data/city_results.csv data/wittgenstein_pop.csv data/wittgenstein_assr.csv \
         "$BATCH_ROOT/collected/object1_dataset.parquet" \
         "$BATCH_ROOT/collected/object2_contributions.parquet" \
         "$BATCH_ROOT/collected/object3_levels.parquet" \
         "$UNCERTAINTY_ROOT/checks/gcms.txt" \
         "$EXTREME_HEAT_DAYS_ROOT/shared/gcms.txt"; do
  require_file "$f"
done
command -v Rscript >/dev/null || fail "Rscript is required"
command -v soffice >/dev/null || fail "LibreOffice (soffice) is required for Table 1 PDFs"
Rscript -e 'p <- c("data.table", "arrow", "dplyr", "dlnm", "ungroup", "DemoDecomp", "ggplot2", "patchwork", "ggrepel", "flextable", "officer", "scales"); missing <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "))' ||
  fail "required R packages are unavailable"

city_text=$(Rscript -e 'x <- data.table::fread("data/city_results.csv"); cat(sort(unique(x$URAU_CODE)), sep="\n")') ||
  fail "could not read city metadata"
country_text=$(Rscript -e 'x <- data.table::fread("data/city_results.csv"); cat(sort(unique(x$CNTR_CODE)), sep="\n")') ||
  fail "could not read country metadata"
mapfile -t cities <<< "$city_text"
mapfile -t countries <<< "$country_text"
mapfile -t gcms < "$UNCERTAINTY_ROOT/checks/gcms.txt"
[ "${#cities[@]}" -eq 854 ] || fail "expected 854 cities in data/city_results.csv"
[ "${#countries[@]}" -eq 30 ] || fail "expected 30 countries in data/city_results.csv"
[ "${#gcms[@]}" -eq 19 ] || fail "expected 19 GCMs in uncertainty manifest"
expected_gcms=(
  ACCESS_CM2 ACCESS_ESM1_5 BCC_CSM2_MR CMCC_ESM2 CanESM5
  EC_Earth3 EC_Earth3_Veg_LR GFDL_ESM4 IITM_ESM INM_CM4_8 INM_CM5_0
  IPSL_CM6A_LR KACE_1_0_G MIROC6 MPI_ESM1_2_HR MPI_ESM1_2_LR
  MRI_ESM2_0 NorESM2_LM NorESM2_MM
)
diff -q <(printf '%s\n' "${gcms[@]}" | LC_ALL=C sort) \
        <(printf '%s\n' "${expected_gcms[@]}" | LC_ALL=C sort) >/dev/null ||
  fail "uncertainty manifest does not contain the expected 19 GCMs"
diff -q <(printf '%s\n' "${gcms[@]}" | LC_ALL=C sort) \
        <(LC_ALL=C sort "$EXTREME_HEAT_DAYS_ROOT/shared/gcms.txt") >/dev/null ||
  fail "uncertainty and extreme-heat-day GCM manifests differ"

for ssp in 1 2 3; do
  require_file "$BATCH_ROOT/ssp${ssp}/run_manifest.txt"
  require_file "$BATCH_ROOT/ssp${ssp}/run_gcms.txt"
  diff -q <(printf '%s\n' "${gcms[@]}" | LC_ALL=C sort) \
          <(LC_ALL=C sort "$BATCH_ROOT/ssp${ssp}/run_gcms.txt") >/dev/null ||
    fail "SSP${ssp} GCM manifest differs from the uncertainty manifest"
  for city in "${cities[@]}"; do
    require_file "$BATCH_ROOT/ssp${ssp}/$city/demography/00_demography_single_age.csv"
    require_file "$ATTRIBUTION_ROOT/ssp${ssp}/$city/ENSEMBLE/01_attribution_grouped.csv"
    if [ "$ssp" -eq 2 ]; then
      require_file "$BATCH_ROOT/ssp2/$city/ENSEMBLE/03_master_table.csv"
    fi
    for gcm in "${gcms[@]}"; do
      require_file "$BATCH_ROOT/ssp${ssp}/$city/$gcm/01_attribution_grouped.csv"
      require_file "$BATCH_ROOT/ssp${ssp}/$city/$gcm/checks/01_moving_heat_thresholds.csv"
      require_file "$EXTREME_HEAT_DAYS_ROOT/ssp${ssp}/$city/$gcm/01_extreme_heat_days.csv"
      for variant in trailing30_p975 trailing15_p95_p05; do
        require_file "$THRESHOLD_ROOT/$variant/ssp${ssp}/$city/$gcm/01_attribution_grouped.csv"
      done
    done
  done
  for spec in "europe/EUROPE" "region/Eastern" "region/Northern" \
              "region/Southern" "region/Western"; do
    for f in 03_master_table.csv 04_le_li_levels.csv \
             04_between_branch_decomposition.csv 04_within_branch_period_decomposition.csv; do
      require_file "$BATCH_ROOT/geographies/ssp${ssp}/$spec/$f"
    done
  done
  for country in "${countries[@]}"; do
    for f in 03_master_table.csv 04_le_li_levels.csv \
             04_between_branch_decomposition.csv 04_within_branch_period_decomposition.csv; do
      require_file "$BATCH_ROOT/geographies/ssp${ssp}/country/$country/$f"
    done
  done
done
require_file "$BATCH_ROOT/ssp3/ES001C/GFDL_ESM4/checks/01_attribution_checks.csv"
for ssp in 2 3; do
  for gcm in "${gcms[@]}"; do
    require_file "$UNCERTAINTY_ROOT/part01/ssp${ssp}/$gcm/experiment_metrics.parquet"
    require_file "$UNCERTAINTY_ROOT/part01/ssp${ssp}/$gcm/annual_risk.parquet"
  done
done

if [ "$check_only" -eq 1 ]; then
  echo "All paper-output inputs are present."
  exit 0
fi

mkdir -p "$FIGURE_ROOT/logs"
export TABLE1_SCRATCH="$FIGURE_ROOT/table1_scratch"
export EXTREME_OUTPUT_DIR="$FIGURE_ROOT/extreme_temperature"
export RESUME_FIGURES=0
analysis/run_paper_figures.sh || fail "paper figures or Table 1 failed; see $FIGURE_ROOT/logs/"

export TABLE1_FULL="$FIGURE_ROOT/table1/table1_full.csv"
export UNCERTAINTY_OUTPUT_ROOT="$FIGURE_ROOT/uncertainty"
for ssp in 2 3; do
  SSP="$ssp" Rscript analysis/summarize_erf_uncertainty.R \
    > "$FIGURE_ROOT/logs/uncertainty_ssp${ssp}.log" 2>&1 || {
      echo "ERROR: uncertainty summary failed; see $FIGURE_ROOT/logs/uncertainty_ssp${ssp}.log" >&2
      exit 1
    }
done

TABLE1_SSP3="$FIGURE_ROOT/table1/table1_ssp3_geographies.csv" \
UNCERTAINTY_SUMMARY="$UNCERTAINTY_OUTPUT_ROOT/tables/ssp3/summary.csv" \
  Rscript analysis/paper_table2.R > "$FIGURE_ROOT/logs/table2.log" 2>&1 || {
    echo "ERROR: Table 2 failed; see $FIGURE_ROOT/logs/table2.log" >&2
    exit 1
  }

export CENTRAL_EXTRAS_ROOT="$FIGURE_ROOT/central_extras"
TABLE1_FILE="$TABLE1_FULL" Rscript analysis/central_missing_results.R \
  > "$FIGURE_ROOT/logs/central_extras.log" 2>&1 || {
    echo "ERROR: central supplementary results failed; see $FIGURE_ROOT/logs/central_extras.log" >&2
    exit 1
  }
CANONICAL_ROOT="$BATCH_ROOT" Rscript analysis/summarize_threshold_sensitivity.R \
  > "$FIGURE_ROOT/logs/threshold_sensitivity.log" 2>&1 || {
    echo "ERROR: threshold sensitivity summary failed; see $FIGURE_ROOT/logs/threshold_sensitivity.log" >&2
    exit 1
  }
EXTREME_HEAT_DAYS_SUMMARY_ROOT="$FIGURE_ROOT/extreme_heat_days" \
  Rscript analysis/summarize_extreme_heat_days.R \
  > "$FIGURE_ROOT/logs/extreme_heat_days.log" 2>&1 || {
    echo "ERROR: extreme-heat-day summary failed; see $FIGURE_ROOT/logs/extreme_heat_days.log" >&2
    exit 1
  }
echo "Paper outputs complete: $FIGURE_ROOT"
