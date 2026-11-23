#!/usr/bin/env bash
################################################################################
# Render Simon's paper figures from completed pipeline outputs.
# Presentation orchestration only; scientific calculations stay in pipeline/.
# Requires: Parts 00-04 and collected objects for the selected SSPs in
# BATCH_ROOT; pooled 03/04 CSVs, city Part 01/00 files and moving-threshold
# checks. Figure 2 sources analysis/decade_decomposition.R.
# Produces: FIGURE_ROOT figure PNGs/PDFs and data, Table 1, threshold outputs
# and logs; TABLE1_SCRATCH holds per-GCM and decade caches.
# How to run: BATCH_ROOT=results/europe/heatsplit_20261013 \
#   FIGURE_ROOT=results/paper_output_new \
#   TABLE1_SCRATCH=results/paper_output_new/table1_scratch GEO_SCOPE=main \
#   analysis/run_paper_figures.sh
# Notes: SSP1 replacement and SSP1/SSP2 shape-check bypasses are exploratory;
#   see README.md. This runner does not itself run Table 2 or uncertainty.
################################################################################
set -euo pipefail

export BATCH_ROOT=${BATCH_ROOT:-results/europe/heatsplit_20261013}
export FIGURE_ROOT=${FIGURE_ROOT:-results/figures_simon_20oct}
export SSP_SET=${SSP_SET:-1,2,3}
export GEO_SCOPE=${GEO_SCOPE:-all}
export RESUME_FIGURES=${RESUME_FIGURES:-0}
mkdir -p "$FIGURE_ROOT/logs"

IFS=',' read -r -a requested_ssps <<< "$SSP_SET"
for ssp in "${requested_ssps[@]}"; do
  [ -f "$BATCH_ROOT/ssp${ssp}/run_manifest.txt" ] || {
    echo "ERROR: missing completed SSP${ssp} manifest under $BATCH_ROOT" >&2
    exit 1
  }
done
for f in object1_dataset.parquet object2_contributions.parquet object3_levels.parquet; do
  [ -s "$BATCH_ROOT/collected/$f" ] || {
    echo "ERROR: missing collected $f under $BATCH_ROOT" >&2
    exit 1
  }
done

run_geo() {
  local level=$1 id=$2 log geo
  log="$FIGURE_ROOT/logs/${level}_${id}.log"
  geo=$(printf '%s' "$id" | tr '[:upper:]' '[:lower:]')
  if [ "$RESUME_FIGURES" = "1" ] &&
     [ -s "$FIGURE_ROOT/figures_pdf/figure02_le65_contributions_$geo.pdf" ] &&
     [ -s "$FIGURE_ROOT/data/figure02_le65_contributions_panel_c_data_$geo.csv" ] &&
     [ -s "$FIGURE_ROOT/figures_pdf/figure03_li65_contributions_$geo.pdf" ] &&
     [ -s "$FIGURE_ROOT/figures_pdf/figure03_supplement_five_ranges_$geo.pdf" ] &&
     [ -s "$FIGURE_ROOT/figures_pdf/figure03_supplement_plane_steps_$geo.pdf" ]; then
    echo "$(date '+%F %T') figures already complete: $level/$id"
    return
  fi
  echo "$(date '+%F %T') figures: $level/$id"
  env GEO_LEVEL="$level" GEO_ID="$id" Rscript analysis/paper_fig2.R > "$log" 2>&1
  env GEO_LEVEL="$level" GEO_ID="$id" Rscript analysis/paper_fig3.R >> "$log" 2>&1
}

# Figure 1 pools exact-band city inputs for Europe and the four regions in one run.
echo "$(date '+%F %T') figure 1: Europe and regions"
env GEO_LEVEL=all Rscript analysis/paper_fig1.R > "$FIGURE_ROOT/logs/figure1.log" 2>&1

# Europe first, then four regions, then 30 countries for supplementary material.
run_geo europe EUROPE
for id in $(find "$BATCH_ROOT/geographies/ssp3/region" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); do
  run_geo region "$id"
done
if [ "$GEO_SCOPE" = "all" ]; then
  for id in $(find "$BATCH_ROOT/geographies/ssp3/country" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); do
    run_geo country "$id"
  done
elif [ "$GEO_SCOPE" != "main" ]; then
  echo "ERROR: GEO_SCOPE must be all or main" >&2
  exit 1
fi

TABLE1_DIR="$FIGURE_ROOT/table1" TABLE1_STAGE=all Rscript analysis/paper_table1.R > "$FIGURE_ROOT/logs/table1.log" 2>&1
Rscript analysis/paper_table_comparison.R > "$FIGURE_ROOT/logs/table_comparison.log" 2>&1
EXTREME_OUTPUT_DIR="${EXTREME_OUTPUT_DIR:-$FIGURE_ROOT/extreme_temperature}" \
  Rscript analysis/extreme_temperature_thresholds.R \
    > "$FIGURE_ROOT/logs/extreme_temperature_thresholds.log" 2>&1
echo "$(date '+%F %T') paper figures and Table 1 complete"
