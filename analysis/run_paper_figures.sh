#!/usr/bin/env bash
################################################################################
# Render Simon's paper figures from completed pipeline outputs.
# Presentation orchestration only; scientific calculations stay in pipeline/.
################################################################################
set -euo pipefail

export BATCH_ROOT=${BATCH_ROOT:-results/europe/heatsplit_20261013}
export FIGURE_ROOT=${FIGURE_ROOT:-results/paper_figures/20261013}
export SSP_SET=${SSP_SET:-1,2,3}
export GEO_SCOPE=${GEO_SCOPE:-all}
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
  local level=$1 id=$2 log
  log="$FIGURE_ROOT/logs/${level}_${id}.log"
  echo "$(date '+%F %T') figures: $level/$id"
  env GEO_LEVEL="$level" GEO_ID="$id" Rscript analysis/paper_fig1.R > "$log" 2>&1
  env GEO_LEVEL="$level" GEO_ID="$id" Rscript analysis/paper_fig2.R >> "$log" 2>&1
  env GEO_LEVEL="$level" GEO_ID="$id" Rscript analysis/paper_fig3.R >> "$log" 2>&1
}

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

Rscript analysis/paper_table1.R > "$FIGURE_ROOT/logs/table1.log" 2>&1
echo "$(date '+%F %T') paper figures and Table 1 complete"
