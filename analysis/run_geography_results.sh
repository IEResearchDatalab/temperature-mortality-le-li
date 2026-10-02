#!/usr/bin/env bash
################################################################################
#
# Auxiliary results runner: Europe, regions and countries
#
#   Usage after run_pipeline.sh:
#     analysis/run_geography_results.sh <SSP> [NCORES]
#
# This script only creates presentation outputs from completed Part 04 files.
# It contains no scientific pipeline transformation.
#
################################################################################
set -uo pipefail

SSP=${1:?'SSP (1, 2 or 3)'}
NCORES=${2:-4}
export SSP
export BATCH_ROOT=${BATCH_ROOT_OVERRIDE:-results/europe}
export ROOT="$BATCH_ROOT/geographies/ssp${SSP}"

[ -d "$ROOT" ] || { echo "ERROR: missing $ROOT (run run_pipeline.sh first)"; exit 1; }
rm -f "$ROOT/failed.txt"

JOBS=$(mktemp)
trap 'rm -f "$JOBS"' EXIT
{
  echo 'europe EUROPE'
  Rscript -e 'x<-unique(data.table::fread("data/city_results.csv")$region); cat(paste("region",sort(x)),sep="\n")'
  Rscript -e 'x<-unique(data.table::fread("data/city_results.csv")$CNTR_CODE); cat(paste("country",sort(x)),sep="\n")'
} > "$JOBS"

results_job() { # results_job <level> <id>
  local level=$1 id=$2 d="$ROOT/$1/$2"
  [ -f "$d/.done" ] || return 0
  mkdir -p "$d/figures"
  if env GEO_LEVEL="$level" GEO_ID="$id" GCM=ENSEMBLE OUT_DIR="$d" \
      CHECK_DIR="$d/checks" FIG_DIR="$d/figures" SSP="$SSP" \
      Rscript analysis/geography_figures.R > "$d/results.log" 2>&1; then
    return 0
  else
    echo "$level $id" >> "$ROOT/failed.txt"
  fi
}
export -f results_job

n=$(wc -l < "$JOBS")
echo "$(date '+%F %T') result figures for $n geographies, $NCORES in parallel"
xargs -P "$NCORES" -L1 bash -c 'results_job "$0" "$1"' < "$JOBS"
if [ -f "$ROOT/failed.txt" ]; then
  nfail=$(wc -l < "$ROOT/failed.txt")
else
  nfail=0
fi
echo "$(date '+%F %T') done; failures: $nfail"
if [ "$nfail" -gt 0 ]; then
  echo "Failed geographies (see <level>/<id>/results.log):"
  cat "$ROOT/failed.txt"
  exit 1
fi

echo "$(date '+%F %T') building cross-geography overview"
if ! env GEOGRAPHY_ROOT="$ROOT" SSP="$SSP" Rscript analysis/geography_overview.R \
    > "$ROOT/overview.log" 2>&1; then
  echo "ERROR: cross-geography overview failed; see $ROOT/overview.log"
  exit 1
fi
echo "$(date '+%F %T') cross-geography overview complete"
