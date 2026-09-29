#!/usr/bin/env bash
################################################################################
#
# Pooled analysis runner: Europe, four regions and all countries
#
#   Usage (from the repository root, after pipeline/06_collect.R):
#     pipeline/run_pooled.sh <SSP> [NCORES]
#
#   Each job runs Parts 07 and 08 for one geography. Counts are read from the
#   collected Object 1 and summed before life tables and decompositions are
#   constructed, following the aggregation decision from 25 Sep 2026.
#
################################################################################
set -uo pipefail

SSP=${1:?'SSP (1, 2 or 3)'}
NCORES=${2:-4}
export SSP N_HORIUCHI=${N_HORIUCHI:-400}
export ROOT=results/europe/pooled/ssp${SSP}
export OBJECT1_FILE=${OBJECT1_FILE:-results/europe/collected/object1_dataset.parquet}

[ -s "$OBJECT1_FILE" ] || { echo "ERROR: missing $OBJECT1_FILE (run pipeline/06_collect.R first)"; exit 1; }
mkdir -p "$ROOT"
rm -f "$ROOT/failed.txt"

JOBS=$(mktemp)
trap 'rm -f "$JOBS"' EXIT
{
  echo 'europe EUROPE'
  Rscript -e 'x<-unique(data.table::fread("data/city_results.csv")$region); cat(paste("region",sort(x)),sep="\n")'
  Rscript -e 'x<-unique(data.table::fread("data/city_results.csv")$CNTR_CODE); cat(paste("country",sort(x)),sep="\n")'
} > "$JOBS"

pooled_job() { # pooled_job <level> <id>
  local level=$1 id=$2 d="$ROOT/$1/$2"
  [ -f "$d/.done" ] && return 0
  mkdir -p "$d/checks" "$d/figures"
  if env GEO_LEVEL="$level" GEO_ID="$id" GCM=ENSEMBLE OUT_DIR="$d" \
      CHECK_DIR="$d/checks" FIG_DIR="$d/figures" OBJECT1_FILE="$OBJECT1_FILE" \
      SSP="$SSP" N_HORIUCHI="$N_HORIUCHI" Rscript pipeline/07_pooled_le_li.R \
      > "$d/log.txt" 2>&1 && \
     env GEO_LEVEL="$level" GEO_ID="$id" GCM=ENSEMBLE OUT_DIR="$d" \
      CHECK_DIR="$d/checks" FIG_DIR="$d/figures" SSP="$SSP" \
      Rscript pipeline/08_pooled_figures.R >> "$d/log.txt" 2>&1; then
    touch "$d/.done"
  else
    echo "$level $id" >> "$ROOT/failed.txt"
  fi
}
export -f pooled_job

n=$(wc -l < "$JOBS")
echo "$(date '+%F %T') pooled analysis for $n geographies, $NCORES in parallel, Horiuchi N=$N_HORIUCHI"
xargs -P "$NCORES" -L1 bash -c 'pooled_job "$0" "$1"' < "$JOBS"
if [ -f "$ROOT/failed.txt" ]; then
  nfail=$(wc -l < "$ROOT/failed.txt")
else
  nfail=0
fi
echo "$(date '+%F %T') done; failures: $nfail"
if [ "$nfail" -gt 0 ]; then
  echo "Failed geographies (see <level>/<id>/log.txt):"
  cat "$ROOT/failed.txt"
  exit 1
fi

echo "$(date '+%F %T') building cross-geography overview"
if ! env POOLED_ROOT="$ROOT" SSP="$SSP" Rscript pipeline/09_pooled_overview.R \
    > "$ROOT/overview.log" 2>&1; then
  echo "ERROR: cross-geography overview failed; see $ROOT/overview.log"
  exit 1
fi
echo "$(date '+%F %T') cross-geography overview complete"

