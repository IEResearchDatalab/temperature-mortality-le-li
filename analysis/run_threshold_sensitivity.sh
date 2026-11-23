#!/usr/bin/env bash
# Rerun Part 01 using an alternative threshold or warning-day-only mode.
# Requires: data/city_results.csv, data/era5series.gz.parquet,
# data/tmeanproj.gz.parquet, data/coefs.csv, and Parts 00-04 under
# CANONICAL_ROOT, including ssp<SSP>/<city>/demography/00_demography_grouped.csv.
# Produces: THRESHOLD_ROOT/<label>/shared/, ssp<SSP>/<city>/<GCM>/, checks and logs.
# How to run: CANONICAL_ROOT=results/europe/heatsplit_20261013 \
#   THRESHOLD_ROOT=results/paper_threshold_new \
#   analysis/run_threshold_sensitivity.sh trailing15_p95_p05 0.95 15 28 0.05
# Notes: SSP1 replacement and SSP1/SSP2 shape-check bypasses are exploratory;
#   the runner reuses their demography. Set EXTREME_HEAT_DAYS_ONLY=1 for the
#   warning-day input run (see README.md).
set -euo pipefail

# Rerun only Part 01 for a moving-window / threshold-percentile sensitivity. The
# canonical demographic outputs are read only. Usage:
#   analysis/run_threshold_sensitivity.sh LABEL HEAT_PROB WINDOW_YEARS [NCORES] [COLD_PROB]
# Example for the paired P95/P5 sensitivity requested for the supplement:
#   analysis/run_threshold_sensitivity.sh trailing15_p95_p05 0.95 15 28 0.05

label=${1:?label}
heat_prob=${2:?heat percentile as probability, e.g. 0.95}
window_years=${3:?trailing window years}
ncores=${4:-28}
cold_prob=${5:-0.025}
canonical=${CANONICAL_ROOT:-results/europe/heatsplit_20261013}
root=${THRESHOLD_ROOT:-results/europe/threshold_sensitivities_20261104}/${label}
pipeline_01=${PIPELINE_01_SCRIPT:-pipeline/01_attribution.R}

[[ "$label" =~ ^[a-z0-9_]+$ ]] || { echo "invalid label"; exit 2; }
[[ "$ncores" =~ ^[1-9][0-9]*$ ]] || { echo "invalid core count"; exit 2; }
mkdir -p "$root/shared/checks" "$root/logs"

export EXTREME_HEAT_PROB=$heat_prob
export EXTREME_COLD_PROB=$cold_prob
export EXTREME_WINDOW_YEARS=$window_years
export PREP_DATA_FILE="$root/shared/prep_data.RData"

if [ ! -s "$PREP_DATA_FILE" ]; then
  echo "$(date '+%F %T') building temperature preparation cache"
  PREP_SCOPE=temperature CHECK_DIR="$root/shared/checks" \
    Rscript pipeline/00_prep_data.R > "$root/logs/temperature.log" 2>&1
fi

Rscript -e 'cat(sort(unique(data.table::fread("data/city_results.csv")$URAU_CODE)), sep="\n")' \
  > "$root/shared/cities.txt"
Rscript -e 'suppressMessages(source("pipeline/00_pkg_params.R")); cat(gcmlist, sep="\n")' \
  > "$root/shared/gcms.txt"

jobs="$root/shared/jobs.tsv"
: > "$jobs"
for ssp in 1 2 3; do
  while read -r city; do
    while read -r gcm; do
      printf '%s\t%s\t%s\n' "$ssp" "$city" "$gcm" >> "$jobs"
    done < "$root/shared/gcms.txt"
  done < "$root/shared/cities.txt"
done

run_threshold_job() {
  local ssp=$1 city=$2 gcm=$3
  local out="$root/ssp${ssp}/${city}/${gcm}"
  local done_file="$out/.done_attribution"
  [ -s "$done_file" ] && return 0
  mkdir -p "$out/checks"
  if CITY_ID="$city" GEO_ID="$city" SSP="$ssp" GCM="$gcm" \
      OUT_DIR="$out" CHECK_DIR="$out/checks" \
      DEMOG_DIR="$canonical/ssp${ssp}/${city}/demography" \
      PREP_DATA_FILE="$PREP_DATA_FILE" \
      EXTREME_COLD_PROB="$EXTREME_COLD_PROB" \
      EXTREME_HEAT_PROB="$EXTREME_HEAT_PROB" \
      EXTREME_WINDOW_YEARS="$EXTREME_WINDOW_YEARS" \
      Rscript "$pipeline_01" > "$out/run.log" 2>&1; then
    printf 'label=%s\ncold_prob=%s\nheat_prob=%s\nwindow_years=%s\nfinished=%s\n' \
      "$label" "$EXTREME_COLD_PROB" "$EXTREME_HEAT_PROB" "$EXTREME_WINDOW_YEARS" "$(date -Iseconds)" \
      > "$done_file"
  else
    printf '%s\t%s\t%s\n' "$ssp" "$city" "$gcm" >> "$root/failed.tsv"
    return 1
  fi
}
export -f run_threshold_job
export label root canonical pipeline_01 PREP_DATA_FILE EXTREME_COLD_PROB EXTREME_HEAT_PROB EXTREME_WINDOW_YEARS

echo "$(date '+%F %T') running $(wc -l < "$jobs") Part 01 jobs with $ncores workers"
set +e
xargs -a "$jobs" -n 3 -P "$ncores" bash -c 'run_threshold_job "$@"' _
status=$?
set -e

done_count=$(find "$root" -name .done_attribution -type f | wc -l)
job_count=$(wc -l < "$jobs")
pipeline_commit=$(git -C "$(dirname "$pipeline_01")/.." rev-parse HEAD 2>/dev/null || git rev-parse HEAD)
printf 'label=%s\ncold_prob=%s\nheat_prob=%s\nwindow_years=%s\njobs=%s\ncompleted=%s\npipeline_commit=%s\n' \
  "$label" "$cold_prob" "$heat_prob" "$window_years" "$job_count" "$done_count" "$pipeline_commit" \
  > "$root/run_manifest.txt"
if [ "$status" -ne 0 ] || [ "$done_count" -ne "$job_count" ]; then
  echo "threshold run incomplete: $done_count/$job_count" >&2
  exit 1
fi
echo "$(date '+%F %T') completed $done_count/$job_count"

