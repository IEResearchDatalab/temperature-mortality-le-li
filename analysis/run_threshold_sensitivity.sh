#!/usr/bin/env bash
set -euo pipefail

# Rerun only Part 01 for a moving-window / heat-percentile sensitivity. The
# canonical demographic outputs are read only. Usage:
#   analysis/run_threshold_sensitivity.sh LABEL HEAT_PROB WINDOW_YEARS [NCORES]

label=${1:?label}
heat_prob=${2:?heat percentile as probability, e.g. 0.95}
window_years=${3:?trailing window years}
ncores=${4:-28}
canonical=${CANONICAL_ROOT:-results/europe/heatsplit_20261013}
root=${THRESHOLD_ROOT:-results/europe/threshold_sensitivities_20261104}/${label}

[[ "$label" =~ ^[a-z0-9_]+$ ]] || { echo "invalid label"; exit 2; }
[[ "$ncores" =~ ^[1-9][0-9]*$ ]] || { echo "invalid core count"; exit 2; }
mkdir -p "$root/shared/checks" "$root/logs"

export EXTREME_HEAT_PROB=$heat_prob
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
      EXTREME_HEAT_PROB="$EXTREME_HEAT_PROB" \
      EXTREME_WINDOW_YEARS="$EXTREME_WINDOW_YEARS" \
      Rscript pipeline/01_attribution.R > "$out/run.log" 2>&1; then
    printf 'label=%s\nheat_prob=%s\nwindow_years=%s\nfinished=%s\n' \
      "$label" "$EXTREME_HEAT_PROB" "$EXTREME_WINDOW_YEARS" "$(date -Iseconds)" \
      > "$done_file"
  else
    printf '%s\t%s\t%s\n' "$ssp" "$city" "$gcm" >> "$root/failed.tsv"
    return 1
  fi
}
export -f run_threshold_job
export label root canonical PREP_DATA_FILE EXTREME_HEAT_PROB EXTREME_WINDOW_YEARS

echo "$(date '+%F %T') running $(wc -l < "$jobs") Part 01 jobs with $ncores workers"
set +e
xargs -a "$jobs" -n 3 -P "$ncores" bash -c 'run_threshold_job "$@"' _
status=$?
set -e

done_count=$(find "$root" -name .done_attribution -type f | wc -l)
job_count=$(wc -l < "$jobs")
printf 'label=%s\nheat_prob=%s\nwindow_years=%s\njobs=%s\ncompleted=%s\ngit_commit=%s\n' \
  "$label" "$heat_prob" "$window_years" "$job_count" "$done_count" "$(git rev-parse HEAD)" \
  > "$root/run_manifest.txt"
if [ "$status" -ne 0 ] || [ "$done_count" -ne "$job_count" ]; then
  echo "threshold run incomplete: $done_count/$job_count" >&2
  exit 1
fi
echo "$(date '+%F %T') completed $done_count/$job_count"

