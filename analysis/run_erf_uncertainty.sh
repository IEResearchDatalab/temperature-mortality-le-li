#!/usr/bin/env bash
set -euo pipefail

# Run one SSP across the 19 GCMs. Each GCM process pools all cities before
# life-table metrics are written. Usage: run_erf_uncertainty.sh SSP [NCORES]

ssp=${1:?SSP 2 or 3}
ncores=${2:-19}
case "$ssp" in 2|3) ;; *) echo "SSP must be 2 or 3" >&2; exit 2 ;; esac
root=${UNCERTAINTY_ROOT:-results/europe/uncertainty_20261104}
mkdir -p "$root/logs/ssp${ssp}" "$root/checks"

gcms=(
  ACCESS_CM2 ACCESS_ESM1_5 BCC_CSM2_MR CMCC_ESM2 CanESM5
  EC_Earth3 EC_Earth3_Veg_LR GFDL_ESM4 IITM_ESM INM_CM4_8 INM_CM5_0
  IPSL_CM6A_LR KACE_1_0_G MIROC6 MPI_ESM1_2_HR MPI_ESM1_2_LR
  MRI_ESM2_0 NorESM2_LM NorESM2_MM
)
printf '%s\n' "${gcms[@]}" > "$root/checks/gcms.txt"

run_gcm() {
  local gcm=$1
  local log="$root/logs/ssp${ssp}/${gcm}.log"
  if [ -s "$root/part01/ssp${ssp}/${gcm}/.done" ]; then
    echo "skip complete SSP${ssp} ${gcm}" >> "$log"
    return 0
  fi
  if SSP="$ssp" GCM="$gcm" DRAW_IDS=1-500 CITY_IDS=all \
      UNCERTAINTY_ROOT="$root" Rscript analysis/erf_uncertainty_pilot.R \
      > "$log" 2>&1; then
    return 0
  fi
  printf '%s\n' "$gcm" >> "$root/checks/ssp${ssp}_failed.txt"
  return 1
}
export -f run_gcm
export ssp root

printf '%s\n' "${gcms[@]}" | xargs -n 1 -P "$ncores" bash -c 'run_gcm "$1"' _
completed=$(find "$root/part01/ssp${ssp}" -mindepth 2 -maxdepth 2 -name .done -type f | wc -l)
if [ "$completed" -ne 19 ]; then
  echo "SSP${ssp} incomplete: $completed/19 GCMs" >&2
  exit 1
fi
printf 'ssp=%s\ngcms=19\ndraws=500\ncompleted=%s\ngit_commit=%s\n' \
  "$ssp" "$(date -Iseconds)" "$(git rev-parse HEAD)" \
  > "$root/ssp${ssp}_manifest.txt"
echo "SSP${ssp} complete: 19 GCMs x 500 draws"

