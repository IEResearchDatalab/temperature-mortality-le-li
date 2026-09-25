#!/usr/bin/env bash
################################################################################
#
# Batch runner: many cities x 19 GCMs x one SSP
#
#   Usage (from the repository root, after 00a_prep_temperature.R):
#     R_pipeline/run_batch.sh <SSP> <cities.txt | all> [NCORES]
#
#   Stage 1  Part 00 (demography), one job per city     -> <city>/demography/
#   Stage 2  Parts 01, 02, 03, 04 (levels only), one job per city x GCM
#                                                        -> <city>/<GCM>/
#            02/03 outputs are deleted afterwards (~10 MB each); the Part 01
#            ANs and Part 04 LE/LI levels are kept (per-GCM uncertainty)
#   Stage 3  Part 01b (ensemble-mean ANs), then 02, 03, 04 (decompositions),
#            one job per city                            -> <city>/ENSEMBLE/
#   Output root: results/europe/ssp<SSP>/. Jobs run in parallel (NCORES) across
#   cities and GCMs. Each finished job leaves a `.done` file, so an interrupted
#   run resumes where it stopped; failures are listed in <root>/failed.txt.
#   Part 04 skips the year-on-year decomposition unless DECOMP_ANNUAL=1, and
#   uses N_HORIUCHI=50 unless set. Collect results with R_pipeline/06_collect.R.
#
################################################################################
set -uo pipefail

SSP=${1:?"SSP (1, 2 or 3)"}
CITIES=${2:?"cities file (URAU code in the first column) or 'all'"}
NCORES=${3:-4}
export ROOT=results/europe/ssp${SSP} SSP
export DECOMP_ANNUAL=${DECOMP_ANNUAL:-0}
# Horiuchi steps: N = 50 as in Lloyd et al. (2024); on Madrid the results equal
# N = 400 to 6 decimals (closure error 8e-8), at 1/8 of the cost
export N_HORIUCHI=${N_HORIUCHI:-50}
mkdir -p "$ROOT"
rm -f "$ROOT/failed.txt"   # failures are re-evaluated on every run

if [ "$CITIES" = "all" ]; then
  CITY_LIST=$(Rscript -e 'cat(sort(unique(data.table::fread("data/city_results.csv")$URAU_CODE)), sep = "\n")')
else
  CITY_LIST=$(cut -d' ' -f1 "$CITIES")
fi
#----- Pre-flight: inputs and packages (fail fast with a clear message)
for f in data/city_results.csv data/coefs.csv data/wittgenstein_pop.csv data/wittgenstein_assr.csv \
         data/tmeanproj.gz.parquet data/era5series.gz.parquet data/prep_data.RData; do
  [ -s "$f" ] || { echo "ERROR: missing $f (see README 'Data'; prep_data.RData comes from 00a_prep_temperature.R)"; exit 1; }
done
Rscript -e 'suppressMessages(source("R_pipeline/00_pkg_params.R")); invisible(arrow::open_dataset("data/tmeanproj.gz.parquet")$schema)' \
  > "$ROOT/preflight.log" 2>&1 || { echo "ERROR: R packages or tmeanproj.gz.parquet not readable; see $ROOT/preflight.log"; cat "$ROOT/preflight.log"; exit 1; }
GCMS=$(Rscript -e 'suppressMessages(source("R_pipeline/00_pkg_params.R")); cat(gcmlist, sep = "\n")')

run_part() {  # run_part <city> <dir name> <gcm> <part> [ENV=value ...]
  local c=$1 dn=$2 g=$3 part=$4; shift 4
  local d=$ROOT/$c/$dn
  mkdir -p "$d/checks" "$d/figures"
  env CITY_ID="$c" SSP="$SSP" GCM="$g" OUT_DIR="$d" CHECK_DIR="$d/checks" FIG_DIR="$d/figures" \
      DEMOG_DIR="$ROOT/$c/demography" "$@" Rscript "R_pipeline/$part.R" >> "$d/log.txt" 2>&1
}

dem_job() {  # dem_job <city>
  local c=$1 d=$ROOT/$1/demography
  [ -f "$d/.done" ] && return 0
  run_part "$c" demography GFDL_ESM4 00_demography && touch "$d/.done" || echo "$c demography" >> "$ROOT/failed.txt"
}

gcm_job() {  # gcm_job <city> <gcm>
  local c=$1 g=$2 d=$ROOT/$1/$2
  [ -f "$d/.done" ] && return 0
  if [ ! -f "$ROOT/$c/demography/.done" ]; then echo "$c $g skipped: demography failed" >> "$ROOT/failed.txt"; return 0; fi
  if run_part "$c" "$g" "$g" 01_attribution && run_part "$c" "$g" "$g" 02_single_age &&
     run_part "$c" "$g" "$g" 03_master_table && run_part "$c" "$g" "$g" 04_le_li_decomposition LEVELS_ONLY=1; then
    rm -f "$d/02_single_age_an.csv" "$d/03_master_table.csv"
    touch "$d/.done"
  else
    echo "$c $g" >> "$ROOT/failed.txt"
  fi
}

ens_job() {  # ens_job <city>
  local c=$1 d=$ROOT/$1/ENSEMBLE
  [ -f "$d/.done" ] && return 0
  local n; n=$(ls "$ROOT/$c"/*/.done 2>/dev/null | grep -v -e /demography/ -e /ENSEMBLE/ | wc -l)
  if [ "$n" -lt 19 ]; then echo "$c ENSEMBLE skipped: only $n/19 GCMs done" >> "$ROOT/failed.txt"; return 0; fi
  if run_part "$c" ENSEMBLE ENSEMBLE 01b_ensemble && run_part "$c" ENSEMBLE ENSEMBLE 02_single_age &&
     run_part "$c" ENSEMBLE ENSEMBLE 03_master_table && run_part "$c" ENSEMBLE ENSEMBLE 04_le_li_decomposition; then
    touch "$d/.done"
  else
    echo "$c ENSEMBLE" >> "$ROOT/failed.txt"
  fi
}
export -f run_part dem_job gcm_job ens_job

echo "$(date '+%F %T') stage 1: demography"
printf '%s\n' $CITY_LIST | xargs -P "$NCORES" -I{} bash -c 'dem_job {}'
echo "$(date '+%F %T') stage 2: city x GCM"
for c in $CITY_LIST; do for g in $GCMS; do echo "$c $g"; done; done | xargs -P "$NCORES" -L1 bash -c 'gcm_job "$0" "$1"'
echo "$(date '+%F %T') stage 3: ensemble"
printf '%s\n' $CITY_LIST | xargs -P "$NCORES" -I{} bash -c 'ens_job {}'
nfail=$(cat "$ROOT/failed.txt" 2>/dev/null | wc -l)
echo "$(date '+%F %T') done; failures: $nfail"
if [ "$nfail" -gt 0 ]; then
  echo "First failures (details in $ROOT/<city>/<dir>/log.txt):"; head -5 "$ROOT/failed.txt"
fi
