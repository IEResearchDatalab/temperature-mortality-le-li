#!/usr/bin/env bash
################################################################################
#
# Runner for the five-step scientific pipeline
#
#   Usage from the repository root:
#     ./run_pipeline.sh <SSP> <cities.txt | all> [NCORES]
#
#   00_prep_data.R             temperature thresholds and demography
#   01_attribution.R           GCM-specific and ensemble attributable deaths
#   02_single_age.R            allocation to single ages
#   03_master_table.R          city and pooled analysis datasets (Object 1)
#   04_le_li_decomposition.R   LE/LI levels and decompositions (Objects 2-3)
#
# The runner only orchestrates these scripts. It contains no scientific
# transformation. Jobs are resumable through `.done` markers.
#
################################################################################
set -uo pipefail

SSP=${1:?"SSP (1, 2 or 3)"}
CITIES=${2:?"cities file (URAU code in the first column) or 'all'"}
NCORES=${3:-4}
export ROOT=results/europe/ssp${SSP} SSP
export BATCH_ROOT=results/europe
export DECOMP_ANNUAL=${DECOMP_ANNUAL:-0}
# Horiuchi steps: N = 50 as in Lloyd et al. (2024); on Madrid the results equal
# N = 400 to 6 decimals (closure error 8e-8), at 1/8 of the cost
export N_HORIUCHI=${N_HORIUCHI:-50}
mkdir -p "$ROOT" "$BATCH_ROOT/shared/checks"
rm -f "$ROOT/failed.txt"   # failures are re-evaluated on every run

if [ "$CITIES" = "all" ]; then
  CITY_LIST=$(Rscript -e 'cat(sort(unique(data.table::fread("data/city_results.csv")$URAU_CODE)), sep = "\n")')
else
  CITY_LIST=$(cut -d' ' -f1 "$CITIES")
fi
#----- Pre-flight: inputs and packages (fail fast with a clear message)
for f in data/city_results.csv data/coefs.csv data/wittgenstein_pop.csv data/wittgenstein_assr.csv \
         data/tmeanproj.gz.parquet data/era5series.gz.parquet; do
  [ -s "$f" ] || { echo "ERROR: missing $f (see README 'Data')"; exit 1; }
done
Rscript -e 'suppressMessages(source("pipeline/00_pkg_params.R")); invisible(arrow::open_dataset("data/tmeanproj.gz.parquet")$schema)' \
  > "$ROOT/preflight.log" 2>&1 || { echo "ERROR: R packages or tmeanproj.gz.parquet not readable; see $ROOT/preflight.log"; cat "$ROOT/preflight.log"; exit 1; }
GCMS=$(Rscript -e 'suppressMessages(source("pipeline/00_pkg_params.R")); cat(gcmlist, sep = "\n")')

run_part() {  # run_part <city> <dir name> <gcm> <part> [ENV=value ...]
  local c=$1 dn=$2 g=$3 part=$4; shift 4
  local d=$ROOT/$c/$dn
  mkdir -p "$d/checks"
  env CITY_ID="$c" GEO_LEVEL=city GEO_ID="$c" SSP="$SSP" GCM="$g" OUT_DIR="$d" CHECK_DIR="$d/checks" \
      DEMOG_DIR="$ROOT/$c/demography" "$@" Rscript "pipeline/$part.R" >> "$d/log.txt" 2>&1
}

dem_job() {  # dem_job <city>
  local c=$1 d=$ROOT/$1/demography
  [ -f "$d/.done" ] && return 0
  run_part "$c" demography GFDL_ESM4 00_prep_data PREP_SCOPE=demography && touch "$d/.done" || echo "$c demography" >> "$ROOT/failed.txt"
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
  if run_part "$c" ENSEMBLE ENSEMBLE 01_attribution && run_part "$c" ENSEMBLE ENSEMBLE 02_single_age &&
     run_part "$c" ENSEMBLE ENSEMBLE 03_master_table && run_part "$c" ENSEMBLE ENSEMBLE 04_le_li_decomposition; then
    touch "$d/.done"
  else
    echo "$c ENSEMBLE" >> "$ROOT/failed.txt"
  fi
}
geo_job() {  # geo_job <level> <id>
  local level=$1 id=$2 label d
  case "$level" in
    europe) label=Europe ;;
    region) label="$id Europe" ;;
    country) label="$id" ;;
  esac
  d="$BATCH_ROOT/geographies/ssp$SSP/$level/$id"
  [ -f "$d/.done" ] && return 0
  mkdir -p "$d/checks"
  if env ANALYSIS_MODE=geography GEO_LEVEL="$level" GEO_ID="$id" GEO_LABEL="$label" \
      SSP="$SSP" GCM=ENSEMBLE OUT_DIR="$d" CHECK_DIR="$d/checks" \
      OBJECT1_FILE="$BATCH_ROOT/collected/object1_dataset.parquet" \
      Rscript pipeline/03_master_table.R > "$d/log.txt" 2>&1 && \
     env ANALYSIS_MODE=run GEO_LEVEL="$level" GEO_ID="$id" GEO_LABEL="$label" \
      SSP="$SSP" GCM=ENSEMBLE OUT_DIR="$d" CHECK_DIR="$d/checks" \
      DECOMP_ANNUAL="$DECOMP_ANNUAL" N_HORIUCHI="$N_HORIUCHI" \
      Rscript pipeline/04_le_li_decomposition.R >> "$d/log.txt" 2>&1; then
    touch "$d/.done"
  else
    echo "$level $id" >> "$ROOT/failed.txt"
  fi
}
export -f run_part dem_job gcm_job ens_job geo_job

if [ ! -s data/prep_data.RData ]; then
  echo "$(date '+%F %T') step 00: observed temperature and thresholds"
  env PREP_SCOPE=temperature CHECK_DIR="$BATCH_ROOT/shared/checks" \
    Rscript pipeline/00_prep_data.R > "$BATCH_ROOT/shared/temperature.log" 2>&1 || {
      echo "ERROR: temperature preparation failed; see $BATCH_ROOT/shared/temperature.log"; exit 1;
    }
fi

echo "$(date '+%F %T') step 00: demography"
printf '%s\n' $CITY_LIST | xargs -P "$NCORES" -I{} bash -c 'dem_job {}'
echo "$(date '+%F %T') steps 01-04: city x GCM"
for c in $CITY_LIST; do for g in $GCMS; do echo "$c $g"; done; done | xargs -P "$NCORES" -L1 bash -c 'gcm_job "$0" "$1"'
echo "$(date '+%F %T') steps 01-04: city ensemble"
printf '%s\n' $CITY_LIST | xargs -P "$NCORES" -I{} bash -c 'ens_job {}'

echo "$(date '+%F %T') step 03: collect Object 1"
mkdir -p "$BATCH_ROOT/collected"
env ANALYSIS_MODE=collect BATCH_ROOT="$BATCH_ROOT" Rscript pipeline/03_master_table.R \
  > "$BATCH_ROOT/collected/03_collect.log" 2>&1 || {
    echo "ERROR: Object 1 collection failed"; exit 1;
  }

if [ "$CITIES" = "all" ]; then
  JOBS=$(mktemp)
  trap 'rm -f "$JOBS"' EXIT
  {
    echo 'europe EUROPE'
    Rscript -e 'x<-unique(data.table::fread("data/city_results.csv")$region); cat(paste("region",sort(x)),sep="\n")'
    Rscript -e 'x<-unique(data.table::fread("data/city_results.csv")$CNTR_CODE); cat(paste("country",sort(x)),sep="\n")'
  } > "$JOBS"
  echo "$(date '+%F %T') steps 03-04: country, region and Europe"
  xargs -P "$NCORES" -L1 bash -c 'geo_job "$0" "$1"' < "$JOBS"
else
  echo "Skipping pooled geographies because this is a partial city run."
fi

echo "$(date '+%F %T') step 04: collect Objects 2 and 3"
env ANALYSIS_MODE=collect BATCH_ROOT="$BATCH_ROOT" Rscript pipeline/04_le_li_decomposition.R \
  > "$BATCH_ROOT/collected/04_collect.log" 2>&1 || {
    echo "ERROR: Objects 2 and 3 collection failed"; exit 1;
  }

nfail=$(cat "$ROOT/failed.txt" 2>/dev/null | wc -l)
echo "$(date '+%F %T') done; failures: $nfail"
if [ "$nfail" -gt 0 ]; then
  echo "First failures (details in $ROOT/<city>/<dir>/log.txt):"; head -5 "$ROOT/failed.txt"
fi
