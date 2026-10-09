#!/usr/bin/env bash
################################################################################
#
# Runner for the five-step scientific pipeline
#
#   Usage from the repository root:
#     ./run_pipeline.sh <SSP> <cities.txt | all> [NCORES]
#     ./run_pipeline.sh --config run_config.csv
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

if [ "${1:-}" = "--config" ]; then
  [ "$#" -eq 2 ] || { echo "ERROR: usage: $0 --config <run_config.csv>"; exit 2; }
  RUN_CONFIG=$2
  [ -s "$RUN_CONFIG" ] || { echo "ERROR: missing or empty run config: $RUN_CONFIG"; exit 2; }

  IFS= read -r header < "$RUN_CONFIG"
  header=${header%$'\r'}
  expected_header='ssp,enabled,cities,ncores,n_horiuchi,decomp_annual'
  [ "$header" = "$expected_header" ] || {
    echo "ERROR: config header must be: $expected_header"
    exit 2
  }

  nruns=0
  while IFS=, read -r run_ssp enabled cities ncores n_horiuchi decomp_annual extra; do
    run_ssp=${run_ssp%$'\r'}
    enabled=${enabled%$'\r'}
    cities=${cities%$'\r'}
    ncores=${ncores%$'\r'}
    n_horiuchi=${n_horiuchi%$'\r'}
    decomp_annual=${decomp_annual%$'\r'}
    [ -z "$run_ssp$enabled$cities$ncores$n_horiuchi$decomp_annual${extra:-}" ] && continue
    [ -z "${extra:-}" ] || { echo "ERROR: too many fields in config row for SSP $run_ssp"; exit 2; }
    case "$enabled" in
      true) ;;
      false) continue ;;
      *) echo "ERROR: enabled must be true or false for SSP $run_ssp"; exit 2 ;;
    esac
    case "$decomp_annual" in
      0|1) ;;
      *) echo "ERROR: decomp_annual must be 0 or 1 for SSP $run_ssp"; exit 2 ;;
    esac

    nruns=$((nruns + 1))
    echo "$(date '+%F %T') config run $nruns: SSP$run_ssp, cities=$cities, cores=$ncores"
    if ! RUN_CONFIG_SOURCE="$RUN_CONFIG" N_HORIUCHI="$n_horiuchi" DECOMP_ANNUAL="$decomp_annual" \
        "$0" "$run_ssp" "$cities" "$ncores"; then
      echo "ERROR: config stopped after SSP$run_ssp failed. Later scenarios were not started."
      exit 1
    fi
  done < <(tail -n +2 "$RUN_CONFIG" | tr -d '\r')

  [ "$nruns" -gt 0 ] || { echo "ERROR: the config contains no enabled runs"; exit 2; }
  echo "$(date '+%F %T') all $nruns configured scenario runs completed"
  exit 0
fi

SSP=${1:?"SSP (1, 2 or 3)"}
CITIES=${2:?"cities file (URAU code in the first column) or 'all'"}
NCORES=${3:-4}
case "$SSP" in 1|2|3) ;; *) echo "ERROR: SSP must be 1, 2 or 3"; exit 2 ;; esac
[[ "$NCORES" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: NCORES must be a positive integer"; exit 2; }
[[ "${N_HORIUCHI:-50}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: N_HORIUCHI must be a positive integer"; exit 2; }
case "${DECOMP_ANNUAL:-0}" in 0|1) ;; *) echo "ERROR: DECOMP_ANNUAL must be 0 or 1"; exit 2 ;; esac
[ "$CITIES" = "all" ] || [ -s "$CITIES" ] || { echo "ERROR: missing or empty cities file: $CITIES"; exit 2; }
case "${SKIP_MX_PLAUSIBILITY_CHECK:-0}" in 0|1) ;; *) echo "ERROR: SKIP_MX_PLAUSIBILITY_CHECK must be 0 or 1"; exit 2 ;; esac
if [ -n "${ASSR_ONE_REPLACEMENT:-}" ]; then
  awk -v x="$ASSR_ONE_REPLACEMENT" 'BEGIN { exit !(x + 0 == x && x > 0 && x < 1) }' || {
    echo "ERROR: ASSR_ONE_REPLACEMENT must be a number strictly between 0 and 1"
    exit 2
  }
fi

# WARNING: Both switches alter a scientific safeguard or source value. They
# are sensitivity analyses only and must use a separate batch root so canonical
# results/europe outputs cannot be overwritten or silently mixed with them.
# Paths are compared after resolving symlinks, "..", "." and trailing slashes.
# The override must resolve outside results/europe, or strictly below
# results/europe/sensitivities/ (the convention used by existing runs).
if [ "${SKIP_MX_PLAUSIBILITY_CHECK:-0}" = "1" ] || [ -n "${ASSR_ONE_REPLACEMENT:-}" ]; then
  [ -n "${BATCH_ROOT_OVERRIDE:-}" ] || {
    echo "ERROR: sensitivity switches require a non-canonical BATCH_ROOT_OVERRIDE"
    exit 2
  }
  canonical_root=$(realpath -m results/europe)
  override_root=$(realpath -m "$BATCH_ROOT_OVERRIDE")
  case "$override_root/" in
    "$canonical_root"/sensitivities/?*) ;;
    "$canonical_root"/*)
      echo "ERROR: BATCH_ROOT_OVERRIDE resolves into the canonical results tree: $override_root"
      echo "       Use a path outside $canonical_root or below $canonical_root/sensitivities/"
      exit 2 ;;
  esac
fi

export BATCH_ROOT=${BATCH_ROOT_OVERRIDE:-results/europe}
export ROOT="$BATCH_ROOT/ssp${SSP}" SSP
export SKIP_MX_PLAUSIBILITY_CHECK=${SKIP_MX_PLAUSIBILITY_CHECK:-0}
export ASSR_ONE_REPLACEMENT=${ASSR_ONE_REPLACEMENT:-}
export DECOMP_ANNUAL=${DECOMP_ANNUAL:-0}
# Horiuchi steps: N = 50 as in Lloyd et al. (2024); on Madrid the results equal
# N = 400 to 6 decimals (closure error 8e-8), at 1/8 of the cost
export N_HORIUCHI=${N_HORIUCHI:-50}

#----- Pre-flight: inputs and packages (fail fast with a clear message)
# Nothing is created or deleted until these checks pass.
for f in data/city_results.csv data/coefs.csv data/wittgenstein_pop.csv data/wittgenstein_assr.csv \
         data/tmeanproj.gz.parquet data/era5series.gz.parquet; do
  [ -s "$f" ] || { echo "ERROR: missing $f (see README 'Data')"; exit 1; }
done
if [ "$CITIES" = "all" ]; then
  CITY_LIST=$(Rscript -e 'cat(sort(unique(data.table::fread("data/city_results.csv")$URAU_CODE)), sep = "\n")')
else
  CITY_LIST=$(tr -d '\r' < "$CITIES" | awk 'NF { print $1 }' | sort -u)
fi
[ -n "$CITY_LIST" ] || { echo "ERROR: the city list is empty"; exit 1; }
# City IDs are later expanded unquoted and passed to bash -c, so accept only
# URAU codes (e.g. ES001C); anything else could be globbed or run as shell code.
bad_cities=$(printf '%s\n' "$CITY_LIST" | grep -v -x -E '[A-Z]{2}[0-9]{3}C')
if [ -n "$bad_cities" ]; then
  echo "ERROR: invalid city IDs in $CITIES (expected URAU codes such as ES001C):"
  printf '%s\n' "$bad_cities" | head -5 | sed 's/^/  /'
  exit 2
fi
mkdir -p "$ROOT" "$BATCH_ROOT/shared/checks"
Rscript -e 'suppressMessages(source("pipeline/00_pkg_params.R")); invisible(arrow::open_dataset("data/tmeanproj.gz.parquet")$schema)' \
  > "$ROOT/preflight.log" 2>&1 || { echo "ERROR: R packages or tmeanproj.gz.parquet not readable; see $ROOT/preflight.log"; cat "$ROOT/preflight.log"; exit 1; }
GCMS=$(Rscript -e 'suppressMessages(source("pipeline/00_pkg_params.R")); cat(gcmlist, sep = "\n")')
[ -n "$GCMS" ] || { echo "ERROR: could not read the GCM list from 00_pkg_params.R"; exit 1; }

#----- Provenance of resumed outputs
# Finished jobs are skipped through .done markers, so rerunning into an
# existing scenario root reuses its outputs. Resume only when the scientific
# settings and the pipeline/ code match the run recorded in run_manifest.txt;
# otherwise stop before anything is deleted or rewritten. The original manifest
# is kept and each resume is appended to it.
manifest="$ROOT/run_manifest.txt"
git_commit=$(git rev-parse HEAD)
if [ -n "$(git status --porcelain)" ]; then git_worktree=dirty; else git_worktree=clean; fi
manifest_value() {  # manifest_value <key> <default if absent>
  local v; v=$(grep -m1 "^$1=" "$manifest" | cut -d= -f2-); echo "${v:-$2}"
}
resume=0
prior_outputs=""
if [ -s "$manifest" ]; then
  resume=1
  mismatch=""
  [ "$(manifest_value ssp '')" = "$SSP" ] || mismatch+=" ssp"
  [ "$(manifest_value n_horiuchi '')" = "$N_HORIUCHI" ] || mismatch+=" n_horiuchi"
  [ "$(manifest_value decomp_annual '')" = "$DECOMP_ANNUAL" ] || mismatch+=" decomp_annual"
  [ "$(manifest_value skip_mx_plausibility_check 0)" = "$SKIP_MX_PLAUSIBILITY_CHECK" ] || mismatch+=" skip_mx_plausibility_check"
  [ "$(manifest_value assr_one_replacement none)" = "${ASSR_ONE_REPLACEMENT:-none}" ] || mismatch+=" assr_one_replacement"
  first_commit=$(manifest_value git_commit '')
  if ! git cat-file -e "${first_commit}^{commit}" 2>/dev/null; then
    mismatch+=" git_commit(unknown)"
  elif ! git diff --quiet "$first_commit" HEAD -- pipeline/; then
    mismatch+=" pipeline_code"
  fi
  [ -z "$(git status --porcelain -- pipeline/)" ] || mismatch+=" uncommitted_pipeline_changes"
  if [ -n "$mismatch" ]; then
    echo "ERROR: $ROOT holds outputs from a run that differs in:$mismatch"
    echo "       Resuming would reuse them under false provenance. Use a new"
    echo "       BATCH_ROOT_OVERRIDE, or deliberately clear that root, before rerunning."
    exit 2
  fi
elif [ -n "$(find "$ROOT" "$BATCH_ROOT/geographies/ssp$SSP" -name .done -print -quit 2>/dev/null)" ]; then
  # Outputs written before run manifests existed (e.g. canonical SSP3).
  if [ "${ACCEPT_UNRECORDED_OUTPUTS:-0}" != "1" ]; then
    echo "ERROR: $ROOT holds finished outputs but no run_manifest.txt, so their"
    echo "       provenance cannot be checked. Set ACCEPT_UNRECORDED_OUTPUTS=1 to reuse"
    echo "       them; the manifest will then record prior_outputs=unrecorded_provenance."
    exit 2
  fi
  prior_outputs=unrecorded_provenance
fi

rm -f "$ROOT/failed.txt"   # failures are re-evaluated on every run

# Record the exact resolved run inputs beside the scenario outputs. The tracked
# config declares the requested runs; these files preserve the expanded city
# and GCM domains and the code revision actually used.
if [ "$resume" = 1 ] && [ -s "$ROOT/run_cities.txt" ]; then
  printf '%s\n' $CITY_LIST | sort -u - "$ROOT/run_cities.txt" -o "$ROOT/run_cities.txt"
else
  printf '%s\n' $CITY_LIST > "$ROOT/run_cities.txt"
fi
printf '%s\n' $GCMS > "$ROOT/run_gcms.txt"
if [ "$resume" = 1 ]; then
  {
    echo "resumed_at=$(date -Iseconds)"
    echo "resume_git_commit=$git_commit"
    echo "resume_git_worktree=$git_worktree"
    echo "resume_run_config=${RUN_CONFIG_SOURCE:-command_line}"
    echo "resume_cities=$CITIES"
  } >> "$manifest"
else
  {
    echo "started_at=$(date -Iseconds)"
    echo "git_commit=$git_commit"
    echo "git_worktree=$git_worktree"
    echo "run_config=${RUN_CONFIG_SOURCE:-command_line}"
    echo "ssp=$SSP"
    echo "cities=$CITIES"
    echo "ncores=$NCORES"
    echo "n_horiuchi=$N_HORIUCHI"
    echo "decomp_annual=$DECOMP_ANNUAL"
    echo "batch_root=$BATCH_ROOT"
    echo "skip_mx_plausibility_check=$SKIP_MX_PLAUSIBILITY_CHECK"
    echo "assr_one_replacement=${ASSR_ONE_REPLACEMENT:-none}"
    echo "n_cities=$(printf '%s\n' $CITY_LIST | wc -l)"
    echo "n_gcms=$(printf '%s\n' $GCMS | wc -l)"
    [ -z "$prior_outputs" ] || echo "prior_outputs=$prior_outputs"
  } > "$manifest"
fi

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
printf '%s\n' $CITY_LIST | xargs -P "$NCORES" -I{} bash -c 'dem_job "$1"' _ {}
echo "$(date '+%F %T') steps 01-04: city x GCM"
for c in $CITY_LIST; do for g in $GCMS; do echo "$c $g"; done; done | xargs -P "$NCORES" -L1 bash -c 'gcm_job "$0" "$1"'
echo "$(date '+%F %T') steps 01-04: city ensemble"
printf '%s\n' $CITY_LIST | xargs -P "$NCORES" -I{} bash -c 'ens_job "$1"' _ {}

if [ -s "$ROOT/failed.txt" ]; then
  nfail=$(wc -l < "$ROOT/failed.txt")
  echo "ERROR: $nfail city pipeline jobs failed; see $ROOT/failed.txt"
  echo "Canonical collected objects and pooled geographies were not updated."
  exit 1
fi

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
  exit 1
fi
