#!/usr/bin/env bash
################################################################################
#
# Figure runner: Part 05 for every city whose ENSEMBLE stage finished
#
#   Usage (from the repository root, after run_batch.sh):
#     pipeline/run_figures.sh <SSP> [NCORES]
#   In the background:
#     nohup pipeline/run_figures.sh 3 32 > figures_ssp3.log 2>&1 &
#
#   Output: results/europe/ssp<SSP>/<city>/ENSEMBLE/figures/05_fig*.png and
#   05_summary.csv; per-city log in ENSEMBLE/fig_log.txt. Cities that already
#   have all three figures are skipped (FORCE=1 redoes them). Failures are
#   listed in <root>/figures_failed.txt.
#
################################################################################
set -uo pipefail

SSP=${1:?"SSP (1, 2 or 3)"}
NCORES=${2:-4}
export ROOT=results/europe/ssp${SSP} SSP FORCE=${FORCE:-0}
[ -d "$ROOT" ] || { echo "ERROR: $ROOT not found (run from the repository root)"; exit 1; }
if grep -q "04_l[ei]_decomposition.csv" pipeline/05_figures.R; then
  echo "ERROR: 05_figures.R still reads the annual decomposition files, which batch runs do not write"; exit 1
fi
rm -f "$ROOT/figures_failed.txt"

fig_job() {  # fig_job <city>
  local c=$1 d=$ROOT/$1/ENSEMBLE
  [ -f "$d/.done" ] || return 0
  if [ "$FORCE" != "1" ] && [ -f "$d/figures/05_fig1_trajectories.png" ] &&
     [ -f "$d/figures/05_fig2_within_branch_by_age_block.png" ] &&
     [ -f "$d/figures/05_fig3_age_profile_cc_effect.png" ]; then return 0; fi
  mkdir -p "$d/figures"
  env CITY_ID="$c" SSP="$SSP" GCM=ENSEMBLE OUT_DIR="$d" CHECK_DIR="$d/checks" FIG_DIR="$d/figures" \
      DEMOG_DIR="$ROOT/$c/demography" Rscript pipeline/05_figures.R > "$d/fig_log.txt" 2>&1 \
    || echo "$c" >> "$ROOT/figures_failed.txt"
}
export -f fig_job

n=$(ls "$ROOT"/*/ENSEMBLE/.done 2>/dev/null | wc -l)
echo "$(date '+%F %T') figures for $n cities, $NCORES in parallel"
ls "$ROOT" | grep -E '^[A-Z]{2}[0-9]{3}C$' | xargs -P "$NCORES" -I{} bash -c 'fig_job {}'
nfail=$(cat "$ROOT/figures_failed.txt" 2>/dev/null | wc -l)
echo "$(date '+%F %T') done; failures: $nfail"
[ "$nfail" -gt 0 ] && { echo "Failed cities (see <city>/ENSEMBLE/fig_log.txt):"; cat "$ROOT/figures_failed.txt"; }
exit 0
