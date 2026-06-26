#!/bin/bash
# Run full pipeline with interleaved validation, background logging
# Usage: bash run_pipeline.sh

cd "$(dirname "$0")/scripts"
LOG_DIR="../results/production"
mkdir -p "$LOG_DIR" "../results/validation_plots"
LOG="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"

echo "=== Pipeline start: $(date) ===" | tee -a "$LOG"
echo "Log: $LOG" | tee -a "$LOG"
echo "Madrid demo mode (city_subset = ES001C)" | tee -a "$LOG"

Rscript --no-save --no-restore 00_RunAll.R 2>&1 | tee -a "$LOG"

echo "=== Pipeline end: $(date) ===" | tee -a "$LOG"