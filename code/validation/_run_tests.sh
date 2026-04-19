#!/usr/bin/env bash
# Sequentially run the three validation test scripts, logging output.
set -euo pipefail

LOG_DIR="output/validation/_test_logs"
mkdir -p "$LOG_DIR"

run_one () {
  local script="$1"
  local name="$2"
  echo "[$(date '+%H:%M:%S')] starting $name"
  time Rscript "$script" > "$LOG_DIR/${name}.log" 2>&1
  echo "[$(date '+%H:%M:%S')] done $name -> $LOG_DIR/${name}.log"
}

run_one code/validation/_test_mundlak.R          mundlak
run_one code/validation/_test_lagged_vote.R      lagged_vote
run_one code/validation/_test_year_window_within.R year_window

echo "[$(date '+%H:%M:%S')] all done"
