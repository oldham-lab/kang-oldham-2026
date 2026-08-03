#!/bin/bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$DIR/logs"
mkdir -p "$LOG_DIR"

SCRIPTS=(
   # "calc_projIndexSE_SEAAD2024_DFC_optimized_loop.py" # done
    "calc_projIndexSE_SEAAD2024_MTG_optimized_loop.py"
    "calc_projIndexSE_SEAAD2024_DFC_ROSMAP_optimized_loop.py"
    "calc_projIndexSE_SEAAD2024_MTG_ROSMAP_optimized_loop.py"
    "calc_projIndexSE_SEAAD2024_DFC_APOE_optimized_loop.py"
)

echo "Starting run: $(date)"
echo "Logs → $LOG_DIR"
echo ""

for script in "${SCRIPTS[@]}"; do
    echo "──────────────────────────────────────────"
    echo "Running: $script"
    echo "Start:   $(date)"
    log="$LOG_DIR/${script%.py}.log"
    python "$DIR/$script" 2>&1 | tee "$log"
    echo "Done:    $(date)"
    echo ""
done

echo "All scripts complete: $(date)"
