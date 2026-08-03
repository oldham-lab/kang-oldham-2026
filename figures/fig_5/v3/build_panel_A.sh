#!/usr/bin/env bash
# Build fig 5 panel A end-to-end from scratch.
#   1. build_panelA_cells.R  -> panelA_cells.csv   (data + colours, seed 23)
#   2. assemble_panel_A.py   -> panel_A_generated.svg (layout, labels, arrows)
# Output is written as panel_A_generated.svg so the hand-built panel_A.svg is
# left untouched; rename it over panel_A.svg once you're happy with it.
set -euo pipefail
cd "$(dirname "$0")"

echo "[1/2] generating cell colours (R)…"
Rscript build_panelA_cells.R

echo "[2/2] assembling SVG (python)…"
python3 assemble_panel_A.py panel_A_generated.svg

echo "done -> $(pwd)/panel_A_generated.svg"
