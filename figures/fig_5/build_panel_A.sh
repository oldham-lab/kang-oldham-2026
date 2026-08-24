#!/usr/bin/env bash
# Build fig 5 panel A end-to-end from scratch.
#   1. build_panelA_cells.R  -> panelA_cells.csv   (data + colours, seed 23)
#   2. assemble_panel_A.py   -> panel_A_generated.svg (layout, labels, arrows)
# Lives at fig_5/ rather than in a version folder: the schematic is illustrative
# (fixed-seed demo matrices), not tied to any figure version, and the current
# version folder ships only the finished panel_A.svg.
# Output is written beside this script as panel_A_generated.svg so the hand-built
# panel stays untouched; copy it over <current version>/panel_A.svg when happy.
set -euo pipefail
cd "$(dirname "$0")"

echo "[1/2] generating cell colours (R)…"
Rscript build_panelA_cells.R

echo "[2/2] assembling SVG (python)…"
python3 assemble_panel_A.py panel_A_generated.svg

echo "done -> $(pwd)/panel_A_generated.svg"
