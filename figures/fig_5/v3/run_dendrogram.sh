#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Input files ---
MOD_EIG="${1:-${REPO_DIR:-/home/gugene/code/git/kang-oldham-2026}/figures/fig_5/v3/mod_eig.csv}"
BP_GENES="${2:-${REPO_DIR:-/home/gugene/code/git/kang-oldham-2026}/figures/fig_5/v2/branchpoint_table_modeig_with_genes.csv}"

# --- Optional args ---
OUTPUT="${3:-dendrogram_only_test.pdf}"
CUT_HEIGHT="${4:-0.3}"
FONT_SIZE="${5:-14}"

python "$SCRIPT_DIR/make_labeled_dendrogram_only_v1.2.py" \
    "$MOD_EIG" \
    "$BP_GENES" \
    -o "$OUTPUT" \
    --cut-height "$CUT_HEIGHT" \
    --font-size "$FONT_SIZE"
