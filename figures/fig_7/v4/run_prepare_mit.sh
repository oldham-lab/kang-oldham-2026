#!/usr/bin/env bash
# Run prepare_mit_multiome_and_pseudobulk.py for both brain regions (PFC and MTC)
# using external Gabitto metacell labels for pseudobulking.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON=python3

H5AD_INPUT="${SHARED_DATA_DIR:-/mnt/bdata/@shared}/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025.h5ad"   
OUTDIR="${REPO_DIR:-/home/gugene/code/git/kang-oldham-2026}/figures/fig_7/v4/mit_pseudobulk"

LABELS_PFC="${REPO_DIR:-/home/gugene/code/git/kang-oldham-2026}/figures/fig_7/liu_gabitto_metacell_labels_DFC_lognorm.csv"
LABELS_MTC="${REPO_DIR:-/home/gugene/code/git/kang-oldham-2026}/figures/fig_7/liu_gabitto_metacell_labels_MTG_lognorm.csv"

echo "=== PFC ==="
$PYTHON "$SCRIPT_DIR/prepare_mit_multiome_and_pseudobulk.py" \
    --input           "$H5AD_INPUT" \
    --outdir          "$OUTDIR" \
    --region-value    PFC \
    --external-labels "$LABELS_PFC"

echo ""
echo "=== MTC ==="
$PYTHON "$SCRIPT_DIR/prepare_mit_multiome_and_pseudobulk.py" \
    --input           "$H5AD_INPUT" \
    --outdir          "$OUTDIR" \
    --region-value    MTC \
    --external-labels "$LABELS_MTC"

echo ""
echo "All done. Output in: $OUTDIR"
