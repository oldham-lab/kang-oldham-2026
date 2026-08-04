"""
adata.X counts are raw

calc_mean_se_per_celltype_MIT_AD.py

Calculates the mean and standard error (SE) of gene expression per cell type
for the MIT AD Multiomic Multiregion dataset, restricted to a specified subset
of genes and to non-AD (normal) samples.

Per-region cell-to-label mappings are supplied via external CSV files
(one per region) that contain Gabitto metacell labels and Pathology info.
This script produces one output CSV for each of the two target regions
(PFC and MTC).

For each (region, cell type) combination:
  1. Load the region-specific label CSV to get cell barcodes and Gabitto labels
  2. Subset the h5ad to those barcodes
  3. Subset cells to the normal/control condition (nonAD)
  4. Subset the expression matrix to the genes of interest AND cells of that type
  5. Flatten (ravel) the resulting submatrix into a 1-D vector
  6. Compute mean and SE (Bessel-corrected: SE = SD / sqrt(n)) over all elements

Edit the CONFIG section below before running.

Output
------
  MIT_AD_PFC_mean_se.csv  -- PFC region
  MIT_AD_MTC_mean_se.csv  -- MTC region

  Each file has columns:
    CellType | Mean | SE | N_observations | N_cells | N_genes
"""

import os
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse

# =============================================================================
# CONFIG — edit these variables before running
# =============================================================================

# Path to the MIT AD Multiomic Multiregion h5ad file
ADATA_PATH = (
    os.path.join(
        os.environ.get("SHARED_DATA_DIR", "/mnt/bdata/@shared"),
        "scsn.expr_data/human_expr/postnatal/"
        "MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025.h5ad",
    )
)

# Per-region label CSVs — first column is the cell barcode (index);
# must contain columns: Gabitto_metacell_labels, Pathology
LABEL_FILES = {
    "PFC": os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s6/liu_gabitto_metacell_labels_DFC_lognorm.csv"),
    "MTC": os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s6/liu_gabitto_metacell_labels_MTG_lognorm.csv"),
}

# Path to gene list file (one gene per line), or a Python list of gene names
GENES = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/genes.txt")
# GENES = ["GENE1", "GENE2", "GENE3", ...]  # alternative: inline list

# Directory to write output CSV files
OUTPUT_DIR = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v2")

# Column in label CSVs containing cell type labels
CELLTYPE_COL = "Gabitto_metacell_labels"

# Column and value in label CSVs used to identify normal/control samples
CONDITION_COL = "Pathology"
NORMAL_VALUE  = "nonAD"

# Set to True to use adata.raw.X instead of adata.X (when available)
USE_RAW = False

# =============================================================================
# END CONFIG
# =============================================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_stats_for_matrix(mat: np.ndarray) -> tuple:
    """
    Ravel a 2-D (cells x genes) array to 1-D and compute mean and
    Bessel-corrected SE over all elements.

    Returns (mean, se, n).
    """
    vec = mat.ravel().astype(np.float64)
    n   = len(vec)
    if n < 2:
        return np.nan, np.nan, n
    mean = vec.mean()
    sd   = vec.std(ddof=1)
    se   = sd / np.sqrt(n)
    return mean, se, n


def process_region(adata, region: str, label_file: str, gene_subset: list) -> pd.DataFrame:
    """
    Load the region label CSV, subset the adata to matching barcodes,
    filter to normal samples, then compute mean and SE per Gabitto cell type.

    Parameters
    ----------
    adata       : AnnData object (all cells x genes)
    region      : region name string (for logging)
    label_file  : path to per-region label CSV
    gene_subset : list of gene names to restrict to

    Returns
    -------
    DataFrame: CellType | Mean | SE | N_observations | N_cells | N_genes
    """
    empty = pd.DataFrame(columns=["CellType", "Mean", "SE",
                                   "N_observations", "N_cells", "N_genes"])

    # -- Load label CSV ---------------------------------------------------------
    label_path = Path(label_file)
    if not label_path.exists():
        log.error("  Label file not found: %s", label_path)
        return empty

    labels = pd.read_csv(label_path, index_col=0)
    log.info("  Loaded label CSV: %d cells, columns: %s",
             len(labels), list(labels.columns))

    for col in [CELLTYPE_COL, CONDITION_COL]:
        if col not in labels.columns:
            log.error("  Column '%s' not found in label CSV. Available: %s",
                      col, list(labels.columns))
            return empty

    # -- Filter label CSV to normal samples ------------------------------------
    before = len(labels)
    labels = labels[labels[CONDITION_COL].astype(str) == str(NORMAL_VALUE)]
    log.info("  Filtered to normal ('%s' == '%s'): %d -> %d cells.",
             CONDITION_COL, NORMAL_VALUE, before, len(labels))
    if len(labels) == 0:
        log.error("  No cells remain after filtering.")
        return empty

    # -- Subset adata to barcodes in label CSV ----------------------------------
    common = labels.index.intersection(adata.obs_names)
    missing_bc = len(labels) - len(common)
    if missing_bc:
        log.warning("  %d barcodes in label CSV not found in adata; skipped.",
                    missing_bc)
    if len(common) == 0:
        log.error("  No matching barcodes between label CSV and adata.")
        return empty

    adata_sub = adata[common].copy()
    labels_sub = labels.loc[common]
    log.info("  Subsetted adata to %d cells.", adata_sub.n_obs)

    # -- Choose expression matrix -----------------------------------------------
    if USE_RAW and adata_sub.raw is not None:
        X   = adata_sub.raw.X
        var = adata_sub.raw.var
        log.info("  Using adata.raw.X  (%d cells x %d genes)", X.shape[0], X.shape[1])
    else:
        if USE_RAW:
            log.warning("  USE_RAW=True but adata.raw is None; falling back to adata.X")
        X   = adata_sub.X
        var = adata_sub.var
        log.info("  Using adata.X  (%d cells x %d genes)", X.shape[0], X.shape[1])

    # -- Intersect requested genes with genes in adata -------------------------
    adata_genes   = list(var.index)
    present_genes = [g for g in gene_subset if g in set(adata_genes)]
    missing       = len(gene_subset) - len(present_genes)
    if missing:
        log.warning("  %d / %d requested genes not found in adata; skipped.",
                    missing, len(gene_subset))
    if not present_genes:
        log.error("  No requested genes found. Returning empty table.")
        return empty
    gene_idx = [adata_genes.index(g) for g in present_genes]
    log.info("  %d genes retained after intersection.", len(present_genes))

    # -- Subset to gene columns first, then iterate over cell types ------------
    X_genes   = X[:, gene_idx]
    celltypes  = labels_sub[CELLTYPE_COL].astype(str)
    unique_cts = sorted(celltypes.unique())
    log.info("  %d cell types found in '%s'.", len(unique_cts), CELLTYPE_COL)

    rows = []
    for ct in unique_cts:
        cell_mask = (celltypes == ct).values
        n_cells   = int(cell_mask.sum())
        sub       = X_genes[cell_mask, :]
        if scipy.sparse.issparse(sub):
            sub = sub.toarray()
        mean, se, n_obs = compute_stats_for_matrix(sub)
        rows.append((ct, mean, se, n_obs, n_cells, len(present_genes)))
        log.debug("    %-30s  n_cells=%5d  n_obs=%10d  mean=%.4g  se=%.4g",
                  ct, n_cells, n_obs, mean, se)

    return pd.DataFrame(rows, columns=["CellType", "Mean", "SE",
                                       "N_observations", "N_cells", "N_genes"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # -- Load gene list ---------------------------------------------------------
    if isinstance(GENES, list):
        gene_subset = [g.strip() for g in GENES if str(g).strip()]
        log.info("Using inline gene list: %d genes.", len(gene_subset))
    else:
        genes_path = Path(GENES)
        if not genes_path.exists():
            log.error("Gene file not found: %s", genes_path)
            sys.exit(1)
        with open(genes_path) as fh:
            gene_subset = [line.strip() for line in fh if line.strip()]
        log.info("Loaded %d genes from %s", len(gene_subset), genes_path)

    if not gene_subset:
        log.error("Gene list is empty.")
        sys.exit(1)

    # -- Load anndata -----------------------------------------------------------
    try:
        import anndata
    except ImportError:
        log.error("anndata is not installed. Run: pip install anndata")
        sys.exit(1)

    h5ad_path = Path(ADATA_PATH)
    if not h5ad_path.exists():
        log.error("File not found: %s", h5ad_path)
        sys.exit(1)

    log.info("=" * 60)
    log.info("Loading %s", h5ad_path.name)
    try:
        adata = anndata.read_h5ad(h5ad_path)
    except Exception as exc:
        log.error("Failed to load %s: %s", h5ad_path, exc)
        sys.exit(1)
    log.info("Loaded: %d cells x %d genes", adata.n_obs, adata.n_vars)

    import scanpy as sc
    sc.pp.log1p(adata)
    log.info("Applied log1p transformation to adata.X.")

    out_dir = Path(OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    # -- Process each region ----------------------------------------------------
    for region, label_file in LABEL_FILES.items():
        log.info("=" * 60)
        log.info("Processing region: %s  |  labels: %s", region, Path(label_file).name)

        result = process_region(
            adata       = adata,
            region      = region,
            label_file  = label_file,
            gene_subset = gene_subset,
        )

        out_path = out_dir / f"MIT_AD_{region}_mean_se.csv"
        result.to_csv(out_path, index=False, float_format="%.6g")
        log.info("  Wrote %s  (%d cell types)", out_path, len(result))

    log.info("=" * 60)
    log.info("Done.")


if __name__ == "__main__":
    main()
