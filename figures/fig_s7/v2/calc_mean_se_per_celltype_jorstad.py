"""
# Counts are already log-transformed in adata.X

calc_mean_se_per_celltype_jorstad.py

Calculates the mean and standard error (SE) of gene expression per cell type
for the two Jorstad 2023 AnnData objects, restricted to a specified subset of genes.

For each (adata, cell type) combination:
  1. Optionally subset cells to those matching the normal/control condition
  2. Subset the expression matrix to the genes of interest AND cells of that type
  3. Flatten (ravel) the resulting submatrix into a 1-D vector
  4. Compute mean and SE (Bessel-corrected: SE = SD / sqrt(n)) over all elements

Edit the CONFIG section below before running.

Output
------
  <adata_basename>_mean_se.csv  -- one file per input adata, with columns:
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

# Paths to the two .h5ad files
ADATA_PATHS = [
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.h5ad",
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_MTG.h5ad",
]

# Path to gene list file (one gene per line), or a Python list of gene names
GENES = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/genes.txt")
# GENES = ["GENE1", "GENE2", "GENE3", ...]  # alternative: inline list

# Directory to write output CSV files
OUTPUT_DIR = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v2")

# Per-dataset configuration — one entry per adata (must match order of ADATA_PATHS)
# Fields:
#   celltype_col  : adata.obs column containing cell type labels
#   filter        : True to subset to normal samples only; False to use all cells
#   condition_col : adata.obs column used to identify normal samples (ignored if filter=False)
#   normal_value  : value in condition_col that means normal/control (ignored if filter=False)
DATASET_CONFIG = [
    {
        "celltype_col":  "Cell_Type",
        "filter":        False,
        "condition_col": "Disorder",
        "normal_value":  "control",
    },
    {
        "celltype_col":  "Cell_Type",
        "filter":        False,
        "condition_col": "Disorder",
        "normal_value":  "control",
    },
]

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


def process_adata(adata, gene_subset: list, celltype_col: str,
                  condition_col: str, normal_value: str) -> pd.DataFrame:
    """
    Optionally filter to normal samples, then compute mean and SE per cell type.

    Parameters
    ----------
    adata         : AnnData object (cells x genes)
    gene_subset   : list of gene names to restrict to
    celltype_col  : adata.obs column with cell type labels
    condition_col : adata.obs column for filtering (None = no filtering)
    normal_value  : value identifying normal samples (None = no filtering)

    Returns
    -------
    DataFrame: CellType | Mean | SE | N_observations | N_cells | N_genes
    """
    # -- Optionally filter to normal samples ----------------------------------
    if condition_col is not None and normal_value is not None:
        if condition_col not in adata.obs.columns:
            log.error("  condition_col '%s' not found in adata.obs. Available: %s",
                      condition_col, list(adata.obs.columns))
            sys.exit(1)
        before = adata.n_obs
        adata  = adata[adata.obs[condition_col].astype(str) == str(normal_value)].copy()
        log.info("  Filtered to normal ('%s' == '%s'): %d -> %d cells.",
                 condition_col, normal_value, before, adata.n_obs)
        if adata.n_obs == 0:
            log.error("  No cells remain after filtering. "
                      "Check condition_col and normal_value in DATASET_CONFIG.")
            return pd.DataFrame(columns=["CellType", "Mean", "SE",
                                         "N_observations", "N_cells", "N_genes"])
    else:
        log.info("  No condition filter applied — using all %d cells.", adata.n_obs)

    # -- Choose expression matrix ---------------------------------------------
    if USE_RAW and adata.raw is not None:
        X   = adata.raw.X
        var = adata.raw.var
        log.info("  Using adata.raw.X  (%d cells x %d genes)", X.shape[0], X.shape[1])
    else:
        if USE_RAW:
            log.warning("  USE_RAW=True but adata.raw is None; falling back to adata.X")
        X   = adata.X
        var = adata.var
        log.info("  Using adata.X  (%d cells x %d genes)", X.shape[0], X.shape[1])

    # -- Intersect requested genes with genes in this adata -------------------
    adata_genes   = list(var.index)
    gene_set      = set(adata_genes)
    present_genes = [g for g in gene_subset if g in gene_set]
    missing       = len(gene_subset) - len(present_genes)
    if missing:
        log.warning("  %d / %d requested genes not found in this adata; skipped.",
                    missing, len(gene_subset))
    if not present_genes:
        log.error("  No requested genes found. Returning empty table.")
        return pd.DataFrame(columns=["CellType", "Mean", "SE",
                                     "N_observations", "N_cells", "N_genes"])
    gene_idx = [adata_genes.index(g) for g in present_genes]
    log.info("  %d genes retained after intersection.", len(present_genes))

    # -- Check celltype column ------------------------------------------------
    if celltype_col not in adata.obs.columns:
        log.error("  celltype_col '%s' not found in adata.obs. Available: %s",
                  celltype_col, list(adata.obs.columns))
        sys.exit(1)

    celltypes  = adata.obs[celltype_col].astype(str)
    unique_cts = sorted(celltypes.unique())
    log.info("  %d cell types found in '%s'.", len(unique_cts), celltype_col)

    # -- Subset to gene columns first, then iterate over cell types -----------
    X_genes = X[:, gene_idx]

    rows = []
    for ct in unique_cts:
        cell_mask = (celltypes == ct).values
        n_cells   = int(cell_mask.sum())
        sub       = X_genes[cell_mask, :]
        if scipy.sparse.issparse(sub):
            sub = sub.toarray()
        mean, se, n_obs = compute_stats_for_matrix(sub)
        rows.append((ct, mean, se, n_obs, n_cells, len(present_genes)))
        log.debug("    %-25s  n_cells=%5d  n_obs=%10d  mean=%.4g  se=%.4g",
                  ct, n_cells, n_obs, mean, se)

    return pd.DataFrame(rows, columns=["CellType", "Mean", "SE",
                                       "N_observations", "N_cells", "N_genes"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # -- Validate config ------------------------------------------------------
    if len(ADATA_PATHS) != 2:
        log.error("ADATA_PATHS must contain exactly 2 paths (got %d).", len(ADATA_PATHS))
        sys.exit(1)
    if len(DATASET_CONFIG) != 2:
        log.error("DATASET_CONFIG must contain exactly 2 entries (got %d).", len(DATASET_CONFIG))
        sys.exit(1)

    # -- Load gene list --------------------------------------------------------
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

    # -- Load anndata ---------------------------------------------------------
    try:
        import anndata
    except ImportError:
        log.error("anndata is not installed. Run: pip install anndata")
        sys.exit(1)

    out_dir = Path(OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    for i, (h5ad_path, cfg) in enumerate(zip(ADATA_PATHS, DATASET_CONFIG)):
        h5ad_path = Path(h5ad_path)
        log.info("=" * 60)
        log.info("[%d/2] %s  |  filter=%s", i + 1, h5ad_path.name, cfg["filter"])

        if not h5ad_path.exists():
            log.error("  File not found: %s -- skipping.", h5ad_path)
            continue

        try:
            adata = anndata.read_h5ad(h5ad_path)
        except Exception as exc:
            log.error("  Failed to load %s: %s -- skipping.", h5ad_path, exc)
            continue

        log.info("  Loaded: %d cells x %d genes", adata.n_obs, adata.n_vars)

        result = process_adata(
            adata         = adata,
            gene_subset   = gene_subset,
            celltype_col  = cfg["celltype_col"],
            condition_col = cfg["condition_col"] if cfg["filter"] else None,
            normal_value  = cfg["normal_value"]  if cfg["filter"] else None,
        )

        out_path = out_dir / f"{h5ad_path.stem}_mean_se.csv"
        result.to_csv(out_path, index=False, float_format="%.6g")
        log.info("  Wrote %s  (%d cell types)", out_path, len(result))

    log.info("=" * 60)
    log.info("Done.")


if __name__ == "__main__":
    main()
