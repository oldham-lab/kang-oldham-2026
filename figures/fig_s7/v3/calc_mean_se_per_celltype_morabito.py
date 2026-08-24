"""
calc_mean_se_per_celltype_morabito.py

Calculates the mean and standard error (SE) of gene expression per cell type
for the Morabito 2021 dataset, restricted to a specified subset of genes and
to control (non-AD) samples.

Data format: Matrix Market (.mtx) + barcodes.tsv + genes.tsv + metadata CSV.

For each cell type:
  1. Subset cells to the normal/control condition
  2. Subset the expression matrix to the genes of interest AND cells of that type
  3. Flatten (ravel) the resulting submatrix into a 1-D vector
  4. Compute mean and SE (Bessel-corrected: SE = SD / sqrt(n)) over all elements

Edit the CONFIG section below before running.

Output
------
  morabito_2021_mean_se.csv  -- columns:
    CellType | Mean | SE | N_observations | N_cells | N_genes
"""

import os
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hgnc_common import select_common_indices

# =============================================================================
# CONFIG — edit these variables before running
# =============================================================================

DATA_DIR = (
    os.path.join(os.environ.get("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/"
    "morabito_2021/all_cells")
)

MATRIX_FILE   = f"{DATA_DIR}/matrix.mtx"
BARCODES_FILE = f"{DATA_DIR}/barcodes.tsv"
GENES_FILE    = f"{DATA_DIR}/genes.tsv"
METADATA_FILE = f"{DATA_DIR}/metadata_withABIanno_filtered.csv"

# Path to gene list file (one gene per line), or a Python list of gene names.
# v3: the common gene universe (intersection across all datasets, canonical HGNC).
GENES = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v3/common_genes.txt")
# GENES = ["GENE1", "GENE2", "GENE3", ...]  # alternative: inline list

# Directory to write output CSV
OUTPUT_DIR = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v3")

# metadata column containing cell type labels
CELLTYPE_COL = "Subclass.x"

# metadata column and value identifying normal/control samples
CONDITION_COL = "Diagnosis"
NORMAL_VALUE  = "Control"

# =============================================================================
# END CONFIG
# =============================================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


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

    # -- Load matrix files ------------------------------------------------------
    log.info("Loading matrix from %s", MATRIX_FILE)
    # MTX is genes x cells; transpose to cells x genes (CSR for row slicing)
    mat = scipy.io.mmread(MATRIX_FILE).T.tocsr()
    mat = mat.astype(np.float64)
    mat.data = np.log1p(mat.data)
    log.info("Applied log1p transformation to matrix.")
    log.info("Matrix shape (cells x genes): %d x %d", mat.shape[0], mat.shape[1])

    barcodes = pd.read_csv(BARCODES_FILE, header=None, names=["Barcode"])["Barcode"].tolist()
    genes    = pd.read_csv(GENES_FILE,    header=None, names=["Gene"])["Gene"].tolist()
    log.info("Barcodes: %d  |  Genes: %d", len(barcodes), len(genes))

    if mat.shape[0] != len(barcodes) or mat.shape[1] != len(genes):
        log.error(
            "Dimension mismatch: matrix %s, barcodes %d, genes %d",
            mat.shape, len(barcodes), len(genes),
        )
        sys.exit(1)

    # -- Load metadata ----------------------------------------------------------
    log.info("Loading metadata from %s", METADATA_FILE)
    meta = pd.read_csv(METADATA_FILE)
    log.info("Metadata: %d rows x %d columns", meta.shape[0], meta.shape[1])

    # Filter matrix to only barcodes present in metadata
    if "Barcode" not in meta.columns:
        log.error("'Barcode' column not found in metadata.")
        sys.exit(1)
    meta = meta.set_index("Barcode")
    meta_barcodes = set(meta.index)
    keep_idx  = [i for i, b in enumerate(barcodes) if b in meta_barcodes]
    discarded = len(barcodes) - len(keep_idx)
    if discarded:
        log.info("Discarding %d / %d barcodes not in metadata; keeping %d.",
                 discarded, len(barcodes), len(keep_idx))
    mat      = mat[keep_idx, :]
    barcodes = [barcodes[i] for i in keep_idx]
    meta     = meta.reindex(barcodes)  # align row order to filtered matrix

    # -- Filter to control samples ---------------------------------------------
    if CONDITION_COL not in meta.columns:
        log.error("CONDITION_COL '%s' not found in metadata. Available: %s",
                  CONDITION_COL, list(meta.columns))
        sys.exit(1)

    ctrl_mask = meta[CONDITION_COL].astype(str) == str(NORMAL_VALUE)
    log.info("Filtered to control ('%s' == '%s'): %d / %d cells.",
             CONDITION_COL, NORMAL_VALUE, ctrl_mask.sum(), len(ctrl_mask))
    if ctrl_mask.sum() == 0:
        log.error("No cells remain after filtering.")
        sys.exit(1)

    mat_ctrl  = mat[ctrl_mask.values, :]
    meta_ctrl = meta[ctrl_mask.values]

    # -- Select common-universe genes (canonical HGNC match) ------------------
    gene_idx, present_genes = select_common_indices(genes, gene_subset)
    missing       = len(gene_subset) - len(present_genes)
    if missing:
        log.warning("%d / %d common-universe genes not found; skipped.",
                    missing, len(gene_subset))
    if not present_genes:
        log.error("No requested genes found in dataset.")
        sys.exit(1)
    log.info("%d genes retained (canonical common-universe match).", len(present_genes))

    mat_genes = mat_ctrl[:, gene_idx]

    # -- Compute stats per cell type -------------------------------------------
    if CELLTYPE_COL not in meta_ctrl.columns:
        log.error("CELLTYPE_COL '%s' not found in metadata. Available: %s",
                  CELLTYPE_COL, list(meta_ctrl.columns))
        sys.exit(1)

    celltypes  = meta_ctrl[CELLTYPE_COL].astype(str)
    unique_cts = sorted(celltypes.unique())
    log.info("%d cell types found in '%s'.", len(unique_cts), CELLTYPE_COL)

    rows = []
    for ct in unique_cts:
        cell_mask = (celltypes == ct).values
        n_cells   = int(cell_mask.sum())
        sub       = mat_genes[cell_mask, :]
        if scipy.sparse.issparse(sub):
            sub = sub.toarray()
        mean, se, n_obs = compute_stats_for_matrix(sub)
        rows.append((ct, mean, se, n_obs, n_cells, len(present_genes)))
        log.debug("  %-25s  n_cells=%5d  n_obs=%10d  mean=%.4g  se=%.4g",
                  ct, n_cells, n_obs, mean, se)

    result = pd.DataFrame(rows, columns=["CellType", "Mean", "SE",
                                         "N_observations", "N_cells", "N_genes"])

    # -- Write output ----------------------------------------------------------
    out_dir  = Path(OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "morabito_2021_mean_se.csv"
    result.to_csv(out_path, index=False, float_format="%.6g")
    log.info("Wrote %s  (%d cell types)", out_path, len(result))


if __name__ == "__main__":
    main()
