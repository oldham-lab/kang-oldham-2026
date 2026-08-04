"""
Metacell mapping: Gabitto (DFC + MTG) → Liu (PFC + MTC)

Processes DFC and MTG sequentially.
Runs with log-normalised counts and saves a CSV per region.

Outputs:
  liu_gabitto_metacell_labels_DFC_lognorm.csv
  liu_gabitto_metacell_labels_MTG_lognorm.csv
"""

import os
import logging
import sys
import numpy as np
import anndata as ad
import scipy.sparse as sp
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
GABITTO_PATH_DFC = os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/"
                                                                  "SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad")

GABITTO_PATH_MTG = (
    os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
)

LIU_PATH = (
    os.path.join(
        os.environ.get("SHARED_DATA_DIR", "/mnt/bdata/@shared"),
        "scsn.expr_data/human_expr/postnatal/"
        "MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025.h5ad",
    )
)
SAVE_DIR = Path(os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s6"))

# ── Regions to process: (Gabitto region name, Gabitto path, Liu region name) ─
# PFC is treated as equivalent to DFC; MTC is treated as equivalent to MTG.
REGIONS = [
    ("DFC", GABITTO_PATH_DFC, "PFC"),
    ("MTG", GABITTO_PATH_MTG, "MTC"),
]

# ── Column names ──────────────────────────────────────────────────────────────
GAB_SUBCLASS_COL = "Subclass"
LIU_SUBCLASS_COL = "RNA.Subclass"
LIU_REGION_COL   = "BrainRegion"

# ── Filters ───────────────────────────────────────────────────────────────────
GAB_EXCLUDE = {"VLMC"}
LIU_EXCLUDE = {"Exc HC", "Exc EC", "T", "SMC", "Fib", "Per"}

LOG_PATH = SAVE_DIR / "metacell_mapping_DFC_MTG.log"

# ── Logging ───────────────────────────────────────────────────────────────────
SAVE_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(LOG_PATH, mode="w"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

# ── Settings ──────────────────────────────────────────────────────────────────
TARGET_SUM = 1e4      # counts per cell for normalisation (used in lognorm only)
CHUNK_SIZE = 50_000   # Liu cells per chunk for correlation

# =============================================================================
# Helpers
# =============================================================================

def to_dense(mat):
    """Convert sparse or dense matrix to dense float32 array."""
    if sp.issparse(mat):
        mat = mat.toarray()
    return mat.astype(np.float32)


def norm_only(mat):
    """Normalise to TARGET_SUM per cell, no log transform → dense float32."""
    mat = to_dense(mat)
    row_sums = mat.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    return mat / row_sums * TARGET_SUM


def lognorm(mat):
    """Normalise to TARGET_SUM per cell then log1p → dense float32."""
    mat = norm_only(mat)
    np.log1p(mat, out=mat)
    return mat


def pearson_rows_vs_matrix(rows, ref):
    """
    Pearson correlation between every row of `rows` (n × g) and every row
    of `ref` (m × g).  Returns (n × m) float32 array.
    """
    rows_c = rows - rows.mean(axis=1, keepdims=True)
    ref_c  = ref  - ref.mean(axis=1, keepdims=True)

    rows_norm = np.linalg.norm(rows_c, axis=1, keepdims=True)
    ref_norm  = np.linalg.norm(ref_c,  axis=1, keepdims=True)

    rows_norm[rows_norm == 0] = 1
    ref_norm[ref_norm   == 0] = 1

    rows_c /= rows_norm
    ref_c  /= ref_norm

    return rows_c @ ref_c.T   # (n × m)


# =============================================================================
# Load Liu once (full) before the region loop
# =============================================================================
log.info("=" * 60)
log.info("Loading Liu (full dataset) ...")
liu_full = ad.read_h5ad(LIU_PATH)
log.info(f"  {liu_full.n_obs:,} nuclei, {liu_full.n_vars:,} genes")

if LIU_REGION_COL not in liu_full.obs.columns:
    raise KeyError(
        f"Region column '{LIU_REGION_COL}' not found. "
        f"Available: {list(liu_full.obs.columns)}"
    )

# =============================================================================
# Main loop: process each region sequentially
# =============================================================================
for gab_region, gabitto_path, liu_region in REGIONS:
    log.info("\n" + "=" * 60)
    log.info(f"Region: {gab_region} (Gabitto) ↔ {liu_region} (Liu)")
    log.info("=" * 60)

    # ── Step 1: Load Gabitto, build metacells for both transforms, free memory ─
    log.info(f"Loading Gabitto {gab_region} (backed/memory-mapped) ...")
    gab = ad.read_h5ad(gabitto_path, backed="r")
    log.info(f"  {gab.n_obs:,} nuclei, {gab.n_vars:,} genes")

    keep_mask = ~gab.obs[GAB_SUBCLASS_COL].isin(GAB_EXCLUDE)
    gab_obs   = gab.obs[keep_mask]
    log.info(f"  {keep_mask.sum():,} nuclei after filtering {GAB_EXCLUDE}")

    subclasses = sorted(gab_obs[GAB_SUBCLASS_COL].unique())
    log.info(f"  {len(subclasses)} subclasses: {subclasses}")

    gabitto_genes = np.array(gab.var_names)

    metacells_both = {}
    # for transform_name, transform_fn in [("raw", norm_only), ("lognorm", lognorm)]:
    for transform_name, transform_fn in [("lognorm", lognorm)]:
        log.info(f"\nBuilding Gabitto metacells ({transform_name}) ...")
        mc = np.zeros((len(subclasses), gab.n_vars), dtype=np.float32)
        for i, sc in enumerate(subclasses):
            mask = (gab.obs[GAB_SUBCLASS_COL] == sc) & keep_mask
            idx = np.where(mask.values)[0]
            subclass_expr = transform_fn(gab.layers["UMIs"][idx])
            mc[i] = subclass_expr.mean(axis=0)
            log.info(f"  {sc}: {mask.sum():,} nuclei")
            del subclass_expr
        metacells_both[transform_name] = mc

    gab.file.close()
    del gab
    log.info(f"\nGabitto {gab_region} freed from memory.")

    # ── Step 2: Filter Liu to region, align genes ─────────────────────────────
    log.info(f"\nFiltering Liu to {liu_region} ...")
    liu = liu_full[liu_full.obs[LIU_REGION_COL] == liu_region].copy()
    log.info(f"  {liu.n_obs:,} nuclei after filtering to {liu_region}")

    liu = liu[~liu.obs[LIU_SUBCLASS_COL].isin(LIU_EXCLUDE)].copy()
    log.info(f"  {liu.n_obs:,} nuclei after excluding {LIU_EXCLUDE}")

    common_genes = np.intersect1d(gabitto_genes, liu.var_names)
    log.info(f"  Common genes: {len(common_genes):,}")

    gab_gene_map = {g: i for i, g in enumerate(gabitto_genes)}
    gab_gene_idx = np.array([gab_gene_map[g] for g in common_genes])
    liu_sub = liu[:, common_genes].copy()
    n_liu = liu_sub.n_obs
    del liu

    # Subset metacells to common genes
    for name in metacells_both:
        metacells_both[name] = metacells_both[name][:, gab_gene_idx]

    # ── Step 3: Run correlation for each transform ────────────────────────────
    # for transform_name, transform_fn in [("raw", norm_only), ("lognorm", lognorm)]:
    for transform_name, transform_fn in [("lognorm", lognorm)]:
        log.info("\n" + "=" * 60)
        log.info(f"Running: {gab_region} / {transform_name}")
        log.info("=" * 60)

        metacells = metacells_both[transform_name]

        log.info(f"Assigning labels in chunks of {CHUNK_SIZE:,} ...")
        best_labels = np.empty(n_liu, dtype=object)
        best_corrs  = np.full(n_liu, -2.0, dtype=np.float32)

        for start in range(0, n_liu, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_liu)
            log.info(f"  Chunk {start:,}–{end:,} ...")
            chunk_expr = transform_fn(liu_sub.X[start:end])
            corr_mat   = pearson_rows_vs_matrix(chunk_expr, metacells)
            best_idx   = corr_mat.argmax(axis=1)
            best_labels[start:end] = [subclasses[i] for i in best_idx]
            best_corrs[start:end]  = corr_mat[np.arange(len(best_idx)), best_idx]
            del chunk_expr, corr_mat

        out_obs = liu_sub.obs.copy()
        out_obs["Gabitto_metacell_labels"]   = best_labels
        out_obs["Gabitto_metacell_max_corr"] = best_corrs

        log.info("\nLabel distribution:")
        log.info(out_obs["Gabitto_metacell_labels"].value_counts().to_string())

        out_path = SAVE_DIR / f"liu_gabitto_metacell_labels_{gab_region}_{transform_name}.csv"
        out_obs.to_csv(out_path)
        log.info(f"\nSaved: {out_path}  ({n_liu:,} rows)")

        del metacells, best_labels, best_corrs

    del liu_sub, metacells_both, gabitto_genes

del liu_full
log.info("\nDone.")
