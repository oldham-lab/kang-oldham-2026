"""
diagnose_seaad_mtg_pseudobulk.py

Spot-checks the SEA-AD MTG pseudobulk CSV against the raw h5ad by:
  1. Checking for donor-ID / group-key collisions (underscore in donor names)
  2. Comparing raw cell-level sums to pseudobulk CSV values for a sample of
     (celltype, donor) pairs — any mismatch indicates a bug in the pipeline
  3. Comparing float32 vs float64 accumulation to quantify precision loss
  4. Reporting raw count statistics (integrality, max value, library sizes)

Usage:
    python diagnose_seaad_mtg_pseudobulk.py

Paths are hardcoded below — edit as needed.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import gzip
import os
import random

# ── Paths ─────────────────────────────────────────────────────────────────────
H5AD_PATH = (
    os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/"
    "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
)
PB_CSV_PATH = (
    os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/"
    "mtg_data_pseudobulk/sea_mtg_pseudobulk.csv.gz")
)
PB_ANNO_PATH = (
    os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/"
    "mtg_data_pseudobulk/sea_mtg_pseudobulk_annotations.csv")
)

# Must match what prepare_seaad_mtg_pseudobulk.py used
RAW_LAYER    = "UMIs"
CELLTYPE_COL = "Subclass"
DONOR_COL    = "Donor ID"
DX_COL       = "Overall AD neuropathological Change"
DX_CONTROL   = "Not AD"
DX_AD_VALUES = ["High", "Intermediate", "Low"]

# Number of (celltype, donor) pairs to spot-check
N_SPOT_CHECK = 5

# ── Helpers ───────────────────────────────────────────────────────────────────

def sep(title=""):
    print("\n" + "=" * 60)
    if title:
        print(title)
        print("=" * 60)


def load_pb_csv(path):
    print(f"Loading pseudobulk CSV: {path} ...")
    df = pd.read_csv(path, index_col=0)
    print(f"  {df.shape[0]:,} genes × {df.shape[1]} samples")
    return df


def load_pb_anno(path):
    print(f"Loading pseudobulk annotations: {path} ...")
    anno = pd.read_csv(path)
    print(f"  {anno.shape[0]} rows, columns: {list(anno.columns)}")
    return anno


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # ── 1. Load pseudobulk outputs ────────────────────────────────────────────
    sep("1. PSEUDOBULK OUTPUT CHECKS")
    pb   = load_pb_csv(PB_CSV_PATH)
    anno = load_pb_anno(PB_ANNO_PATH)

    print(f"\n  Pseudobulk sample columns match annotation rows: "
          f"{pb.shape[1]} vs {len(anno)}")
    if pb.shape[1] != len(anno):
        print("  WARNING: column count mismatch between CSV and annotation file")

    # Group balance
    if "Dx" in anno.columns:
        print("\n  Dx distribution:")
        print(anno["Dx"].value_counts().to_string())
    else:
        print(f"  WARNING: 'Dx' column not found. Columns: {list(anno.columns)}")

    # Celltype distribution
    if CELLTYPE_COL in anno.columns:
        print(f"\n  Celltype distribution ({CELLTYPE_COL}):")
        print(anno[CELLTYPE_COL].value_counts().to_string())

    # ── 2. Donor ID collision check ───────────────────────────────────────────
    sep("2. DONOR ID COLLISION CHECK")
    print("  Checking whether any Donor IDs contain '___' (group key separator) ...")
    if DONOR_COL in anno.columns:
        collision = anno[DONOR_COL].astype(str).str.contains("___")
        if collision.any():
            print(f"  WARNING: {collision.sum()} donor IDs contain '___'")
            print("  Affected:", anno[DONOR_COL][collision].unique().tolist())
        else:
            print("  OK: no donor IDs contain '___'")

        # Also check for single underscore that could cause ambiguous splits
        multi_under = anno[DONOR_COL].astype(str).str.contains("_")
        n_multi = multi_under.sum()
        if n_multi > 0:
            print(f"\n  Note: {n_multi} donor IDs contain '_' (single underscore).")
            print("  These are fine — the separator is '___' (triple underscore).")
            sample_ids = anno[DONOR_COL][multi_under].unique()[:5]
            print(f"  Sample donor IDs with '_': {sample_ids.tolist()}")

    # ── 3. Count integrality in pseudobulk CSV ────────────────────────────────
    sep("3. PSEUDOBULK COUNT INTEGRALITY")
    sample_vals = pb.values[:200, :20].ravel()
    frac_int    = np.mean(sample_vals == np.floor(sample_vals))
    print(f"  Sampled {len(sample_vals)} values")
    print(f"  Fraction that are whole numbers: {frac_int:.4f}")
    print(f"  Value range: [{pb.values.min():.2f}, {pb.values.max():,.0f}]")
    if frac_int < 0.99:
        print("  WARNING: many non-integer values — data may not be raw counts")
    if pb.values.max() < 50:
        print("  WARNING: max value < 50 — data is likely normalised/log-transformed")

    float32_limit = 2 ** 24  # 16,777,216
    max_val = pb.values.max()
    print(f"\n  float32 exact-integer limit: {float32_limit:,}")
    print(f"  Max pseudobulk value:        {max_val:,.0f}")
    if max_val > float32_limit:
        print("  WARNING: max count exceeds float32 precision limit!")
        print("  Counts were accumulated as float32 in the Python script.")
        print("  Values above 2^24 may be rounded, inflating/deflating counts.")
    else:
        print("  OK: max count is within float32 exact-integer range")

    # ── 4. Load raw h5ad and spot-check sums ─────────────────────────────────
    sep("4. RAW h5ad vs PSEUDOBULK SPOT-CHECK")
    print(f"Loading {H5AD_PATH} (backed mode) ...")
    adata = ad.read_h5ad(H5AD_PATH, backed="r")
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Filter to AD/Control cells only (same as pipeline)
    keep_dx   = set([DX_CONTROL] + DX_AD_VALUES)
    mask      = adata.obs[DX_COL].isin(keep_dx)
    adata_flt = adata[mask].to_memory()
    print(f"  After dx filter: {adata_flt.n_obs:,} cells")

    genes = adata_flt.var_names.tolist()

    # Get the raw count matrix (CSR, cells × genes)
    X = adata_flt.layers[RAW_LAYER]
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    elif not isinstance(X, sp.csr_matrix):
        X = X.tocsr()

    # Sample integrality of raw counts
    sample_raw = X.data[:10_000] if X.nnz > 10_000 else X.data
    frac_int_raw = np.mean(sample_raw == np.floor(sample_raw))
    print(f"\n  Raw count integrality (sampled {len(sample_raw)} nonzeros): "
          f"{frac_int_raw:.4f}")
    print(f"  Raw count range: [{sample_raw.min():.2f}, {sample_raw.max():.2f}]")
    if frac_int_raw < 0.99 or sample_raw.max() < 50:
        print("  WARNING: raw values do not look like integer counts")

    # ── Build group keys (same logic as Python script) ────────────────────────
    obs       = adata_flt.obs.copy()
    obs["_gk"] = (obs[CELLTYPE_COL].astype(str) + "___" +
                  obs[DONOR_COL].astype(str))
    unique_grps = np.unique(obs["_gk"].values)
    print(f"\n  Total (celltype, donor) pseudobulk groups: {len(unique_grps)}")
    print(f"  Pseudobulk CSV columns: {pb.shape[1]}")
    if len(unique_grps) != pb.shape[1]:
        print("  WARNING: group count mismatch!")

    # ── Spot-check N groups ────────────────────────────────────────────────────
    print(f"\n  Spot-checking {N_SPOT_CHECK} random groups "
          "(float64 sum vs CSV value vs float32 sum) ...")

    rng      = random.Random(42)
    to_check = rng.sample(list(unique_grps), min(N_SPOT_CHECK, len(unique_grps)))

    gene_idx = {g: i for i, g in enumerate(genes)}
    pb_gene_idx = {g: i for i, g in enumerate(pb.index.tolist())}

    # Pick a few high-expression genes to make precision differences visible
    pb_gene_max   = pb.values.max(axis=1)
    top_gene_pos  = np.argsort(pb_gene_max)[::-1][:10]
    check_genes   = [pb.index[i] for i in top_gene_pos if pb.index[i] in gene_idx][:5]

    all_ok = True
    for gk in to_check:
        ct, donor = gk.split("___", 1)
        cell_mask = obs["_gk"] == gk
        cell_rows = np.where(cell_mask)[0]

        # Expected label in pseudobulk CSV
        pb_label = gk.replace("___", "_", 1)

        if pb_label not in pb.columns:
            print(f"  [{gk}]: label '{pb_label}' NOT FOUND in pseudobulk CSV")
            all_ok = False
            continue

        print(f"\n  Group: {ct} / {donor}  ({len(cell_rows)} cells)")

        for gene in check_genes:
            gi_raw = gene_idx[gene]
            gi_pb  = pb_gene_idx.get(gene)
            if gi_pb is None:
                continue

            # float64 sum of raw counts
            raw_sum_f64 = float(np.asarray(
                X[cell_rows, gi_raw].sum(), dtype=np.float64
            ).ravel()[0])

            # float32 sum (simulating pipeline)
            raw_sum_f32 = float(np.float32(np.asarray(
                X[cell_rows, gi_raw].sum(dtype=np.float32)
            ).ravel()[0]))

            # Value in pseudobulk CSV
            pb_val = float(pb.iloc[gi_pb][pb_label])

            match_f64 = abs(raw_sum_f64 - pb_val) < 0.5
            f32_loss  = abs(raw_sum_f32 - raw_sum_f64)

            status = "OK" if match_f64 else "MISMATCH"
            if not match_f64:
                all_ok = False

            print(f"    {gene:20s}  raw_f64={raw_sum_f64:>10.0f}  "
                  f"csv={pb_val:>10.0f}  f32_loss={f32_loss:>6.1f}  [{status}]")

    if all_ok:
        print("\n  All spot-checked values match (within rounding).")
    else:
        print("\n  WARNING: mismatches detected — investigate further.")

    sep("DONE")


if __name__ == "__main__":
    main()
