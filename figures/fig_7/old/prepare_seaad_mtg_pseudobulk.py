"""
prepare_seaad_mtg_pseudobulk.py

Loads the SEA-AD MTG .h5ad object and saves:
  1. A pseudobulked (genes x samples) compressed CSV — recommended default
  2. A pseudobulk annotation CSV with a recoded Dx column
  OR (with --no-pseudobulk):
  1. Raw sparse matrix in HDF5 format
  2. Fallback to gzipped .mtx if HDF5 is unavailable

Key differences from prepare_mit_multiome_and_pseudobulk.py:
  - No brain-region filtering (file is already MTG-only)
  - Raw counts live in the 'UMIs' layer, not .X
  - Diagnosis column: 'Overall AD neuropathological Change'
      "Not AD"                         → Control   (configurable via --dx-control-label)
      "High" / "Intermediate" / "Low"  → Alzheimers (configurable via --dx-ad-label)
      anything else (e.g. Reference)   → excluded before pseudobulking
  - Cell type column: "Subclass"
  - Donor column:     "Donor ID"
  - Output prefix:    sea_mtg_

Usage:
    # Recommended: pseudobulk in Python, load direct in R
    python prepare_seaad_mtg_pseudobulk.py \
        --input  /path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
        --outdir /path/to/output/

    # Export raw sparse matrix instead
    python prepare_seaad_mtg_pseudobulk.py \
        --input  /path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
        --outdir /path/to/output/ \
        --no-pseudobulk

Dependencies:
    pip install anndata scipy pandas numpy h5py
"""

import argparse
import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


# ── Placeholders ─────────────────────────────────────────────────────────────
# Raw counts are in the 'UMIs' layer, not .X
RAW_LAYER = "UMIs"

DX_COL           = "Overall AD neuropathological Change"
DX_CONTROL_VALUE = "Not AD"
DX_AD_VALUES     = ["High", "Intermediate", "Low"]

CELLTYPE_COL = "Subclass"
DONOR_COL    = "Donor ID"

# Columns from .obs to keep in the cell-level metadata CSV
OBS_COLS_TO_KEEP = [
    "Subclass",
    "Supertype",
    "Donor ID",
    "Overall AD neuropathological Change",
    "Cognitive Status",
    "Brain Region",
    "Sex",
    "Age at Death",
    "APOE Genotype",
    "Primary Study Name",
]
# ── End placeholders ──────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(
        description="Pseudobulk SEA-AD MTG .h5ad for R DE pipeline"
    )
    p.add_argument("--input",             required=True,
                   help="Path to input .h5ad file")
    p.add_argument("--outdir",            required=True,
                   help="Output directory")
    p.add_argument("--layer",             default=RAW_LAYER,
                   help=f"AnnData layer with raw counts (default: '{RAW_LAYER}')")
    p.add_argument("--celltype-col",      default=CELLTYPE_COL,
                   help=f"obs column for cell type (default: '{CELLTYPE_COL}')")
    p.add_argument("--donor-col",         default=DONOR_COL,
                   help=f"obs column for donor ID (default: '{DONOR_COL}')")
    p.add_argument("--dx-col",            default=DX_COL,
                   help=f"obs column for diagnosis (default: '{DX_COL}')")
    p.add_argument("--dx-control-value",  default=DX_CONTROL_VALUE,
                   help=f"Value in dx-col indicating controls (default: '{DX_CONTROL_VALUE}')")
    p.add_argument("--dx-ad-values",      default=DX_AD_VALUES, nargs="+",
                   help=f"Values in dx-col indicating AD cases (default: {DX_AD_VALUES})")
    p.add_argument("--dx-control-label",  default="Control",
                   help="Label written to output Dx column for controls (default: 'Control')")
    p.add_argument("--dx-ad-label",       default="Alzheimers",
                   help="Label written to output Dx column for AD cases (default: 'Alzheimers')")
    p.add_argument("--no-pseudobulk",     action="store_true",
                   help="Skip pseudobulking; save raw sparse matrix instead")
    return p.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_counts(adata, layer):
    """Return count matrix as CSR sparse (cells x genes)."""
    X = adata.layers[layer] if layer else adata.X
    if not sp.issparse(X):
        print("  Converting dense matrix to sparse CSR ...")
        X = sp.csr_matrix(X)
    elif not isinstance(X, sp.csr_matrix):
        X = X.tocsr()
    return X


def warn_if_normalised(X):
    sample = X.data[:10_000] if X.nnz > 10_000 else X.data
    if len(sample) == 0:
        return
    if sample.max() < 50 or not np.allclose(sample, sample.astype(int)):
        print(
            f"  WARNING: max sampled value = {sample.max():.3f}. "
            "Matrix may be normalised/log-transformed. DE tools require raw integer counts."
        )


def pseudobulk_sparse(X_csr, obs, celltype_col, donor_col):
    """
    Sum counts per (celltype, donor) group in Python, keeping data sparse
    until the final group sums.

    Returns
    -------
    pb_matrix : np.ndarray float32, shape (n_genes, n_pseudobulk_samples)
    pb_meta   : pd.DataFrame
    """
    t0 = time.time()

    group_keys  = (obs[celltype_col].astype(str) + "___" +
                   obs[donor_col].astype(str)).values
    unique_grps, inv = np.unique(group_keys, return_inverse=True)

    n_genes   = X_csr.shape[1]
    n_samples = len(unique_grps)
    pb_matrix = np.zeros((n_genes, n_samples), dtype=np.float32)

    sort_order = np.argsort(inv)
    sorted_inv = inv[sort_order]
    boundaries = np.concatenate(
        [[0], np.where(np.diff(sorted_inv))[0] + 1, [len(sorted_inv)]]
    )

    for k in range(n_samples):
        rows = sort_order[boundaries[k]:boundaries[k + 1]]
        pb_matrix[:, k] = np.asarray(X_csr[rows, :].sum(axis=0)).ravel()

    print(f"  Pseudobulked {n_samples} samples in {time.time() - t0:.1f}s")

    pb_meta = pd.DataFrame({
        "label":       [g.replace("___", "_", 1) for g in unique_grps],
        celltype_col:  [g.split("___")[0] for g in unique_grps],
        donor_col:     [g.split("___")[1] for g in unique_grps],
    })
    return pb_matrix, pb_meta


def save_pseudobulk_csv(pb_matrix, pb_meta, genes, outdir):
    t0 = time.time()
    df = pd.DataFrame(
        pb_matrix.astype(np.int32),
        index   = genes,
        columns = pb_meta["label"].values,
    )
    df.index.name = "Gene"
    out_path = os.path.join(outdir, "sea_mtg_pseudobulk.csv.gz")
    df.to_csv(out_path, compression="gzip")
    print(f"  Saved pseudobulk CSV ({df.shape[0]:,} genes × {df.shape[1]} samples) "
          f"in {time.time() - t0:.1f}s: {out_path}")
    return out_path


def save_sparse_h5(X_csr, genes, barcodes, outdir):
    """Write 10x-style HDF5; readable in R via BPCells or HDF5Array."""
    try:
        import h5py
    except ImportError:
        print("  h5py not installed — falling back to .mtx.gz")
        return None

    t0 = time.time()
    X_csc   = X_csr.T.tocsc()
    h5_path = os.path.join(outdir, "sea_mtg_counts.h5")

    with h5py.File(h5_path, "w") as f:
        grp = f.create_group("matrix")
        grp.create_dataset("data",     data=X_csc.data,
                           compression="gzip", compression_opts=4, chunks=True)
        grp.create_dataset("indices",  data=X_csc.indices,
                           compression="gzip", compression_opts=4, chunks=True)
        grp.create_dataset("indptr",   data=X_csc.indptr,
                           compression="gzip", compression_opts=4, chunks=True)
        grp.create_dataset("shape",    data=np.array(X_csc.shape, dtype=np.int32))
        dt = h5py.string_dtype(encoding="ascii")
        grp.create_dataset("features", data=np.array(genes,    dtype="S"), dtype=dt)
        grp.create_dataset("barcodes", data=np.array(barcodes, dtype="S"), dtype=dt)

    print(f"  Saved HDF5 sparse matrix in {time.time() - t0:.1f}s: {h5_path}")
    return h5_path


def save_sparse_mtx_fast(X_csr, genes, barcodes, outdir):
    """Fast gzipped MTX using numpy — avoids slow scipy.io.mmwrite."""
    import gzip
    t0 = time.time()
    X_coo    = X_csr.T.tocoo()
    mtx_path = os.path.join(outdir, "sea_mtg_counts.mtx.gz")

    with gzip.open(mtx_path, "wt") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{X_coo.shape[0]} {X_coo.shape[1]} {X_coo.nnz}\n")
        triplets = np.column_stack([
            X_coo.row + 1,
            X_coo.col + 1,
            X_coo.data.astype(np.int32),
        ])
        np.savetxt(f, triplets, fmt="%d %d %d")

    genes_path    = os.path.join(outdir, "sea_mtg_genes.txt")
    barcodes_path = os.path.join(outdir, "sea_mtg_barcodes.txt")
    with open(genes_path,    "w") as f: f.write("\n".join(genes)    + "\n")
    with open(barcodes_path, "w") as f: f.write("\n".join(barcodes) + "\n")

    print(f"  Saved MTX.gz in {time.time() - t0:.1f}s: {mtx_path}")
    return mtx_path, genes_path, barcodes_path


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    t_start = time.time()

    # 1. Load — backed="r" memory-maps the file, only filtered cells hit RAM
    print(f"Loading {args.input} ...")
    adata = ad.read_h5ad(args.input, backed="r")
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # 2. Validate diagnosis column and report value counts
    if args.dx_col not in adata.obs.columns:
        sys.exit(
            f"ERROR: '{args.dx_col}' not in .obs.\n"
            f"Available columns: {list(adata.obs.columns)}"
        )
    print(f"\n  '{args.dx_col}' value counts (before filtering):")
    print(adata.obs[args.dx_col].value_counts().to_string())

    # 3. Filter to control and AD cells — excludes Reference and any other values
    keep_values = set([args.dx_control_value] + list(args.dx_ad_values))
    values_found = set(adata.obs[args.dx_col].unique().tolist())
    missing = keep_values - values_found
    if missing:
        sys.exit(
            f"ERROR: the following dx values were not found in '{args.dx_col}': {missing}\n"
            f"Values present: {sorted(values_found)}"
        )

    mask       = adata.obs[args.dx_col].isin(keep_values)
    n_excluded = int((~mask).sum())
    adata      = adata[mask].to_memory()
    print(f"\n  Excluded {n_excluded:,} cells (outside control/AD groups)")
    print(f"  Retained {adata.n_obs:,} cells  ({time.time()-t_start:.1f}s)")

    # 4. Extract counts from UMIs layer as CSR
    X = get_counts(adata, args.layer)
    warn_if_normalised(X)

    genes    = adata.var_names.tolist()
    barcodes = adata.obs_names.tolist()

    # 5. Save cell-level metadata
    missing_cols = [c for c in OBS_COLS_TO_KEEP if c not in adata.obs.columns]
    if missing_cols:
        print(f"  WARNING: OBS_COLS_TO_KEEP columns not found, skipped: {missing_cols}")
    keep_cols = [c for c in OBS_COLS_TO_KEEP if c in adata.obs.columns]
    meta = adata.obs[keep_cols].copy()
    meta.index.name = "Cell_ID"
    meta = meta.reset_index()

    meta_path = os.path.join(args.outdir, "sea_mtg_metadata.csv")
    meta.to_csv(meta_path, index=False)
    print(f"  Saved cell metadata ({meta.shape[0]:,} cells): {meta_path}")

    print(f"\n  Diagnosis breakdown:\n"
          f"{adata.obs[args.dx_col].value_counts().to_string()}")
    print(f"\n  Cell types ({adata.obs[args.celltype_col].nunique()} unique):\n"
          f"{adata.obs[args.celltype_col].value_counts().to_string()}")

    # 6a. Default: pseudobulk in Python — R receives a small dense CSV
    if not args.no_pseudobulk:
        print("\nPseudobulking in Python ...")
        pb_matrix, pb_meta = pseudobulk_sparse(
            X, adata.obs, args.celltype_col, args.donor_col
        )

        # Build donor → original dx mapping, then recode to output labels
        recode = {args.dx_control_value: args.dx_control_label}
        for v in args.dx_ad_values:
            recode[v] = args.dx_ad_label

        dx_map = (adata.obs[[args.donor_col, args.dx_col]]
                  .drop_duplicates()
                  .set_index(args.donor_col)[args.dx_col]
                  .to_dict())
        pb_meta["Dx"] = pb_meta[args.donor_col].map(dx_map).map(recode)

        pb_path  = save_pseudobulk_csv(pb_matrix, pb_meta, genes, args.outdir)
        ann_path = os.path.join(args.outdir, "sea_mtg_pseudobulk_annotations.csv")
        pb_meta.to_csv(ann_path, index=False)
        print(f"  Saved pseudobulk annotations: {ann_path}")

        print("\n--- R loading snippet ---")
        print(f"""\
  library(data.table)
  rawcounts   <- fread("{pb_path}") |> tibble::column_to_rownames("Gene")
  sea_anno_pb <- fread("{ann_path}")
  # Dx column contains '{args.dx_control_label}' / '{args.dx_ad_label}'
  # Pass directly to run_edger_dx_de() / run_deseq2_dx_de()
  # Update SEA_EXPR_PATH / SEA_ANNO_PATH in the R pipeline to point to these files.
""")

    # 6b. Alternative: export raw sparse matrix
    else:
        print("\nSaving raw sparse matrix ...")
        h5_path = save_sparse_h5(X, genes, barcodes, args.outdir)
        if h5_path is None:
            save_sparse_mtx_fast(X, genes, barcodes, args.outdir)

        print("\n--- R loading snippet (BPCells) ---")
        print(f"""\
  library(BPCells)
  cell_expr <- BPCells::open_matrix_10x_hdf5("{h5_path or os.path.join(args.outdir, 'sea_mtg_counts.mtx.gz')}")
""")

    print(f"\nTotal time: {time.time() - t_start:.1f}s")
    print("Done.")


if __name__ == "__main__":
    main()
