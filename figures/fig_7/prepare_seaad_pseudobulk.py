"""
prepare_seaad_pseudobulk.py

Generic pseudobulking script for SEA-AD .h5ad files (any region).
Replaces prepare_seaad_mtg_pseudobulk.py.

Loads a SEA-AD .h5ad object (optionally filtered to one brain region) and saves:
  1. A pseudobulked (genes x samples) compressed CSV — recommended default
  2. A pseudobulk annotation CSV with a recoded Dx column
  OR (with --no-pseudobulk):
  1. Raw sparse matrix in HDF5 format
  2. Fallback to gzipped .mtx if HDF5 is unavailable

Key defaults (all overridable via CLI):
  - Raw counts layer:  'UMIs'
  - Diagnosis column:  'Overall AD neuropathological Change'
      "Not AD"                         → Control
      "High" / "Intermediate" / "Low"  → Alzheimers
      anything else (e.g. Reference)   → excluded before pseudobulking
  - Cell type column:  "Subclass"
  - Donor column:      "Donor ID"
  - Output prefix:     "seaad"  (set --outprefix to e.g. "sea_mtg", "sea_dfc")

Usage:
    # MTG (region already isolated in the file):
    python prepare_seaad_pseudobulk.py \\
        --input     /path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \\
        --outdir    /path/to/output/ \\
        --outprefix sea_mtg

    # DFC from a multi-region file, filter on the fly:
    python prepare_seaad_pseudobulk.py \\
        --input        /path/to/SEAAD_all_regions.h5ad \\
        --outdir       /path/to/output/ \\
        --outprefix    sea_dfc \\
        --region-col   "Brain Region" \\
        --region-value DFC

    # Export raw sparse matrix instead of pseudobulking:
    python prepare_seaad_pseudobulk.py \\
        --input     /path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \\
        --outdir    /path/to/output/ \\
        --outprefix sea_mtg \\
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
RAW_LAYER        = "UMIs"
DX_COL           = "Overall AD neuropathological Change"
DX_CONTROL_VALUE = "Not AD"
DX_AD_VALUES     = ["High", "Intermediate", "Low"]
CELLTYPE_COL     = "Subclass"
DONOR_COL        = "Donor ID"

# Columns from .obs to keep in the cell-level metadata CSV.
# Missing columns are skipped with a warning — safe to use on any SEA-AD file.
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
        description="Pseudobulk a SEA-AD .h5ad file for the R DE pipeline"
    )
    p.add_argument("--input",             required=True,
                   help="Path to input .h5ad file")
    p.add_argument("--outdir",            required=True,
                   help="Output directory")
    p.add_argument("--outprefix",         default="seaad",
                   help="Prefix for all output filenames (default: 'seaad'). "
                        "Set to e.g. 'sea_mtg' or 'sea_dfc' to match region.")
    # Optional region filter (only applied when both args are provided)
    p.add_argument("--region-col",        default=None,
                   help="obs column containing brain region labels. "
                        "Required together with --region-value to filter.")
    p.add_argument("--region-value",      default=None,
                   help="Brain region value to retain (e.g. 'DFC', 'MTG'). "
                        "Required together with --region-col to filter.")
    p.add_argument("--layer",             default=RAW_LAYER,
                   help=f"AnnData layer with raw counts (default: '{RAW_LAYER}')")
    p.add_argument("--celltype-col",      default=CELLTYPE_COL,
                   help=f"obs column for cell type (default: '{CELLTYPE_COL}')")
    p.add_argument("--donor-col",         default=DONOR_COL,
                   help=f"obs column for donor ID (default: '{DONOR_COL}')")
    p.add_argument("--dx-col",            default=DX_COL,
                   help=f"obs column for diagnosis (default: '{DX_COL}')")
    p.add_argument("--dx-control-value",  default=DX_CONTROL_VALUE,
                   help=f"Value in dx-col for controls (default: '{DX_CONTROL_VALUE}')")
    p.add_argument("--dx-ad-values",      default=DX_AD_VALUES, nargs="+",
                   help=f"Values in dx-col for AD cases (default: {DX_AD_VALUES})")
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
    Sum counts per (celltype, donor) group in Python.

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


def save_pseudobulk_csv(pb_matrix, pb_meta, genes, outdir, prefix):
    t0 = time.time()
    df = pd.DataFrame(
        pb_matrix.astype(np.int32),
        index   = genes,
        columns = pb_meta["label"].values,
    )
    df.index.name = "Gene"
    out_path = os.path.join(outdir, f"{prefix}_pseudobulk.csv.gz")
    df.to_csv(out_path, compression="gzip")
    print(f"  Saved pseudobulk CSV ({df.shape[0]:,} genes × {df.shape[1]} samples) "
          f"in {time.time() - t0:.1f}s: {out_path}")
    return out_path


def save_sparse_h5(X_csr, genes, barcodes, outdir, prefix):
    """Write 10x-style HDF5; readable in R via BPCells or HDF5Array."""
    try:
        import h5py
    except ImportError:
        print("  h5py not installed — falling back to .mtx.gz")
        return None

    t0 = time.time()
    X_csc   = X_csr.T.tocsc()
    h5_path = os.path.join(outdir, f"{prefix}_counts.h5")

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


def save_sparse_mtx_fast(X_csr, genes, barcodes, outdir, prefix):
    """Fast gzipped MTX using numpy — avoids slow scipy.io.mmwrite."""
    import gzip
    t0 = time.time()
    X_coo    = X_csr.T.tocoo()
    mtx_path = os.path.join(outdir, f"{prefix}_counts.mtx.gz")

    with gzip.open(mtx_path, "wt") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{X_coo.shape[0]} {X_coo.shape[1]} {X_coo.nnz}\n")
        triplets = np.column_stack([
            X_coo.row + 1,
            X_coo.col + 1,
            X_coo.data.astype(np.int32),
        ])
        np.savetxt(f, triplets, fmt="%d %d %d")

    genes_path    = os.path.join(outdir, f"{prefix}_genes.txt")
    barcodes_path = os.path.join(outdir, f"{prefix}_barcodes.txt")
    with open(genes_path,    "w") as f: f.write("\n".join(genes)    + "\n")
    with open(barcodes_path, "w") as f: f.write("\n".join(barcodes) + "\n")

    print(f"  Saved MTX.gz in {time.time() - t0:.1f}s: {mtx_path}")
    return mtx_path, genes_path, barcodes_path


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args    = parse_args()
    prefix  = args.outprefix
    os.makedirs(args.outdir, exist_ok=True)
    t_start = time.time()

    # Validate region filter args: must provide both or neither
    if bool(args.region_col) != bool(args.region_value):
        sys.exit("ERROR: --region-col and --region-value must be provided together.")

    # 1. Load
    print(f"Loading {args.input} ...")
    adata = ad.read_h5ad(args.input, backed="r")
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # 2. Optional region filter
    if args.region_col and args.region_value:
        if args.region_col not in adata.obs.columns:
            sys.exit(
                f"ERROR: --region-col '{args.region_col}' not in .obs.\n"
                f"Available columns: {list(adata.obs.columns)}"
            )
        regions_found = adata.obs[args.region_col].unique().tolist()
        if args.region_value not in regions_found:
            sys.exit(
                f"ERROR: --region-value '{args.region_value}' not found in "
                f"'{args.region_col}'.\nValues present: {sorted(regions_found)}"
            )
        region_mask = adata.obs[args.region_col] == args.region_value
        n_before    = adata.n_obs
        adata       = adata[region_mask].to_memory()
        print(f"  Region filter '{args.region_value}': "
              f"{n_before:,} → {adata.n_obs:,} cells  ({time.time()-t_start:.1f}s)")
    else:
        adata = adata.to_memory()

    # 3. Validate required columns and layer
    obs_cols = list(adata.obs.columns)
    errors = []
    if args.dx_col not in obs_cols:
        errors.append(f"  --dx-col '{args.dx_col}' not found in .obs")
    if args.celltype_col not in obs_cols:
        errors.append(f"  --celltype-col '{args.celltype_col}' not found in .obs")
    if args.donor_col not in obs_cols:
        errors.append(f"  --donor-col '{args.donor_col}' not found in .obs")
    if args.layer and args.layer not in adata.layers:
        errors.append(f"  --layer '{args.layer}' not found in .layers")
    if errors:
        sys.exit(
            "ERROR: the following inputs were not found:\n" +
            "\n".join(errors) +
            f"\n\nAvailable .obs columns: {obs_cols}" +
            (f"\nAvailable layers: {list(adata.layers)}" if args.layer else "")
        )

    print(f"\n  '{args.dx_col}' value counts (before dx filter):")
    print(adata.obs[args.dx_col].value_counts().to_string())

    # 4. Filter to control and AD cells
    keep_values  = set([args.dx_control_value] + list(args.dx_ad_values))
    values_found = set(adata.obs[args.dx_col].unique().tolist())
    missing      = keep_values - values_found
    if missing:
        sys.exit(
            f"ERROR: dx values not found in '{args.dx_col}': {missing}\n"
            f"Values present: {sorted(values_found)}"
        )

    mask       = adata.obs[args.dx_col].isin(keep_values)
    n_excluded = int((~mask).sum())
    adata      = adata[mask]
    print(f"\n  Excluded {n_excluded:,} cells (outside control/AD groups)")
    print(f"  Retained {adata.n_obs:,} cells  ({time.time()-t_start:.1f}s)")

    # 5. Extract counts
    X = get_counts(adata, args.layer)
    warn_if_normalised(X)

    genes    = adata.var_names.tolist()
    barcodes = adata.obs_names.tolist()

    # 6. Save cell-level metadata
    missing_cols = [c for c in OBS_COLS_TO_KEEP if c not in adata.obs.columns]
    if missing_cols:
        print(f"  WARNING: OBS_COLS_TO_KEEP columns not found, skipped: {missing_cols}")
    keep_cols = [c for c in OBS_COLS_TO_KEEP if c in adata.obs.columns]
    meta = adata.obs[keep_cols].copy()
    meta.index.name = "Cell_ID"
    meta = meta.reset_index()

    meta_path = os.path.join(args.outdir, f"{prefix}_metadata.csv")
    meta.to_csv(meta_path, index=False)
    print(f"  Saved cell metadata ({meta.shape[0]:,} cells): {meta_path}")

    print(f"\n  Diagnosis breakdown:\n"
          f"{adata.obs[args.dx_col].value_counts().to_string()}")
    print(f"\n  Cell types ({adata.obs[args.celltype_col].nunique()} unique):\n"
          f"{adata.obs[args.celltype_col].value_counts().to_string()}")

    # 7a. Default: pseudobulk in Python
    if not args.no_pseudobulk:
        print("\nPseudobulking in Python ...")
        pb_matrix, pb_meta = pseudobulk_sparse(
            X, adata.obs, args.celltype_col, args.donor_col
        )

        recode = {args.dx_control_value: args.dx_control_label}
        for v in args.dx_ad_values:
            recode[v] = args.dx_ad_label

        dx_map = (adata.obs[[args.donor_col, args.dx_col]]
                  .drop_duplicates()
                  .set_index(args.donor_col)[args.dx_col]
                  .to_dict())
        pb_meta["Dx"] = pb_meta[args.donor_col].map(dx_map).map(recode)

        # Add per-donor metadata columns to the pseudobulk annotations
        donor_meta_cols = [
            c for c in OBS_COLS_TO_KEEP
            if c in adata.obs.columns
            and c not in (args.celltype_col, args.donor_col, args.dx_col)
        ]
        if donor_meta_cols:
            donor_meta = (
                adata.obs[[args.donor_col] + donor_meta_cols]
                .drop_duplicates(subset=args.donor_col)
                .reset_index(drop=True)
            )
            pb_meta = pb_meta.merge(donor_meta, on=args.donor_col, how="left")

        pb_path  = save_pseudobulk_csv(pb_matrix, pb_meta, genes, args.outdir, prefix)
        ann_path = os.path.join(args.outdir, f"{prefix}_pseudobulk_annotations.csv")
        pb_meta.to_csv(ann_path, index=False)
        print(f"  Saved pseudobulk annotations: {ann_path}")

        print("\n--- R loading snippet ---")
        print(f"""\
  library(data.table)
  rawcounts   <- fread("{pb_path}") |> tibble::column_to_rownames("Gene")
  sea_anno_pb <- fread("{ann_path}")
  # Dx column contains '{args.dx_control_label}' / '{args.dx_ad_label}'
  # Pass directly to run_edger_dx_de() / run_deseq2_dx_de()
""")

    # 7b. Alternative: export raw sparse matrix
    else:
        print("\nSaving raw sparse matrix ...")
        h5_path = save_sparse_h5(X, genes, barcodes, args.outdir, prefix)
        if h5_path is None:
            save_sparse_mtx_fast(X, genes, barcodes, args.outdir, prefix)

        print("\n--- R loading snippet (BPCells) ---")
        fallback = os.path.join(args.outdir, f"{prefix}_counts.mtx.gz")
        print(f"""\
  library(BPCells)
  cell_expr <- BPCells::open_matrix_10x_hdf5("{h5_path or fallback}")
""")

    print(f"\nTotal time: {time.time() - t_start:.1f}s")
    print("Done.")


if __name__ == "__main__":
    main()
