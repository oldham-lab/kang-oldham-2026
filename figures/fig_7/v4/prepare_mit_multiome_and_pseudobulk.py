"""
prepare_mit_multiome_and_pseudobulk.py  (v4)

Loads a MIT_Multiome .h5ad object, filters to a specified brain region, and saves:
  1. A pseudobulked (genes x samples) compressed CSV — recommended default,
     R receives a dense matrix already in the right shape for edgeR/DESeq2
     without ever needing to densify a large sparse matrix
  2. Cell-level metadata as a .csv
  OR (with --no-pseudobulk):
  1. Raw sparse matrix in HDF5 format (fast binary, readable by R's BPCells)
  2. Fallback to a fast gzipped .mtx if HDF5 is unavailable

v4 change:
  Accepts --external-labels pointing to an external CSV (index = cell barcode)
  that contains a "Gabitto_metacell_labels" column.  When provided, these labels
  replace the h5ad's own celltype column for pseudobulking.  Cells whose barcodes
  are absent from the external table are excluded before pseudobulking.

Key bottlenecks fixed vs the previous version:
  - backed="r" on read_h5ad: memory-maps the file, only loads filtered cells
    into RAM rather than the full dataset
  - Never calls sp.csc_matrix(X.T): that triggers a full reindex of all nnz
    values (O(nnz) copy). We keep CSR and operate on rows directly.
  - Pseudobulking in Python (numpy grouped row sums) instead of densifying
    the full sparse matrix in R via as.matrix() — avoids the single biggest
    memory and time cost in the R pipeline
  - Replaces scipy.io.mmwrite (pure Python, extremely slow) with either HDF5
    (h5py, fast binary) or a numpy-based MTX writer as fallback

Usage:
    # Recommended: pseudobulk in Python using external celltype labels
    python prepare_mit_multiome_and_pseudobulk.py \
        --input           /path/to/mit_multiome.h5ad \
        --outdir          /path/to/output/ \
        --region-value    PFC \
        --external-labels /path/to/liu_gabitto_metacell_labels_DFC_lognorm.csv

    # Without external labels (uses h5ad's own celltype column)
    python prepare_mit_multiome_and_pseudobulk.py \
        --input           /path/to/mit_multiome.h5ad \
        --outdir          /path/to/output/ \
        --region-value    PFC

    # Export raw sparse matrix instead of pseudobulking
    python prepare_mit_multiome_and_pseudobulk.py \
        --input           /path/to/mit_multiome.h5ad \
        --outdir          /path/to/output/ \
        --region-value    PFC \
        --external-labels /path/to/liu_gabitto_metacell_labels_DFC_lognorm.csv \
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
BRAIN_REGION_COL = "BrainRegion"
# BRAIN_REGION_VALUE is intentionally not set here — pass --region-value on
# the command line (e.g. --region-value PFC) to avoid silent wrong-region runs.

# Set to a layer name (e.g. "counts") if raw counts are not in .X
RAW_LAYER = None

# Column in the external labels CSV that holds the celltype labels
EXTERNAL_LABEL_COL = "Gabitto_metacell_labels"

# Columns from .obs to keep in the metadata CSV
OBS_COLS_TO_KEEP = [
    "RNA.Subclass",         # original celltype label
    "ROSMAP_IndividualID",  # donor ID           → MIT_DONOR_COL in R
    "Pathology",            # AD / Control label → MIT_DX_COL    in R
    "BrainRegion"
]
# ── End placeholders ──────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(
        description="Filter MIT_Multiome .h5ad to a specified brain region and export for R DE pipeline"
    )
    p.add_argument("--input",            required=True,  help="Path to input .h5ad file")
    p.add_argument("--outdir",           required=True,  help="Output directory")
    p.add_argument("--region-col",       default=BRAIN_REGION_COL,
                   help=f"obs column containing brain region labels (default: '{BRAIN_REGION_COL}')")
    p.add_argument("--region-value",     required=True,
                   help="Brain region value to filter on (e.g. 'PFC', 'MTC')")
    p.add_argument("--layer",            default=RAW_LAYER,
                   help="AnnData layer with raw counts; omit to use .X")
    p.add_argument("--celltype-col",     default="RNA.Subclass",
                   help="obs column for celltypes (used when --external-labels is not provided)")
    p.add_argument("--donor-col",        default="ROSMAP_IndividualID")
    p.add_argument("--dx-col",           default="Pathology")
    p.add_argument("--external-labels",  default=None,
                   help=(
                       "Path to external CSV (cell barcode as index) containing "
                       f"'{EXTERNAL_LABEL_COL}' column.  When provided, these labels "
                       "replace the h5ad celltype column for pseudobulking.  Cells "
                       "absent from this table are excluded."
                   ))
    p.add_argument("--no-pseudobulk",    action="store_true",
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


def apply_external_labels(adata, labels_path):
    """
    Load external celltype labels CSV, join onto adata by cell barcode index,
    and return a filtered AnnData with the external labels added as
    EXTERNAL_LABEL_COL in .obs.  Cells not present in the labels table are dropped.
    """
    print(f"\nLoading external celltype labels from:\n  {labels_path}")
    ext = pd.read_csv(labels_path, index_col=0)

    if EXTERNAL_LABEL_COL not in ext.columns:
        sys.exit(
            f"ERROR: '{EXTERNAL_LABEL_COL}' not found in external labels CSV.\n"
            f"Available columns: {list(ext.columns)}"
        )

    ext_labels = ext[[EXTERNAL_LABEL_COL]]

    # Match on cell barcode (obs_names)
    common = adata.obs_names.intersection(ext_labels.index)
    n_before = adata.n_obs
    n_missing = n_before - len(common)

    if len(common) == 0:
        sys.exit(
            "ERROR: No cell barcodes overlap between h5ad and external labels CSV. "
            "Check that the index column of the CSV matches adata.obs_names."
        )
    if n_missing > 0:
        print(f"  WARNING: {n_missing:,} cells not found in external labels table — excluded.")

    adata = adata[common].copy()
    adata.obs[EXTERNAL_LABEL_COL] = ext_labels.loc[common, EXTERNAL_LABEL_COL].values
    print(f"  {len(common):,} cells retained after joining external labels "
          f"({n_before - len(common):,} dropped).")
    print(f"  External celltype breakdown ({adata.obs[EXTERNAL_LABEL_COL].nunique()} types):\n"
          f"{adata.obs[EXTERNAL_LABEL_COL].value_counts().to_string()}")
    return adata


def pseudobulk_sparse(X_csr, obs, celltype_col, donor_col):
    """
    Sum counts per (celltype, donor) group in Python, keeping data sparse
    until the final group sums — avoids any full densification.

    Uses sorted index groups + CSR row slicing so the inner loop is a
    sparse row-sum (handled by scipy's C backend) rather than a Python loop
    over individual cells.

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

    # Sort cells by group so we can slice contiguous blocks from the CSR matrix
    sort_order = np.argsort(inv)
    sorted_inv = inv[sort_order]
    boundaries = np.concatenate(
        [[0], np.where(np.diff(sorted_inv))[0] + 1, [len(sorted_inv)]]
    )

    for k in range(n_samples):
        rows = sort_order[boundaries[k]:boundaries[k + 1]]
        # X_csr[rows, :].sum(axis=0) is a sparse row sum — fast C-level op
        pb_matrix[:, k] = np.asarray(X_csr[rows, :].sum(axis=0)).ravel()

    print(f"  Pseudobulked {n_samples} samples in {time.time() - t0:.1f}s")

    pb_meta = pd.DataFrame({
        "label":       [g.replace("___", "_", 1) for g in unique_grps],
        celltype_col:  [g.split("___")[0] for g in unique_grps],
        donor_col:     [g.split("___")[1] for g in unique_grps],
    })
    return pb_matrix, pb_meta


def save_pseudobulk_csv(pb_matrix, pb_meta, genes, outdir, region_tag):
    t0 = time.time()
    df = pd.DataFrame(
        pb_matrix.astype(np.int32),
        index   = genes,
        columns = pb_meta["label"].values,
    )
    df.index.name = "Gene"
    out_path = os.path.join(outdir, f"mit_{region_tag}_pseudobulk.csv.gz")
    df.to_csv(out_path, compression="gzip")
    print(f"  Saved pseudobulk CSV ({df.shape[0]:,} genes × {df.shape[1]} samples) "
          f"in {time.time() - t0:.1f}s: {out_path}")
    return out_path


def save_sparse_h5(X_csr, genes, barcodes, outdir, region_tag):
    """
    Write 10x-style HDF5. Fastest format; readable in R via BPCells or
    HDF5Array without loading into memory.
    """
    try:
        import h5py
    except ImportError:
        print("  h5py not installed — falling back to .mtx.gz")
        return None

    t0 = time.time()
    # tocsc() is the one necessary transpose/reindex, and it's unavoidable
    # for the standard genes-x-cells column-major HDF5 layout.
    # We do it once here rather than multiple times.
    X_csc   = X_csr.T.tocsc()
    h5_path = os.path.join(outdir, f"mit_{region_tag}_counts.h5")

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


def save_sparse_mtx_fast(X_csr, genes, barcodes, outdir, region_tag):
    """
    Fast gzipped MTX using numpy savetxt — avoids scipy.io.mmwrite which is
    pure Python and orders of magnitude slower for large matrices.
    """
    import gzip
    t0 = time.time()
    X_coo    = X_csr.T.tocoo()   # genes x cells COO
    mtx_path = os.path.join(outdir, f"mit_{region_tag}_counts.mtx.gz")

    with gzip.open(mtx_path, "wt") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{X_coo.shape[0]} {X_coo.shape[1]} {X_coo.nnz}\n")
        triplets = np.column_stack([
            X_coo.row + 1,
            X_coo.col + 1,
            X_coo.data.astype(np.int32),
        ])
        np.savetxt(f, triplets, fmt="%d %d %d")

    genes_path    = os.path.join(outdir, f"mit_{region_tag}_genes.txt")
    barcodes_path = os.path.join(outdir, f"mit_{region_tag}_barcodes.txt")
    with open(genes_path,    "w") as f: f.write("\n".join(genes)    + "\n")
    with open(barcodes_path, "w") as f: f.write("\n".join(barcodes) + "\n")

    print(f"  Saved MTX.gz in {time.time() - t0:.1f}s: {mtx_path}")
    return mtx_path, genes_path, barcodes_path


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args   = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    # Derive a safe filename tag from the region value (e.g. "PFC" → "pfc")
    region_tag = args.region_value.lower().replace(" ", "_")
    t_start = time.time()

    # 1. Load — backed="r" memory-maps the file so we only pull filtered
    #    cells into RAM rather than loading the entire dataset first
    print(f"Loading {args.input} ...")
    adata = ad.read_h5ad(args.input, backed="r")
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # 2. Validate region column
    if args.region_col not in adata.obs.columns:
        sys.exit(
            f"ERROR: '{args.region_col}' not in .obs.\n"
            f"Available columns: {list(adata.obs.columns)}"
        )
    regions_found = adata.obs[args.region_col].unique().tolist()
    if args.region_value not in regions_found:
        sys.exit(
            f"ERROR: '{args.region_value}' not found in '{args.region_col}'.\n"
            f"Values present: {regions_found}"
        )

    # 3. Filter to region — slice obs index first (cheap), then .to_memory()
    #    loads only the selected rows from disk
    mask  = adata.obs[args.region_col] == args.region_value
    adata = adata[mask].to_memory()
    print(f"  After region filter: {adata.n_obs:,} cells  ({time.time()-t_start:.1f}s)")

    # 4. Optionally replace celltype labels with external Gabitto metacell labels
    if args.external_labels:
        adata = apply_external_labels(adata, args.external_labels)
        celltype_col = EXTERNAL_LABEL_COL
    else:
        celltype_col = args.celltype_col

    # 5. Extract counts as CSR
    X = get_counts(adata, args.layer)
    warn_if_normalised(X)

    genes    = adata.var_names.tolist()
    barcodes = adata.obs_names.tolist()

    # 6. Save metadata
    missing = [c for c in OBS_COLS_TO_KEEP if c not in adata.obs.columns]
    if missing:
        print(f"  WARNING: OBS_COLS_TO_KEEP not found, skipped: {missing}")
    keep_cols = [c for c in OBS_COLS_TO_KEEP if c in adata.obs.columns]
    # Always include the active celltype column in the metadata CSV
    if celltype_col not in keep_cols:
        keep_cols = [celltype_col] + keep_cols
    meta = adata.obs[keep_cols].copy()
    meta.index.name = "Cell_ID"
    meta = meta.reset_index()

    meta_path = os.path.join(args.outdir, f"mit_{region_tag}_metadata.csv")
    meta.to_csv(meta_path, index=False)
    print(f"  Saved cell metadata ({meta.shape[0]:,} cells): {meta_path}")

    if args.dx_col in meta.columns:
        print(f"  Diagnosis breakdown:\n{meta[args.dx_col].value_counts().to_string()}")
    if celltype_col in meta.columns:
        print(f"  Cell types ({meta[celltype_col].nunique()} unique):\n"
              f"{meta[celltype_col].value_counts().to_string()}")

    # 7a. Default: pseudobulk in Python — R loads a small dense CSV
    if not args.no_pseudobulk:
        print("\nPseudobulking in Python ...")
        pb_matrix, pb_meta = pseudobulk_sparse(
            X, adata.obs, celltype_col, args.donor_col
        )

        # Attach diagnosis to pseudobulk metadata
        dx_map = (adata.obs[[args.donor_col, args.dx_col]]
                  .drop_duplicates()
                  .set_index(args.donor_col)[args.dx_col]
                  .to_dict())
        pb_meta["Diagnosis"] = pb_meta[args.donor_col].map(dx_map)

        pb_path  = save_pseudobulk_csv(pb_matrix, pb_meta, genes, args.outdir, region_tag)
        ann_path = os.path.join(args.outdir, f"mit_{region_tag}_pseudobulk_annotations.csv")
        pb_meta.to_csv(ann_path, index=False)
        print(f"  Saved pseudobulk annotations: {ann_path}")

        print("\n--- R loading snippet ---")
        print(f"""\
  library(data.table)
  rawcounts   <- fread("{pb_path}") |> tibble::column_to_rownames("Gene")
  mit_anno_pb <- fread("{ann_path}")
  # rawcounts is now a genes x samples dense matrix — pass directly to
  # run_edger_dx_de() / run_deseq2_dx_de() in full_pipeline_ADvsCon.R
""")

    # 7b. Alternative: export raw sparse matrix (for manual pseudobulking in R)
    else:
        print("\nSaving raw sparse matrix ...")
        h5_path = save_sparse_h5(X, genes, barcodes, args.outdir, region_tag)
        if h5_path is None:
            save_sparse_mtx_fast(X, genes, barcodes, args.outdir, region_tag)

        print("\n--- R loading snippet (BPCells) ---")
        print(f"""\
  library(BPCells)
  cell_expr <- BPCells::open_matrix_10x_hdf5("{h5_path or os.path.join(args.outdir, f'mit_{region_tag}_counts.mtx.gz')}")
  # proceed with pseudobulking in R as in full_pipeline_ADvsCon.R
""")

    print(f"\nTotal time: {time.time() - t_start:.1f}s")
    print("Done.")


if __name__ == "__main__":
    main()
