"""
diagnose_pseudobulk.py

Loads pseudobulk CSV + annotation CSV pairs and prints summary statistics
for each dataset/region combination.

Datasets analysed:
  MIT  — PFC  : mit_pfc_pseudobulk.csv.gz  + mit_pfc_pseudobulk_annotations.csv
  MIT  — MTC  : mit_mtc_pseudobulk.csv.gz  + mit_mtc_pseudobulk_annotations.csv
  SEA  — DFC  : sea_dfc_pseudobulk.csv.gz  + sea_dfc_pseudobulk_annotations.csv
  SEA  — MTG  : sea_mtg_pseudobulk.csv.gz  + sea_mtg_pseudobulk_annotations.csv

Usage:
    python diagnose_pseudobulk.py
    python diagnose_pseudobulk.py --out diagnose_pseudobulk_report.txt
"""

import argparse
import os
import sys
import textwrap
import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────

MIT_DIR = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/mit_pseudobulk")
SEA_DFC_DIR = os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/dfc_data_pseudobulk")
SEA_MTG_DIR = os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/mtg_data_pseudobulk")

DATASETS = [
    dict(label="MIT — PFC",
         expr=os.path.join(MIT_DIR, "mit_pfc_pseudobulk.csv.gz"),
         anno=os.path.join(MIT_DIR, "mit_pfc_pseudobulk_annotations.csv"),
         celltype_col="Gabitto_metacell_labels",
         donor_col="ROSMAP_IndividualID",
         dx_col="Diagnosis"),
    dict(label="MIT — MTC",
         expr=os.path.join(MIT_DIR, "mit_mtc_pseudobulk.csv.gz"),
         anno=os.path.join(MIT_DIR, "mit_mtc_pseudobulk_annotations.csv"),
         celltype_col="Gabitto_metacell_labels",
         donor_col="ROSMAP_IndividualID",
         dx_col="Diagnosis"),
    dict(label="SEA-AD — DFC",
         expr=os.path.join(SEA_DFC_DIR, "sea_dfc_pseudobulk.csv.gz"),
         anno=os.path.join(SEA_DFC_DIR, "sea_dfc_pseudobulk_annotations.csv"),
         celltype_col="Subclass",
         donor_col="Donor ID",
         dx_col="Dx"),
    dict(label="SEA-AD — MTG",
         expr=os.path.join(SEA_MTG_DIR, "sea_mtg_pseudobulk.csv.gz"),
         anno=os.path.join(SEA_MTG_DIR, "sea_mtg_pseudobulk_annotations.csv"),
         celltype_col="Subclass",
         donor_col="Donor ID",
         dx_col="Dx"),
]

# ── Helpers ───────────────────────────────────────────────────────────────────

SEP = "=" * 70

def fmt_num(x):
    return f"{x:,.0f}"

def five_num(arr, label=""):
    """Return a formatted string with min/p25/median/p75/max/mean."""
    a = np.asarray(arr, dtype=float)
    return (
        f"{label}"
        f"  min={fmt_num(a.min())}  p25={fmt_num(np.percentile(a,25))}"
        f"  median={fmt_num(np.median(a))}  p75={fmt_num(np.percentile(a,75))}"
        f"  max={fmt_num(a.max())}  mean={fmt_num(a.mean())}"
    )

def value_counts_str(series, indent=4):
    vc = series.value_counts().sort_index()
    pad = " " * indent
    return "\n".join(f"{pad}{k}: {v}" for k, v in vc.items())

def analyse(ds, out):
    label        = ds["label"]
    expr_path    = ds["expr"]
    anno_path    = ds["anno"]
    celltype_col = ds["celltype_col"]
    donor_col    = ds["donor_col"]
    dx_col       = ds["dx_col"]

    def p(s=""):
        print(s, file=out)

    p(SEP)
    p(f"  {label}")
    p(SEP)

    # ── Check files exist ──────────────────────────────────────────────────
    for path in (expr_path, anno_path):
        if not os.path.exists(path):
            p(f"  [MISSING] {path}")
            p()
            return

    # ── Load ──────────────────────────────────────────────────────────────
    p(f"  expr : {expr_path}")
    p(f"  anno : {anno_path}")
    p()

    mat  = pd.read_csv(expr_path, index_col=0)     # genes × samples
    anno = pd.read_csv(anno_path)

    n_genes   = mat.shape[0]
    n_samples = mat.shape[1]
    p(f"  Shape: {fmt_num(n_genes)} genes × {fmt_num(n_samples)} pseudobulk samples")

    # ── Counts per sample (column sums) ───────────────────────────────────
    col_sums = mat.sum(axis=0)
    p()
    p("  Counts per pseudobulk sample (library size):")
    p("  " + five_num(col_sums))

    # Samples with very low total counts
    low_thresh = 1_000
    n_low = (col_sums < low_thresh).sum()
    if n_low:
        p(f"  WARNING: {n_low} samples have < {fmt_num(low_thresh)} total counts")

    # ── Counts per gene (row sums) ─────────────────────────────────────────
    row_sums = mat.sum(axis=1)
    p()
    p("  Counts per gene (across all samples):")
    p("  " + five_num(row_sums))

    # Genes with zero counts across all samples
    n_zero_genes = (row_sums == 0).sum()
    p(f"  Genes with 0 counts across all samples: {fmt_num(n_zero_genes)} "
      f"({100 * n_zero_genes / n_genes:.1f}%)")

    # ── Sparsity ──────────────────────────────────────────────────────────
    n_zeros     = (mat == 0).values.sum()
    total_cells = mat.size
    sparsity    = 100 * n_zeros / total_cells
    p()
    p(f"  Matrix sparsity: {sparsity:.1f}% zeros  "
      f"({fmt_num(n_zeros)} / {fmt_num(total_cells)} entries)")

    # ── Annotation summary ────────────────────────────────────────────────
    p()
    p(f"  Annotation columns: {list(anno.columns)}")

    # Diagnosis
    if dx_col in anno.columns:
        p()
        p(f"  Diagnosis breakdown ({dx_col}):")
        p(value_counts_str(anno[dx_col]))
    else:
        p(f"  WARNING: dx column '{dx_col}' not found in annotations")

    # Donors
    if donor_col in anno.columns:
        n_donors = anno[donor_col].nunique()
        p()
        p(f"  Unique donors: {n_donors}")
        if dx_col in anno.columns:
            donor_dx = (anno.drop_duplicates(subset=donor_col)
                        .groupby(dx_col)[donor_col].count())
            p(f"  Donors per diagnosis group:")
            for dx, n in donor_dx.items():
                p(f"    {dx}: {n}")
    else:
        p(f"  WARNING: donor column '{donor_col}' not found in annotations")

    # Celltypes
    if celltype_col in anno.columns:
        n_cts = anno[celltype_col].nunique()
        p()
        p(f"  Celltypes ({n_cts} unique):")
        p(value_counts_str(anno[celltype_col]))

        # Samples per celltype × diagnosis
        if dx_col in anno.columns:
            p()
            p(f"  Samples per celltype × diagnosis:")
            ct_dx = (anno.groupby([celltype_col, dx_col])
                     .size()
                     .unstack(fill_value=0))
            p(ct_dx.to_string(index=True))

        # Library size stats broken down by celltype
        if "label" in anno.columns:
            p()
            p("  Median library size per celltype:")
            sample_totals = col_sums.rename("total_counts")
            merged = anno.set_index("label").join(sample_totals, how="left")
            if celltype_col in merged.columns:
                ct_lib = (merged.groupby(celltype_col)["total_counts"]
                          .median()
                          .sort_values(ascending=False))
                for ct, med in ct_lib.items():
                    p(f"    {ct}: {fmt_num(med)}")
    else:
        p(f"  WARNING: celltype column '{celltype_col}' not found in annotations")

    # Samples with no matching annotation
    if "label" in anno.columns:
        anno_labels = set(anno["label"])
        mat_labels  = set(mat.columns)
        unmatched   = mat_labels - anno_labels
        if unmatched:
            p()
            p(f"  WARNING: {len(unmatched)} matrix columns have no matching annotation row:")
            for u in sorted(unmatched):
                p(f"    {u}")

    p()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Diagnose pseudobulk output files")
    parser.add_argument("--out", default=None,
                        help="Write report to this file in addition to stdout")
    args = parser.parse_args()

    outputs = [sys.stdout]
    fh = None
    if args.out:
        fh = open(args.out, "w")
        outputs.append(fh)

    class MultiOut:
        def __init__(self, streams):
            self.streams = streams
        def write(self, s):
            for st in self.streams:
                st.write(s)
        def flush(self):
            for st in self.streams:
                st.flush()

    out = MultiOut(outputs)

    for ds in DATASETS:
        analyse(ds, out)

    if fh:
        fh.close()
        print(f"\nReport also written to: {args.out}")


if __name__ == "__main__":
    main()
