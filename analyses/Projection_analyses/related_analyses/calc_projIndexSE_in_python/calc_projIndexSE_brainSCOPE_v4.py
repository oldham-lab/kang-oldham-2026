"""
calculate_se_by_cohort_celltype_disorder.py

Calculates the standard error (SE) of gene expression across all nuclei
(cells) sharing the same Cohort, Disorder, CellType, and gene Module, for the
CMC and SZBDMulti-Seq cohorts.

Four output tables are produced (one per cohort x disorder combination):
  Rows    = gene modules  (keys of the module JSON file, sorted numerically)
  Columns = cell types    (derived from column headers of each matrix)
  Values  = SE of expression across all cells in that (cohort, disorder, celltype, module) group

Data sources:
  Metadata : https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt
  Matrices : <matrix_dir>/<Cohort>/<Individual_ID>-annotated_matrix.txt.gz
             Rows = genes, Columns = <barcode>_<CellType>

Module JSON format (--modules):
  {
    "module_1": ["GENE_A", "GENE_B", ...],
    "module_2": ["GENE_C", ...],
    ...
  }
  Keys are module names (become row labels); values are gene lists.
  Genes absent from the expression matrices are silently ignored per module.

Usage
-----
  python calc_projIndexSE_brainSCOPE.py \\
      --metadata   PEC2_sample_metadata.txt \\
      --matrix_dir /path/to/annotated_matrices \\
      --modules    gene_modules.json \\
      --output_dir ./results/ \\
      [--normalize] [--verbose]

Output files (written to --output_dir):
  CMC_control_SE.csv            -- rows: modules, columns: cell types
  CMC_Schizophrenia_SE.csv      -- rows: modules, columns: cell types
  SZBDMulti-Seq_control_SE.csv          -- rows: modules, columns: cell types
  SZBDMulti-Seq_Schizophrenia_SE.csv    -- rows: modules, columns: cell types

Notes
-----
  - Only CMC and SZBDMulti-Seq cohorts are processed.
  - Only samples whose Disorder is 'control' or 'Schizophrenia' are used.
  - Raw counts are used by default; pass --normalize for log1p(CPM).
  - SE is computed with Bessel correction: SD = sqrt(Var_sample), SE = SD/sqrt(n).
  - Each (gene, cell) pair within a module is treated as an individual
    observation. For a module with G genes and C cells the accumulator
    receives G × C values per sample. SE is then computed across all such
    values pooled over every sample in the (cohort, disorder, celltype) group.
    This is equivalent to ravelling (flattening) the gene × cell sub-matrix
    and computing SE on the resulting 1-D vector.
  - Gene lists in the module JSON with zero overlap with any matrix will
    produce NaN SE values for that module (flagged in the log).
  - Module rows are sorted by the numeric portion of their name (e.g.
    module_1 < module_2 < ... < module_10), falling back to lexicographic
    order for names with no embedded integer.
"""

import argparse
import gzip
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# Hard-coded per requirements
ALLOWED_COHORTS   = {"CMC", "SZBDMulti-Seq"}
ALLOWED_DISORDERS = {"control", "Schizophrenia"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_celltype(barcode: str, sep: str = "_") -> str:
    """
    Extract cell-type label from a barcode column header.

    Column format: <barcode><sep><CellType>  e.g. AAACCCAAGCATCACG-1_Ex-L2-3
    Everything after the first occurrence of `sep` is the cell-type label.
    Returns "Unknown" if the separator is not found.

    R automatically appends decimal suffixes (.1, .2, ...) to deduplicate
    repeated column names when a matrix is written from R. Any such suffix
    is stripped so that e.g. "Ex-L2-3.1" and "Ex-L2-3.14" are both
    normalised to "Ex-L2-3".
    """
    idx = barcode.find(sep)
    if idx == -1:
        return "Unknown"
    celltype = barcode[idx + 1:]
    # Strip trailing R-style deduplication suffix: one or more ".digits" groups
    celltype = re.sub(r"(\.[0-9]+)+$", "", celltype)
    return celltype


def detect_celltype_sep(columns, candidates=("_", " ", "\t", "|", ":", ";")):
    """
    Inspect actual column headers and return the first candidate separator
    that appears in the majority of them.  Falls back to '_' if none match.
    """
    for sep in candidates:
        hits = sum(1 for c in columns if sep in c)
        if hits >= len(columns) * 0.5:
            return sep
    return "_"


def find_matrix_file(matrix_dir: Path, cohort: str, sample_id: str) -> "Path | None":
    """Return path to a sample's gzipped annotated matrix, or None if not found."""
    candidates = [
        matrix_dir / cohort / f"{sample_id}-annotated_matrix.txt.gz",
        matrix_dir / f"{sample_id}-annotated_matrix.txt.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def read_matrix(path: Path) -> pd.DataFrame:
    """
    Read a gzipped gene x cell expression matrix.
    Returns DataFrame: index = genes, columns = cell barcodes (with cell-type suffix).
    Always logs the first 5 column headers so separator issues are immediately visible.
    """
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)
    log.info("  Shape: %d genes x %d cells | first 5 cols: %s",
             df.shape[0], df.shape[1], list(df.columns[:5]))
    return df


def log1p_cpm_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Apply log1p(CPM) normalization column-wise (per cell)."""
    col_sums = df.sum(axis=0).replace(0, 1)
    return np.log1p(df.div(col_sums, axis=1) * 1_000_000)


def module_sort_key(name: str):
    """
    Sort key that orders module names by their embedded integer (numeric sort),
    falling back to the full string for names without an integer.

    Examples:  module_1 < module_2 < module_10  (not module_1 < module_10 < module_2)
    """
    m = re.search(r"\d+", name)
    return (int(m.group()), name) if m else (float("inf"), name)


# ---------------------------------------------------------------------------
# Streaming statistics accumulator
# ---------------------------------------------------------------------------

class StreamingStats:
    """
    Accumulates sufficient statistics (sum_x, sum_x2, n) for later SE
    calculation, keyed by (cohort, disorder, celltype, module).

    Every individual (gene, cell) value in the module sub-matrix is treated
    as a separate observation. For a module with G genes and C cells the
    accumulator receives G × C values per sample — equivalent to ravelling
    (np.ravel) the gene × cell sub-matrix into a 1-D array and accumulating
    stats on that. No per-gene or per-cell averaging is performed before
    the statistics are recorded.
    """

    def __init__(self):
        # key: (cohort, disorder, celltype, module_name)
        # value: {"sum_x": float, "sum_x2": float, "n": int}
        self._data: dict = defaultdict(lambda: {"sum_x": 0.0, "sum_x2": 0.0, "n": 0})

    def add_sample(
        self,
        df: pd.DataFrame,
        cohort: str,
        disorder: str,
        modules: dict,
        celltype_sep: str,
    ) -> None:
        """
        Ingest one sample's expression matrix into the running statistics.

        For each (cohort, disorder, celltype, module) group the gene × cell
        sub-matrix is ravelled into a flat array of individual expression
        values, and sum_x / sum_x2 / n are updated from that array.

        Parameters
        ----------
        df       : genes (rows) x cells (columns), optionally normalized.
        cohort   : 'CMC' or 'SZBDMulti-Seq'
        disorder : 'control' or 'Schizophrenia'
        modules  : {module_name: [gene, ...]}
        """
        all_genes = set(df.index)

        # Auto-detect the separator actually present in this matrix's column headers.
        # The user-supplied celltype_sep is used as a hint/first candidate, but we
        # verify it appears in the data and fall back to auto-detection if not.
        effective_sep = detect_celltype_sep(df.columns, candidates=(celltype_sep, "_", " ", "|", ":", ";"))
        if effective_sep != celltype_sep:
            log.warning(
                "  celltype_sep %r not found in column headers; auto-detected %r instead. "
                "First 3 columns: %s",
                celltype_sep, effective_sep, list(df.columns[:3]),
            )

        # Group cell-barcode columns by cell type
        ct_groups: dict = defaultdict(list)
        for col in df.columns:
            ct_groups[parse_celltype(col, effective_sep)].append(col)

        for mod_name, gene_list in modules.items():
            present_genes = [g for g in gene_list if g in all_genes]
            if not present_genes:
                continue  # no overlap for this module in this sample

            # Sub-matrix for module genes: shape (n_genes x n_cells_in_sample)
            mod_expr = df.loc[present_genes]

            for ct, cols in ct_groups.items():
                # Ravel the (n_genes x n_cells_in_ct) sub-matrix into a
                # flat 1-D array of individual expression values.
                # n = n_present_genes * n_cells_in_ct
                values = mod_expr[cols].to_numpy(dtype=np.float64).ravel()

                key   = (cohort, disorder, ct, mod_name)
                entry = self._data[key]
                entry["sum_x"]  += float(values.sum())
                entry["sum_x2"] += float((values ** 2).sum())
                entry["n"]      += len(values)

    def compute_se(self) -> pd.DataFrame:
        """
        Finalize accumulated statistics and return a long-form DataFrame:
          Cohort | Disorder | CellType | Module | N_observations | Mean | SD | SE

        N_observations = total number of (gene, cell) pairs accumulated,
        i.e. sum over all samples of (n_present_module_genes × n_cells_in_ct).
        """
        rows = []
        for (cohort, disorder, ct, mod), entry in self._data.items():
            n = entry["n"]
            if n < 2:
                log.warning(
                    "Group (%s, %s, %s, %s) has only n=%d observation; SE set to NaN.",
                    cohort, disorder, ct, mod, n,
                )
                rows.append((cohort, disorder, ct, mod, n, np.nan, np.nan, np.nan))
                continue

            s    = entry["sum_x"]
            s2   = entry["sum_x2"]
            mean = s / n
            # Bessel-corrected sample variance
            var  = max((s2 - s**2 / n) / (n - 1), 0.0)
            sd   = np.sqrt(var)
            se   = sd / np.sqrt(n)
            rows.append((cohort, disorder, ct, mod, n, mean, sd, se))

        return pd.DataFrame(
            rows,
            columns=["Cohort", "Disorder", "CellType", "Module",
                     "N_observations", "Mean", "SD", "SE"],
        )


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def make_wide_table(long_df: pd.DataFrame, cohort: str, disorder: str) -> pd.DataFrame:
    """
    Pivot the long SE table for one (cohort, disorder) pair into a wide matrix:
      Index   = Module  (rows, sorted numerically by embedded integer)
      Columns = CellType
      Values  = SE
    """
    subset = long_df[
        (long_df["Cohort"] == cohort) & (long_df["Disorder"] == disorder)
    ][["Module", "CellType", "SE"]]

    if subset.empty:
        return pd.DataFrame()

    wide = subset.pivot(index="Module", columns="CellType", values="SE")
    wide.index.name   = "Module"
    wide.columns.name = None

    # Sort rows numerically by embedded integer in module name
    sorted_index = sorted(wide.index, key=module_sort_key)
    wide = wide.loc[sorted_index]

    return wide


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=(
            "Compute module-level SE of expression per (Cohort x Disorder x CellType), "
            "for CMC and SZBDMulti-Seq cohorts. "
            "Outputs four CSV tables: one per cohort x disorder combination."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--metadata", required=True,
                   help="Path to PEC2_sample_metadata.txt (tab-separated).")
    p.add_argument("--matrix_dir", required=True,
                   help="Root directory containing cohort sub-dirs with "
                        "*-annotated_matrix.txt.gz files.")
    p.add_argument("--modules", required=True,
                   help="Path to gene module JSON file. "
                        "Format: {module_name: [gene, ...], ...}")
    p.add_argument("--output_dir", default=".",
                   help="Directory for output CSV files (default: current dir).")
    p.add_argument("--sample_col",   default="Individual_ID",
                   help="Metadata column matching sample file names (default: Individual_ID).")
    p.add_argument("--cohort_col",   default="Cohort",
                   help="Metadata column for cohort (default: Cohort).")
    p.add_argument("--disorder_col", default="Disorder",
                   help="Metadata column for disorder (default: Disorder).")
    p.add_argument("--celltype_sep", default="_",
                   help="Separator between barcode and cell-type label in column "
                        "headers (default: '_').")
    p.add_argument("--normalize", action="store_true",
                   help="Apply log1p(CPM) normalization per cell before computing SE.")
    p.add_argument("--verbose", action="store_true",
                   help="Enable DEBUG-level logging.")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    # ------------------------------------------------------------------
    # 1. Load gene modules
    # ------------------------------------------------------------------
    log.info("Loading gene modules from %s", args.modules)
    with open(args.modules) as fh:
        modules: dict = json.load(fh)
    log.info("  %d modules loaded.", len(modules))

    # ------------------------------------------------------------------
    # 2. Load and filter sample metadata
    # ------------------------------------------------------------------
    log.info("Loading metadata from %s", args.metadata)
    meta = pd.read_csv(args.metadata, sep="\t")

    required_cols = {args.sample_col, args.cohort_col, args.disorder_col}
    missing = required_cols - set(meta.columns)
    if missing:
        log.error("Metadata is missing required columns: %s", missing)
        sys.exit(1)

    meta = meta[
        meta[args.cohort_col].isin(ALLOWED_COHORTS) &
        meta[args.disorder_col].isin(ALLOWED_DISORDERS)
    ].copy()

    log.info(
        "  %d samples retained after filtering to cohorts %s and disorders %s.",
        len(meta), sorted(ALLOWED_COHORTS), sorted(ALLOWED_DISORDERS),
    )

    if meta.empty:
        log.error(
            "No samples remain after filtering. "
            "Verify that cohort and disorder labels in the metadata match "
            "exactly: cohorts=%s, disorders=%s.",
            sorted(ALLOWED_COHORTS), sorted(ALLOWED_DISORDERS),
        )
        sys.exit(1)

    # ------------------------------------------------------------------
    # 3. Stream through matrices and accumulate statistics
    # ------------------------------------------------------------------
    matrix_dir = Path(args.matrix_dir)
    stats   = StreamingStats()
    n_loaded = n_skipped = 0

    for _, row in meta.iterrows():
        sample_id = str(row[args.sample_col])
        cohort    = str(row[args.cohort_col])
        disorder  = str(row[args.disorder_col])

        mat_path = find_matrix_file(matrix_dir, cohort, sample_id)
        if mat_path is None:
            log.warning(
                "Matrix file not found for sample '%s' (cohort '%s'); skipping.",
                sample_id, cohort,
            )
            n_skipped += 1
            continue

        log.info("[%-18s | %-15s | %s]  Reading %s",
                 cohort, disorder, sample_id, mat_path.name)

        try:
            df = read_matrix(mat_path)
        except Exception as exc:
            log.error("  Failed to read %s: %s", mat_path, exc)
            n_skipped += 1
            continue

        if args.normalize:
            df = log1p_cpm_normalize(df)

        stats.add_sample(
            df=df,
            cohort=cohort,
            disorder=disorder,
            modules=modules,
            celltype_sep=args.celltype_sep,
        )
        n_loaded += 1

    log.info("Matrix loading complete: %d processed, %d skipped.", n_loaded, n_skipped)

    if n_loaded == 0:
        log.error("No matrices were successfully loaded. Cannot compute SE.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 4. Compute SE and write wide output tables (one per cohort x disorder)
    # ------------------------------------------------------------------
    log.info("Computing standard errors across all groups ...")
    long_df = stats.compute_se()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    output_files = []
    for cohort in sorted(ALLOWED_COHORTS):
        for disorder in sorted(ALLOWED_DISORDERS):
            wide = make_wide_table(long_df, cohort, disorder)
            if wide.empty:
                log.warning(
                    "No data accumulated for cohort '%s' disorder '%s'; skipping output.",
                    cohort, disorder,
                )
                continue

            # Sanitize cohort name for use in filename (replace spaces/special chars)
            cohort_safe = cohort.replace(" ", "_")
            out_path = out_dir / f"{cohort_safe}_{disorder}_SE.csv"
            wide.to_csv(out_path, float_format="%.6g")
            output_files.append(str(out_path))
            log.info(
                "  Wrote %-50s  (%d modules x %d cell types)",
                str(out_path), wide.shape[0], wide.shape[1],
            )

    # ------------------------------------------------------------------
    # 5. Print summary to stdout
    # ------------------------------------------------------------------
    n_null = long_df["SE"].isna().sum()
    summary = (
        long_df.groupby(["Cohort", "Disorder", "CellType"])
        .agg(N_modules=("Module", "count"), Total_observations=("N_observations", "sum"))
        .reset_index()
    )
    print("\n=== Summary: observations and modules per (Cohort x Disorder x CellType) ===")
    print("(N_observations = total gene x cell pairs accumulated per group)")
    print(summary.to_string(index=False))
    if n_null:
        print(f"\nWARNING: {n_null} (cohort, disorder, celltype, module) groups had n<2 and SE=NaN.")
    print(f"\nOutput files: {output_files}")
    log.info("Done.")


if __name__ == "__main__":
    main()
