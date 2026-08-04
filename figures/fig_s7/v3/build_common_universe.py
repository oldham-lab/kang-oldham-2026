"""
build_common_universe.py — compute the common gene universe for Fig. S7 v3.

The seven plotted series draw genes from six raw sources. Each source's gene symbols
are canonicalized to current HGNC, then intersected across ALL sources. The result is
the exact gene set every series is restricted to before pooled per-celltype means are
computed (so all series are compared on an identical gene basis).

genes.txt (the analyzed/bulk gene list) is NOT used as a filter here — the universe is
the pure intersection of the datasets themselves.

Output: v3/common_genes.txt  (one canonical symbol per line, sorted; no header)
"""
import os
import logging
from pathlib import Path

from hgnc_common import load_alias_map, canonicalize

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

OUT = Path(__file__).resolve().parent / "common_genes.txt"

# Raw gene-symbol sources for every plotted series.
H5AD_SOURCES = {
    "jorstad_DFC": os.path.join(os.environ.get("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.h5ad"),
    "jorstad_MTG": os.path.join(os.environ.get("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_MTG.h5ad"),
    "seaad_A9":    os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"),
    "seaad_MTG":   os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"),
    "mit_ad":      os.path.join(os.environ.get("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025.h5ad"),
}
MORABITO_GENES = os.path.join(os.environ.get("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/morabito_2021/all_cells/genes.tsv")


def h5ad_var_names(path: str) -> list:
    """Read var_names (gene symbols) from an h5ad by reading only the var index
    dataset via h5py — avoids anndata's full obs/var load (seconds, not minutes)."""
    import h5py
    with h5py.File(path, "r") as f:
        var = f["var"]
        k = var.attrs.get("_index", "_index")
        if isinstance(k, bytes):
            k = k.decode()
        arr = var[k][:]
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


def main():
    alias = load_alias_map()
    log.info("Loaded HGNC alias map: %d prev->current entries.", len(alias))

    raw_sets = {}
    for name, path in H5AD_SOURCES.items():
        syms = h5ad_var_names(path)
        raw_sets[name] = syms
        log.info("%-12s : %6d genes (raw)", name, len(syms))

    import pandas as pd
    mora = pd.read_csv(MORABITO_GENES, header=None, names=["Gene"])["Gene"].astype(str).tolist()
    raw_sets["morabito"] = mora
    log.info("%-12s : %6d genes (raw)", "morabito", len(mora))

    # Canonicalize each source to current HGNC, dedupe within source.
    canon_sets = {}
    for name, syms in raw_sets.items():
        cset = set(canonicalize(syms, alias))
        canon_sets[name] = cset
        log.info("%-12s : %6d genes (canonical, unique)", name, len(cset))

    common = set.intersection(*canon_sets.values())
    log.info("=" * 60)
    log.info("Common universe (intersection of all %d sources): %d genes",
             len(canon_sets), len(common))

    common_sorted = sorted(common)
    OUT.write_text("\n".join(common_sorted) + "\n")
    log.info("Wrote %s (%d genes)", OUT, len(common_sorted))


if __name__ == "__main__":
    main()
