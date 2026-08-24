"""
hgnc_common.py — shared HGNC symbol canonicalization + common-universe selection
for the Fig. S7 v3 pipeline.

Gene symbols across the seven snRNA-seq series come from different HGNC vintages
(e.g. Jorstad uses current symbols, SEA-AD an older vintage). Canonicalizing every
dataset's symbols to current HGNC before intersecting recovers genes that would
otherwise be dropped as spurious mismatches. This mirrors the logic in
fig_1/full_DE_pipeline_optimized_v2.4.R (canonicalize_symbols).
"""
import os
from pathlib import Path

import pandas as pd

# prev_symbol -> current_symbol map (HGNC complete set), same file the DE pipeline uses.
ALIAS_MAP_PATH = os.environ.get("HGNC_MAP", "/home/gugene/code/git/CoPA/inst/python/dataset_processing_unified/hgnc_prev_to_current.tsv")


def load_alias_map(path: str = ALIAS_MAP_PATH) -> dict:
    """Load the prev_symbol -> current_symbol lookup (unique prev symbols)."""
    t = pd.read_csv(path, sep="\t")
    return dict(zip(t["prev_symbol"].astype(str), t["current_symbol"].astype(str)))


def canonicalize(symbols, alias: dict = None) -> list:
    """Map symbols to current HGNC; leave unmapped symbols unchanged."""
    if alias is None:
        alias = load_alias_map()
    return [alias.get(s, s) for s in symbols]


def select_common_indices(dataset_symbols, common_genes, alias: dict = None):
    """
    Given a dataset's raw gene symbols and the canonical common-universe gene list,
    return (gene_idx, present_genes):

      gene_idx      : for each present common gene, the *first* dataset column index
                      whose canonical symbol equals that gene (dedupes collapses so no
                      gene is counted twice).
      present_genes : common genes actually found in this dataset, in common_genes order.

    Because common_genes is the intersection across all datasets, present_genes should
    equal common_genes for every dataset (so N_genes is identical everywhere).
    """
    if alias is None:
        alias = load_alias_map()
    common_set = set(common_genes)
    first_idx = {}
    for i, s in enumerate(dataset_symbols):
        c = alias.get(s, s)
        if c in common_set and c not in first_idx:
            first_idx[c] = i
    present_genes = [g for g in common_genes if g in first_idx]
    gene_idx = [first_idx[g] for g in present_genes]
    return gene_idx, present_genes
