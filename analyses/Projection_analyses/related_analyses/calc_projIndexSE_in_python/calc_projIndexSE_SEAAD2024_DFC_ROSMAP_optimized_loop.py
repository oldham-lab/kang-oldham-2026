import os
import json
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.stats import sem

# ─── Load dataset once ───────────────────────────────────────────────────────
adata = ad.read_h5ad(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
adata.X = adata.layers["UMIs"].copy()
sc.pp.log1p(adata)

with open(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.json"), "r") as f:
    modules = json.load(f)

mod_type = "topmodposbc"
ct_string = "Subclass"
gene_string = "gene_ids"

# ─── Define runs ─────────────────────────────────────────────────────────────
# Each dict specifies one run:
#   case_control_names : the two group labels to compare
#   out_dir            : where to write results
#   cc_col_map         : dict mapping adata.obs column values -> group label (or None to skip)
runs = [
    {
        "case_control_names": ["Con", "AD"],
        "out_dir": os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/"),
        "cc_col_map": {
            "source_col": "Overall AD neuropathological Change",
            "mapping": {
                "Not AD": "Con",
                # anything else -> "AD", handled by default below
            },
            "default": "AD",
            "exclude": ["Reference"],  # set these to NA / exclude from analysis
        }
    },
    {
        "case_control_names": ["Con", "Early"],
        "out_dir": os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_earlyVsCon_DFC/"),
        "cc_col_map": {
            "source_col": "Overall AD neuropathological Change",
            "mapping": {
                "Not AD": "Con",
                "Low": "Early",
            },
            "default": None,   # None means rows not in mapping or exclude get set to NA
            "exclude": ["Reference"]
        }
    },
    {
        "case_control_names": ["Late", "Early"],
        "out_dir": os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_lateVsEarly_DFC/"),
        "cc_col_map": {
            "source_col": "Overall AD neuropathological Change",
            "mapping": {
                "High": "Late",
                "Low": "Early",
            },
            "default": None,   # None means rows not in mapping or exclude get set to NA
            "exclude": ["Reference", "Not AD"]
        }
    },
]

# ─── One-time gene index precomputation ──────────────────────────────────────
gene_to_col = {g: i for i, g in enumerate(adata.var[gene_string])}

print("Precomputing module gene indices...")
module_col_indices = {
    str(i): [gene_to_col[g] for g in modules[str(i)] if g in gene_to_col]
    for i in range(1, len(modules) + 1)
}

celltypes = adata.obs[ct_string].unique()

# ─── Run loop ────────────────────────────────────────────────────────────────
for run in runs:
    case_control_names = run["case_control_names"]
    out_dir = run["out_dir"]
    cc_map = run["cc_col_map"]

    print(f"\n{'='*60}")
    print(f"Starting run: {case_control_names} -> {out_dir}")
    print(f"{'='*60}")

    # Build case_control_col for this run
    source_col = cc_map["source_col"]
    adata.obs["case_control_col"] = "NA"
    for obs_val, label in cc_map["mapping"].items():
        adata.obs.loc[adata.obs[source_col] == obs_val, "case_control_col"] = label
    if cc_map.get("default"):
        # Apply default to any rows not yet assigned and not excluded
        unassigned = adata.obs["case_control_col"] == "NA"
        not_excluded = ~adata.obs[source_col].isin(cc_map.get("exclude", []))
        adata.obs.loc[unassigned & not_excluded, "case_control_col"] = cc_map["default"]
    for exc_val in cc_map.get("exclude", []):
        adata.obs.loc[adata.obs[source_col] == exc_val, "case_control_col"] = "NA"

    # Precompute celltype means for this run's case/control grouping
    print("Precomputing celltype means...")
    allmeanlist = []
    for ct in celltypes:
        temp = adata[(adata.obs[ct_string] == ct) & (adata.obs["case_control_col"].isin(case_control_names))]
        allmeanlist.append(float(temp.X.mean()))

    # Pre-slice dense arrays for this run
    print("Pre-slicing data per celltype/condition...")
    ct_cc_arrays = {}
    for c in case_control_names:
        for ct in celltypes:
            mask = (adata.obs[ct_string] == ct) & (adata.obs["case_control_col"] == c)
            ct_cc_arrays[(ct, c)] = adata[mask].X.toarray()
            print(f"  Sliced ({ct}, {c}) | {datetime.now()}")

    # Main computation loop
    for c in case_control_names:
        varlist_native = []
        varlist_normMean = []
        varlist_REI = []

        for i in range(1, len(modules) + 1):
            if i % 20 == 0:
                print(f"  Module {i}/{len(modules)} | {datetime.now()}")
            col_idx = module_col_indices[str(i)]
            modvec = []
            modvec_normMean = []
            REIlist = []

            for ct_idx, ct in enumerate(celltypes):
                arr = ct_cc_arrays[(ct, c)][:, col_idx]
                flat = arr.ravel()
                flat_normMean = flat / allmeanlist[ct_idx]
                modvec.append(sem(flat))
                modvec_normMean.append(sem(flat_normMean))
                REIlist.append(flat_normMean)

            REImeans = np.array([np.mean(f) for f in REIlist])
            modvec_REI = [sem(f / REImeans.max()) for f in REIlist]

            varlist_native.append(modvec)
            varlist_normMean.append(modvec_normMean)
            varlist_REI.append(modvec_REI)

        for varlist, subfolder in [
            (varlist_native, "log_native"),
            (varlist_normMean, "log_normByMean"),
            (varlist_REI, "log_REI"),
        ]:
            combdf = pd.DataFrame(np.vstack(varlist), columns=celltypes).sort_index(axis=1)
            fname = f"{out_dir}/sn_proj_indices/{subfolder}/indices_se_{ct_string}{c}{mod_type}.csv"
            combdf.to_csv(fname, index=False)
            print(f"  Saved: {fname}")

    # Free memory before next run
    del ct_cc_arrays

print("\nAll runs complete.")