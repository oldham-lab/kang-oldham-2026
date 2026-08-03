import os
import json
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.stats import sem

# Load dataset
adata = ad.read_h5ad(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
sc.pp.log1p(adata, layer="UMIs", copy=False, chunked=False)

# with open(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/seed_table.json"), "r") as f:
#     modules = json.load(f)

# mod_type = "seed"

with open(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/topmodposbc_table.json"), "r") as f:
    modules = json.load(f)

mod_type = "topmodposbc"
ct_string = "Subclass"
gene_string = "gene_ids"
case_control_names = ["Con", "AD"]
out_dir = os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/")

# Fixed: use .loc to avoid chained assignment warnings/failures
adata.obs['case_control_col'] = 'NA'
adata.obs.loc[adata.obs['Overall AD neuropathological Change'] == 'Not AD', 'case_control_col'] = 'Con'
adata.obs.loc[adata.obs['Overall AD neuropathological Change'] != 'Not AD', 'case_control_col'] = 'AD'
adata.obs.loc[adata.obs['Overall AD neuropathological Change'] == 'Reference', 'case_control_col'] = 'Reference'

celltypes = adata.obs[ct_string].unique()
all_var_genes = set(adata.var[gene_string])

# Precompute global celltype means (unchanged logic)
print("Precomputing celltype means...")
allmeanlist = []
for ct in celltypes:
    temp = adata[(adata.obs[ct_string] == ct) & (adata.obs['case_control_col'].isin(case_control_names))]
    allmeanlist.append(float(temp.X.mean()))

# Build gene -> integer column index map for fast lookup
gene_to_col = {g: i for i, g in enumerate(adata.var[gene_string])}

# Precompute module gene column indices once (not inside the loop)
print("Precomputing module gene indices...")
module_col_indices = {
    str(i): [gene_to_col[g] for g in modules[str(i)] if g in gene_to_col]
    for i in range(1, len(modules) + 1)
}

# KEY SPEEDUP: Pre-extract dense arrays per (celltype, condition) once.
# All module iterations then do cheap numpy column slicing instead of sparse adata subsetting.
print("Pre-slicing data per celltype/condition...")
ct_cc_arrays = {}  # (ct, condition) -> dense numpy array shape (n_cells, n_genes)
for c in case_control_names:
    for ct in celltypes:
        mask = (adata.obs[ct_string] == ct) & (adata.obs['case_control_col'] == c)
        ct_cc_arrays[(ct, c)] = adata[mask].X.toarray()  # pay this cost once
        current_timestamp = datetime.now()
        print(f"Current timestamp: {current_timestamp}")

for c in case_control_names:
    varlist_native = []
    varlist_normMean = []
    varlist_REI = []
    for i in range(1, len(modules) + 1):
        if i % 20 == 0:
            print(f"Module {i}/{len(modules)} | {datetime.now()}")
        col_idx = module_col_indices[str(i)]
        modvec = []
        modvec_normMean = []
        REIlist = []
        for ct_idx, ct in enumerate(celltypes):
            arr = ct_cc_arrays[(ct, c)][:, col_idx]  # (n_cells, n_mod_genes)
            flat = arr.ravel()
            flat_normMean = flat / allmeanlist[ct_idx]
            modvec.append(sem(flat))
            modvec_normMean.append(sem(flat_normMean))
            REIlist.append(flat_normMean)  # REI builds on normMean, matching original
        REImeans = np.array([np.mean(f) for f in REIlist])
        max_REImean = REImeans.max()
        modvec_REI = [sem(f / max_REImean) for f in REIlist]
        varlist_native.append(modvec)
        varlist_normMean.append(modvec_normMean)
        varlist_REI.append(modvec_REI)
    for varlist, subfolder in [
        (varlist_native, 'log_native'),
        (varlist_normMean, 'log_normByMean'),
        (varlist_REI, 'log_REI'),
    ]:
        combdf = pd.DataFrame(np.vstack(varlist), columns=celltypes).sort_index(axis=1)
        fname = f"{out_dir}/sn_proj_indices/{subfolder}/indices_se_{ct_string}{c}{mod_type}.csv"
        combdf.to_csv(fname, index=False)