import os
import json
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.io import mmread
from scipy.sparse import csc_matrix
from scipy.stats import sem

# Load dataset (SEAAD2024, all donors)
adata = ad.read_h5ad(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
sc.pp.log1p(adata, layer="UMIs", copy=False, chunked=False) 

# Load module list, set variables
# with open(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/topmodposbc_table.json"), "r") as f:
#     modules = json.load(f)

# mod_type = "topmodposbc"

with open(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/seed_table.json"), "r") as f:
    modules = json.load(f)

mod_type = "seed"

ct_string = "Subclass"
gene_string = "gene_ids"
case_control_names = ["Con", "AD"]
out_dir = os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/")

# Create column for case/control
adata.obs['case_control_col'] = 'NA'
adata.obs['case_control_col'][adata.obs['Overall AD neuropathological Change'] == 'Not AD'] = 'Con'
adata.obs['case_control_col'][adata.obs['Overall AD neuropathological Change'] != 'Not AD'] = 'AD'
adata.obs['case_control_col'][adata.obs['Overall AD neuropathological Change'] == 'Reference'] = 'Reference'

# For each subclass, calculate means over all genes
allmeanlist = []
for i in adata.obs[ct_string].unique():
    print(i)
    temp = adata[(adata.obs[ct_string] == i) & (adata.obs['case_control_col'].isin(case_control_names))]
    #temp2 = np.mean(temp.X, axis=0)
    temp2 = temp.X.mean(axis=0)
    allmeanlist.append(np.mean(temp2))

# For each module, calculate variance across each subclass
# v2 nested loop: mod first, then celltype
for c in case_control_names:
    varlist_native = []
    varlist_normMean = []
    varlist_REI = []
    for i in range(1, len(modules) + 1):
        print(i)
        current_timestamp = datetime.now()
        print(f"Current timestamp: {current_timestamp}")
        modvec = []
        modvec_normMean = []
        REIlist = []
        count = 0
        for j in adata.obs[ct_string].unique():
            if c != "":
                temp = adata[(adata.obs[ct_string] == j) & (adata.obs['case_control_col'] == c)]
            else:
                temp = adata[adata.obs[ct_string] == j]
            present_genes = list(set(modules[str(i)]) & set(temp.var[gene_string]))
            temp2 = temp[:, present_genes]
            tempx = temp2.X 
            tempx = tempx.reshape((1, -1))
            tempx = tempx.toarray()
            tempx_normMean = tempx / allmeanlist[count]
            modvec.append(sem(tempx.reshape(-1)))
            modvec_normMean.append(sem(tempx_normMean.reshape(-1)))
            REIlist.append(tempx_normMean)
            count = count + 1
        REImeans = [np.mean(sublist) for sublist in REIlist]
        REIlist_REI = [sublist / max(REImeans) for sublist in REIlist]
        modvec_REI = [sem(sublist.reshape(-1)) for sublist in REIlist_REI]
        varlist_native.append(modvec)
        varlist_normMean.append(modvec_normMean)
        varlist_REI.append(modvec_REI)
    flattened_data = np.vstack(varlist_native)
    combdf = pd.DataFrame(flattened_data)
    combdf.columns = adata.obs[ct_string].unique()
    combdf = combdf.sort_index(axis=1)
    filename = out_dir + '/sn_proj_indices/log_native/indices_se_' + ct_string + c + mod_type + '.csv'
    combdf.to_csv(filename, index = False)
    flattened_data = np.vstack(varlist_normMean)
    combdf = pd.DataFrame(flattened_data)
    combdf.columns = adata.obs[ct_string].unique()
    combdf = combdf.sort_index(axis=1)
    filename = out_dir + '/sn_proj_indices/log_normByMean/indices_se_' + ct_string + c + mod_type + '.csv'
    combdf.to_csv(filename, index = False)
    flattened_data = np.vstack(varlist_REI)
    combdf = pd.DataFrame(flattened_data)
    combdf.columns = adata.obs[ct_string].unique()
    combdf = combdf.sort_index(axis=1)
    filename = out_dir + '/sn_proj_indices/log_REI/indices_se_' + ct_string + c + mod_type + '.csv'
    combdf.to_csv(filename, index = False)

