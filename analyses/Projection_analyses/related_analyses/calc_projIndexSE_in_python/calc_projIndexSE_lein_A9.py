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

# # Load dataset (Lein DFC) from mtx and convert to anndata object with genes and barcodes
# scipy_csc_matrix = mmread("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x//expression_DFC/matrix.mtx")
# if not isinstance(scipy_csc_matrix, csc_matrix):
#     scipy_csc_matrix = scipy_csc_matrix.tocsc()
# adata = ad.AnnData(scipy_csc_matrix.T)
# genetsv = pd.read_csv("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC/genes.tsv", sep='\t', header=None)
# genetsv.columns = ['gene_names']
# genetsv = genetsv.set_index('gene_names', drop=False)
# sn_anno = pd.read_csv(os.path.join(os.environ.get("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"))
# sn_anno = sn_anno.set_index('Cell_ID', drop=False)
# adata.obs = sn_anno
# adata.var = genetsv
# adata.raw = adata
# sc.pp.log1p(adata)
# adata.write(filename="/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x//expression_DFC.h5ad")

# # This method doesn't work for some reason ("Keyerror: 1")
# adata = sc.read_10x_mtx("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC/")
# adata.write(filename = "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.h5ad")

# Read data (Lein DFC)
adata = ad.read_h5ad(os.path.join(os.environ.get("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.h5ad"))

# Load module list

# # Topmodposbc
# with open(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/topmodposbc_table.json"), "r") as f:
#     modules = json.load(f)
#
# out_dir = os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/")
# mod_type = "topmodposbc"

# Seed
with open(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/seed_table.json"), "r") as f:
    modules = json.load(f)

out_dir = os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/")
mod_type = "seed"

# For each subclass, calculate means over all genes
allmeanlist = []
for i in adata.obs["Cell_Type"].unique():
    print(i)
    temp = adata[adata.obs["Cell_Type"] == i]
    temp = temp.X.toarray()
    allmeanlist.append(np.mean(np.mean(temp, axis = 0)))

# # For each module, calculate variance across each subclass
# # v1 nested loop: cell type first, then mod. this makes calculating REI awkward.
# varlist_native = []
# varlist_normMean = []
# varlist_REI = []
# count = 0
# for i in adata.obs["Cell_Type"].unique():
#     print(i)
#     current_timestamp = datetime.now()
#     print(f"Current timestamp: {current_timestamp}")
#     temp = adata[adata.obs["Cell_Type"] == i]
#     modvec = []
#     modvec_normMean = []
#     REIlist = []
#     for j in range(1, len(modules) + 1):
#         present_genes = list(set(modules[str(j)]) & set(temp.var['gene_names']))
#         temp2 = temp[:, present_genes]
#         temp2 = temp2.X.toarray()
#         temp2 = temp2.flatten()
#         temp2_normMean = temp2 / allmeanlist[count]
#         modvec.append(sem(temp2))
#         modvec_normMean.append(sem(temp2_normMean))
#         REIlist.append(temp2_normMean)

#     varlist_native.append(modvec)
#     varlist_normMean.append(modvec_normMean)
#     count = count + 1

# v2 nested loop: mod first, then celltype
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
    for j in adata.obs["Cell_Type"].unique():
        temp = adata[adata.obs["Cell_Type"] == j]
        present_genes = list(set(modules[str(i)]) & set(temp.var['gene_names']))
        temp2 = temp[:, present_genes]
        # temp2 = temp2.X.toarray()
        # temp2 = temp2.flatten()
        # temp2_normMean = temp2 / allmeanlist[count]
        # modvec.append(np.var(temp2))
        # modvec_normMean.append(np.var(temp2_normMean))
        # REIlist.append(temp2_normMean)
        tempx = temp2.X / allmeanlist[count]
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
combdf.columns = adata.obs["Cell_Type"].unique()
combdf = combdf.sort_index(axis=1)
filename = out_dir + '/sn_proj_indices/log_native/indices_se_' + mod_type + '.csv'
combdf.to_csv(filename, index = False)
#combdf.to_csv(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_proj_indices/log_native/indices_se_Cell_Type.csv"), index = False)


flattened_data = np.vstack(varlist_normMean)
combdf = pd.DataFrame(flattened_data)
combdf.columns = adata.obs["Cell_Type"].unique()
combdf = combdf.sort_index(axis=1)
filename = out_dir + '/sn_proj_indices/log_normByMean/indices_se_' + mod_type + '.csv'
combdf.to_csv(filename, index = False)
#combdf.to_csv(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_proj_indices/log_normByMean/indices_se_Cell_Type.csv"), index = False)


flattened_data = np.vstack(varlist_REI)
combdf = pd.DataFrame(flattened_data)
combdf.columns = adata.obs["Cell_Type"].unique()
combdf = combdf.sort_index(axis=1)
filename = out_dir + '/sn_proj_indices/log_REI/indices_se_' + mod_type + '.csv'
combdf.to_csv(filename, index = False)
#combdf.to_csv(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_proj_indices/log_REI/indices_se_Cell_Type.csv"), index = False)

