# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html

#conda create -n scvi-env python=3.12  # any python 3.10 to 3.12
#conda activate scvi-env
#pip install -U scvi-tools[tutorials]

import os
import tempfile

import scanpy as sc
import scvi
import seaborn as sns
import torch

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)


########## Tutorial
save_dir = tempfile.TemporaryDirectory()
adata = scvi.data.heart_cell_atlas_subsampled(save_path=save_dir.name)
adata
sc.pp.filter_genes(adata, min_counts=3)
adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # freeze the state in `.raw`
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="cell_source",
)
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["cell_source", "donor"],
    continuous_covariate_keys=["percent_mito", "percent_ribo"],
)

SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

denoised = model.get_normalized_expression(adata_subset, library_size=1e4)
denoised.iloc[:5, :5]
bla=adata.to_df(layer='counts')
bla.iloc[:5, :5]
df1 = denoised.set_index('row_names')
df2 = bla.set_index('row_names_column')

# 2. Reindex df1 to match the order of df2
df2 = bla.reindex(denoised.index)