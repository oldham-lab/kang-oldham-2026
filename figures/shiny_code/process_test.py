# Single-cell multiregion epigenomic rewiring in Alzheimer’s disease progression and cognitive resilience (Cell 2025)
# DOI: 10.1016/j.cell.2025.06.031 
# https://compbio.mit.edu/AD_Multiomic_MultiRegion/

import json
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy import io
print(ad.__version__)

os.chdir('/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/MTC/')
adata = ad.read_h5ad('/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025.h5ad', backed='r')

# Load module list, set variables
with open(os.path.join(os.environ.get("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/topmodposbc_table.json"), "r") as f:
    modules = json.load(f)

ct_string = "Subclass"
gene_string = "gene_ids"
case_control_names = ["Con", "AD"]

# Testing code
testlist = []
for i in adata.obs["RNA.Subclass"].unique():
    print(i)
    temp = adata[(adata.obs["RNA.Subclass"] == i) & (adata.obs["Pathology"] != 'nonAD') & (adata.obs["BrainRegion"] == "MTC")].to_memory()
    for j in modules:
        #temp = sc.pp.log1p(adata[(adata.obs["RNA.Subclass"] == i) & (adata.obs["Pathology"] != 'nonAD') & (adata.obs["BrainRegion"] == "MTC"), j].to_memory()) # slow
        temp2 = temp[:,j].copy()
        sc.pp.log1p(temp2)
        temp3 = temp2.X.toarray()
        sc.tl.rank_genes_groups(temp2, groupby="RNA.Subclass", method="wilcoxon")
        os.chdir("/home/gugene/test")
        sc.pl.violin(temp2, n_genes = 10, jitter = False, save = "test.png")

        blatest = sns.load_dataset("titanic")
        bla = pd.DataFrame(temp2.X.toarray(), index=temp2.obs['RNA.Subclass'], columns=temp2.var_names)
        bla2 = bla.stack().reset_index()
        sns.catplot(data=bla2, x="RNA.Subclass", y=0, kind="violin")
        plt.savefig('my_seaborn_plot.pdf')

        #temp2 = np.log(temp + 1)
        testlist.append(np.mean(temp2, axis = 0))

# Make a plot per module
#for j in modules:
temp = adata[(adata.obs["Pathology"] != 'nonAD') & (adata.obs["BrainRegion"] == "MTC")].to_memory() # slow
temp2 = temp[:,j].copy()
sc.pp.log1p(temp2)

os.chdir(os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/python_exploratory/"))

row_multiindex_tuples = list(zip(temp2.obs['RNA.Subclass'], temp2.obs['Pathology']))
row_multiindex = pd.MultiIndex.from_tuples(row_multiindex_tuples, names=['Subclass', 'Dx'])
bla = pd.DataFrame(temp2.X.toarray(), index=row_multiindex, columns=temp2.var_names)
bla2 = bla.stack().reset_index()
sns.catplot(data=bla2, x="Subclass", hue="Dx", y=0, kind="violin", aspect = 5)
plt.savefig('module1_catplot_violin.pdf')
sns.catplot(data=bla2, x="Subclass", hue="Dx", y=0, kind="violin", aspect = 5, split="True")
plt.savefig('module1_catplot_violin_split.pdf')
sns.catplot(data=bla2, x="Subclass", hue="Dx", y=0, kind="box", aspect = 5)
plt.savefig('module1_catplot_box.pdf')
sns.catplot(data=bla2, x="Subclass", hue="Dx", y=0, kind="boxen", aspect = 5)
plt.savefig('module1_catplot_boxen.pdf')