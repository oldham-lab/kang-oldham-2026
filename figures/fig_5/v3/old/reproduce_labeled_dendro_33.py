# ============================================================
# Reproduce Dendrogram33
# Required input files (update paths as needed):
#   - mod_eig.csv                            : 24 x 101 eigengene matrix (rows=cell types, cols=modules)
#   - branchpoint_table_modeig_with_genes.csv: branchpoint table with Pct_of_Total and Pct_of_Total_genes
#   - indiv_leaf_modgene_per.csv             : per-module mod_per and gene_per values
# Output:
#   - dendrogram33.pdf
# ============================================================

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# ---- File paths ----
MOD_EIG_FILE   = 'mod_eig.csv'
BP_GENES_FILE  = 'branchpoint_table_modeig_with_genes.csv'
LEAF_PER_FILE  = 'indiv_leaf_modgene_per.csv'
OUTPUT_FILE    = 'dendrogram33.pdf'

# ---- 1. Load data ----
eig      = pd.read_csv(MOD_EIG_FILE, header=0, index_col=0)
modules  = list(eig.columns)
n        = len(modules)

# ---- 2. Compute Pearson correlation and complete-linkage clustering ----
corr_mat = eig.corr(method='pearson').values
corr_mat = (corr_mat + corr_mat.T) / 2
np.fill_diagonal(corr_mat, 1.0)
dist_mat = np.clip(1 - corr_mat, 0, None)
np.fill_diagonal(dist_mat, 0)
Z = linkage(squareform(dist_mat), method='complete')

cut_height = 0.3

# ---- 3. Identify branchpoints above cut ----
def get_leaves(node_id, Z, n):
    if node_id < n:
        return {int(node_id)}
    idx = int(node_id - n)
    return get_leaves(int(Z[idx, 0]), Z, n) | get_leaves(int(Z[idx, 1]), Z, n)

branchpoints = []
for i in range(len(Z)):
    if Z[i, 2] >= cut_height:
        node_id = n + i
        leaves  = get_leaves(node_id, Z, n)
        branchpoints.append({
            'node_id' : node_id,
            'height'  : Z[i, 2],
            'leaves'  : leaves,
            'n_leaves': len(leaves),
        })

branchpoints_sorted = sorted(branchpoints, key=lambda x: (x['height'], x['node_id']))
for j, bp in enumerate(branchpoints_sorted):
    bp['label'] = f"BP{j+1:02d}"

# ---- 4. Load annotation data ----
genes_df         = pd.read_csv(BP_GENES_FILE)
pct_mods_lookup  = dict(zip(genes_df['Branchpoint_ID'], genes_df['Pct_of_Total']))
pct_genes_lookup = dict(zip(genes_df['Branchpoint_ID'], genes_df['Pct_of_Total_genes']))
leaf_per         = pd.read_csv(LEAF_PER_FILE, index_col=0)

# ---- 5. Build figure layout ----
fig = plt.figure(figsize=(40, 24))
gs  = gridspec.GridSpec(4, 1, height_ratios=[1, 0.12, 0.12, 2],
                        hspace=0.0, top=0.90, bottom=0.18)
ax_dend   = fig.add_subplot(gs[0])
ax_modper = fig.add_subplot(gs[1])
ax_genper = fig.add_subplot(gs[2])
ax_heat   = fig.add_subplot(gs[3])

# ---- 6. Draw dendrogram ----
dend = dendrogram(
    Z, labels=modules, color_threshold=cut_height, ax=ax_dend,
    leaf_rotation=90, leaf_font_size=0, above_threshold_color='gray'
)
for collection in ax_dend.collections:
    collection.set_linewidth(3)
for line in ax_dend.lines:
    line.set_linewidth(3)

ax_dend.axhline(y=cut_height, color='red', linestyle='--', linewidth=2.5,
                label=f'Cut height = {cut_height}')
ax_dend.set_ylabel('Distance (1 − r)', fontsize=16)
ax_dend.legend(fontsize=13)
ax_dend.tick_params(axis='x', bottom=False, labelbottom=False)

# ---- 7. Compute node x-positions for branchpoint annotation ----
leaf_order = dend['leaves']
leaf_x     = {leaf: 10 * pos + 5 for pos, leaf in enumerate(leaf_order)}
node_x     = {i: leaf_x[i] for i in range(n)}
for i in range(len(Z)):
    node_id       = n + i
    node_x[node_id] = (node_x[int(Z[i, 0])] + node_x[int(Z[i, 1])]) / 2

# ---- 8. Annotate branchpoints ----
for bp in branchpoints_sorted:
    x         = node_x[bp['node_id']]
    y         = bp['height']
    pct_mods  = pct_mods_lookup.get(bp['label'], np.nan)
    pct_genes = pct_genes_lookup.get(bp['label'], np.nan)
    text      = f"{bp['label']}\n{pct_mods:.1f}% (m)\n{pct_genes:.1f}% (g)"
    ax_dend.text(x, y, text, fontsize=13, ha='center', va='bottom',
                 color=(0.55, 0.0, 0.0, 0.75), fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow',
                           alpha=0.5, edgecolor='none'),
                 zorder=3)

# ---- 9. Bar plot rows (Mod % and Gene %) ----
ordered_modules = [modules[i] for i in leaf_order]
ordered_modper  = [leaf_per.loc[m, 'mod_per']  for m in ordered_modules]
ordered_geneper = [leaf_per.loc[m, 'gene_per'] for m in ordered_modules]

for ax_bar, values, color in [
    (ax_modper, ordered_modper,  '#4878cf'),
    (ax_genper, ordered_geneper, '#6acc65'),
]:
    ax_bar.bar(range(n), values, width=1.0, color=color, linewidth=0)
    ax_bar.set_xlim(-0.5, n - 0.5)
    ax_bar.set_xticks([])
    ax_bar.set_yticks([])
    ax_bar.set_ylabel('')
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_bar.spines['bottom'].set_visible(False)
    for i, v in enumerate(values):
        ax_bar.text(i, v, f'{v:.1f}', ha='center', va='bottom',
                    fontsize=8, color='black')

# ---- 10. Heatmap ----
eig_ordered = eig[ordered_modules]
vmax = np.percentile(np.abs(eig_ordered.values), 98)
im   = ax_heat.imshow(
    eig_ordered.values, aspect='auto', cmap='RdBu_r',
    vmin=-vmax, vmax=vmax, interpolation='nearest'
)
ax_heat.set_xticks(range(n))
ax_heat.set_xticklabels(ordered_modules, rotation=90, fontsize=12)
ax_heat.set_ylabel('Cell type', fontsize=16)
ax_heat.set_yticks(range(eig.shape[0]))
ax_heat.set_yticklabels(list(eig.index), fontsize=13)

# ---- 11. Align all axes widths to heatmap ----
fig.canvas.draw()
heat_pos = ax_heat.get_position()
for ax in [ax_dend, ax_modper, ax_genper]:
    pos = ax.get_position()
    ax.set_position([heat_pos.x0, pos.y0, heat_pos.width, pos.height])
fig.canvas.draw()
heat_xlim = ax_heat.get_xlim()
ax_dend.set_xlim((heat_xlim[0] + 0.5) * 10, (heat_xlim[1] + 0.5) * 10)
ax_modper.set_xlim(heat_xlim)
ax_genper.set_xlim(heat_xlim)

# ---- 12. Mod % / Gene % axis labels ----
fig.canvas.draw()
axes_left_fig = ax_heat.get_position().x0
for ax_bar, label in [(ax_modper, 'Mod %'), (ax_genper, 'Gene %')]:
    bar_pos = ax_bar.get_position()
    fig.text(axes_left_fig - 0.005,
             bar_pos.y0 + bar_pos.height / 2,
             label, fontsize=8, ha='right', va='center',
             transform=fig.transFigure)

# ---- 13. Colorbar ----
heat_pos   = ax_heat.get_position()
cbar_width = heat_pos.width / 10
ax_cbar    = fig.add_axes([
    heat_pos.x0 + heat_pos.width / 2 - cbar_width / 2,
    0.06, cbar_width, 0.015
])
cbar = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
cbar.set_label('Eigengene value', fontsize=11)
ax_cbar.tick_params(labelsize=9)

# ---- 14. Title and save ----
fig.suptitle('Jorstad metamodules, consensusMin', fontsize=15, y=0.97)
plt.savefig(OUTPUT_FILE, bbox_inches='tight')
print(f"Saved {OUTPUT_FILE}")
