import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

df = pd.read_csv('/mnt/user-data/uploads/test.csv', header=0, index_col=None)
labels = list(df.columns)
n = len(labels)

corr = df.values.astype(float)
corr = (corr + corr.T) / 2
np.fill_diagonal(corr, 1.0)
dist = np.clip(1 - corr, 0, None)
np.fill_diagonal(dist, 0)
condensed = squareform(dist)
Z = linkage(condensed, method='complete')

cut_height = 0.3

def get_leaves(node_id, Z, n):
    if node_id < n:
        return {int(node_id)}
    idx = int(node_id - n)
    return get_leaves(int(Z[idx, 0]), Z, n) | get_leaves(int(Z[idx, 1]), Z, n)

branchpoints = []
for i in range(len(Z)):
    if Z[i, 2] >= cut_height:
        node_id = n + i
        leaves = get_leaves(node_id, Z, n)
        branchpoints.append({
            'node_id': node_id,
            'height': Z[i, 2],
            'leaves': leaves,
            'n_leaves': len(leaves),
            'pct': len(leaves) / n * 100
        })

branchpoints_sorted = sorted(branchpoints, key=lambda x: (x['height'], x['node_id']))
for j, bp in enumerate(branchpoints_sorted):
    bp['label'] = f"BP{j+1:02d}"

genes_df = pd.read_csv('/mnt/user-data/uploads/branchpoint_table3_with_genes.csv')
genes_lookup = dict(zip(genes_df['Branchpoint_ID'], genes_df['Pct_of_Total_genes']))

eig = pd.read_csv('/mnt/user-data/uploads/mod_eig1.csv', header=0, index_col=None)

fig = plt.figure(figsize=(40, 22))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2], hspace=0.02,
                       top=0.90, bottom=0.18)  # leave room at bottom for colorbar and top for title
ax_dend = fig.add_subplot(gs[0])
ax_heat = fig.add_subplot(gs[1])

# Draw dendrogram
dend = dendrogram(
    Z, labels=labels, color_threshold=cut_height, ax=ax_dend,
    leaf_rotation=90, leaf_font_size=0, above_threshold_color='gray'
)

for collection in ax_dend.collections:
    collection.set_linewidth(3)
for line in ax_dend.lines:
    line.set_linewidth(3)

ax_dend.axhline(y=cut_height, color='red', linestyle='--', linewidth=2.5, label=f'Cut height = {cut_height}')
ax_dend.set_ylabel('Distance (1 − r)', fontsize=13)
ax_dend.legend(fontsize=11)
ax_dend.tick_params(axis='x', bottom=False, labelbottom=False)

leaf_order = dend['leaves']
leaf_x = {leaf: 10 * pos + 5 for pos, leaf in enumerate(leaf_order)}
node_x = {}
for i in range(n):
    node_x[i] = leaf_x[i]
for i in range(len(Z)):
    node_id = n + i
    node_x[node_id] = (node_x[int(Z[i, 0])] + node_x[int(Z[i, 1])]) / 2

font_size = 13
for bp in branchpoints_sorted:
    pct_genes = genes_lookup.get(bp['label'], np.nan)
    x = node_x[bp['node_id']]
    y = bp['height']
    text = f"{bp['label']}\n{bp['pct']:.1f}% (m)\n{pct_genes:.1f}% (g)"
    ax_dend.text(x, y, text, fontsize=font_size, ha='center', va='bottom',
                 color=(0.55, 0.0, 0.0, 0.75), fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow', alpha=0.5, edgecolor='none'),
                 zorder=3)

ordered_labels = [labels[i] for i in leaf_order]
eig_ordered = eig[ordered_labels]

vmax = np.percentile(np.abs(eig_ordered.values), 98)
im = ax_heat.imshow(
    eig_ordered.values, aspect='auto', cmap='RdBu_r',
    vmin=-vmax, vmax=vmax, interpolation='nearest'
)

ax_heat.set_xticks(range(n))
ax_heat.set_xticklabels(ordered_labels, rotation=90, fontsize=8)
ax_heat.set_ylabel('Sample', fontsize=13)
ax_heat.set_yticks(range(eig.shape[0]))
ax_heat.set_yticklabels([f'S{i+1}' for i in range(eig.shape[0])], fontsize=9)

# Align dendrogram axes width to heatmap
fig.canvas.draw()
heat_pos = ax_heat.get_position()
dend_pos = ax_dend.get_position()
ax_dend.set_position([heat_pos.x0, dend_pos.y0, heat_pos.width, dend_pos.height])
fig.canvas.draw()

heat_xlim = ax_heat.get_xlim()
ax_dend.set_xlim((heat_xlim[0] + 0.5) * 10, (heat_xlim[1] + 0.5) * 10)

# Small horizontal colorbar: centered, narrow, well below the heatmap
# Place it in figure coordinates: centered horizontally, 1/10th width of heatmap
heat_pos = ax_heat.get_position()
cbar_width = heat_pos.width / 10
cbar_height = 0.015
cbar_x = heat_pos.x0 + heat_pos.width / 2 - cbar_width / 2  # centered
cbar_y = 0.06  # well below the heatmap x-tick labels

ax_cbar = fig.add_axes([cbar_x, cbar_y, cbar_width, cbar_height])
cbar = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
cbar.set_label('Eigengene value', fontsize=11)
ax_cbar.tick_params(labelsize=9)

fig.suptitle('Jorstad metamodules, consensusMean', fontsize=15, y=0.97)

plt.savefig('/home/claude/dendrogram20.png', dpi=150, bbox_inches='tight')
print("Done.")

