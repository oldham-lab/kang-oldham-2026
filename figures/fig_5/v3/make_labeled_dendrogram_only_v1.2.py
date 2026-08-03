# added argparse through cli

# ============================================================
# Dendrogram-only plot based on Dendrogram33
# Changes from original:
#   - Only the dendrogram panel (no heatmap, no bar plots)
#   - Leaves below the cut line are all grey (no cluster coloring)
#   - Y-axis title is "1 - cor"
#   - Refreshed label aesthetics: dark slate background boxes,
#     white text, monospaced font, minimal connector lines
# Required inputs:
#   - mod_eig.csv                            : 24 x 101 eigengene matrix
#   - branchpoint_table_modeig_with_genes.csv: branchpoint annotations
# Output:
#   - dendrogram_only.pdf
# ============================================================

# Examples of input:
# # Basic — positional args required
# python make_labeled_dendrogram_only.py mod_eig.csv branchpoint_table_modeig_with_genes.csv

# # With custom output name
# python make_labeled_dendrogram_only.py mod_eig.csv branchpoint_table_modeig_with_genes.csv -o my_output.pdf

# # With optional tweaks
# python make_labeled_dendrogram_only.py mod_eig.csv branchpoint_table_modeig_with_genes.csv \
#     -o my_output.pdf --cut-height 0.25 --font-size 18

# # Help
# python make_labeled_dendrogram_only.py --help

import argparse
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.font_manager as fm
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# ---- Command-line arguments ----
parser = argparse.ArgumentParser(
    description='Produce a labeled dendrogram from a module eigengene matrix.'
)
parser.add_argument(
    'mod_eig_file',
    help='Path to mod_eig CSV (rows=cell types, cols=modules)'
)
parser.add_argument(
    'bp_genes_file',
    help='Path to branchpoint annotation CSV (must contain Branchpoint_ID, '
         'Pct_of_Total, Pct_of_Total_genes columns)'
)
parser.add_argument(
    '-o', '--output',
    default='dendrogram_only.pdf',
    help='Output file path (default: dendrogram_only.pdf). '
         'Both a PDF and SVG will be saved with the same stem.'
)
parser.add_argument(
    '--cut-height',
    type=float,
    default=0.3,
    help='Cut height for branchpoint labeling (default: 0.3)'
)
parser.add_argument(
    '--font-size',
    type=float,
    default=20.23,
    help='Branchpoint label font size in points (default: 20.23)'
)
args = parser.parse_args()

# ---- Font setup: use Arial if available, otherwise closest alternatives ----
def find_font(candidates):
    """Return the name of the first font from candidates that matplotlib can find."""
    available = {f.name for f in fm.fontManager.ttflist}
    for name in candidates:
        if name in available:
            return name
    return candidates[-1]  # final fallback

FONT_FAMILY = find_font([
    'Arial',           # Windows / macOS (with Office)
    'Helvetica',       # macOS system font
    'Nimbus Sans',     # Linux (Arial substitute in many distros)
    'Liberation Sans', # Linux open-source Arial metric-compatible font
    'FreeSans',        # another Linux fallback
    'DejaVu Sans',     # matplotlib's built-in fallback
])
print(f"Using font: {FONT_FAMILY}")

# Apply globally so all text elements use it
matplotlib.rcParams['font.family'] = FONT_FAMILY

# ---- File paths (from command line) ----
MOD_EIG_FILE  = args.mod_eig_file
BP_GENES_FILE = args.bp_genes_file
OUTPUT_FILE   = args.output

# ---- 1. Load data and cluster ----
eig     = pd.read_csv(MOD_EIG_FILE, header=0, index_col=0)
modules = list(eig.columns)
n       = len(modules)

corr_mat = eig.corr(method='pearson').values
corr_mat = (corr_mat + corr_mat.T) / 2
np.fill_diagonal(corr_mat, 1.0)
dist_mat = np.clip(1 - corr_mat, 0, None)
np.fill_diagonal(dist_mat, 0)
Z = linkage(squareform(dist_mat), method='complete')

cut_height = args.cut_height

# ---- 2. Identify branchpoints above cut ----
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

genes_df         = pd.read_csv(BP_GENES_FILE)
pct_mods_lookup  = dict(zip(genes_df['Branchpoint_ID'], genes_df['Pct_of_Total']))
pct_genes_lookup = dict(zip(genes_df['Branchpoint_ID'], genes_df['Pct_of_Total_genes']))

# ---- 3. Draw dendrogram ----
fig, ax = plt.subplots(figsize=(36, 10))
fig.patch.set_facecolor('#FAFAFA')
ax.set_facecolor('#FAFAFA')

dend = dendrogram(
    Z,
    labels=modules,
    color_threshold=0,
    above_threshold_color='#4A4A4A',
    ax=ax,
    leaf_rotation=90,
    leaf_font_size=0,           # hide leaf labels
    link_color_func=lambda k: '#4A4A4A',
)

# Thicken lines
for line in ax.get_lines():
    line.set_linewidth(3.5)
    line.set_color('#4A4A4A')

# Cut line
ax.axhline(y=cut_height, color='#E63946', linestyle='--', linewidth=2.0,
           label=f'Cut height = {cut_height}', zorder=5)

ax.set_ylabel('1 - cor', fontsize=28, labelpad=14, fontfamily=FONT_FAMILY)
ax.set_xlabel('')
ax.tick_params(axis='y', labelsize=22)
ax.tick_params(axis='x', bottom=False, labelbottom=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.legend(fontsize=20, framealpha=0.6, prop={'family': FONT_FAMILY, 'size': 20})
# No title

# ---- 4. Compute node x-positions ----
leaf_order = dend['leaves']
leaf_x     = {leaf: 10 * pos + 5 for pos, leaf in enumerate(leaf_order)}
node_x     = {i: leaf_x[i] for i in range(n)}
for i in range(len(Z)):
    node_id         = n + i
    node_x[node_id] = (node_x[int(Z[i, 0])] + node_x[int(Z[i, 1])]) / 2

# ---- 5. Annotate branchpoints — dark slate style ----
BOX_COLOR  = '#1D3557'
TEXT_COLOR = '#1A1A2E'      # dark charcoal
LINE_COLOR = '#457B9D'

ax_xlim = ax.get_xlim()
ax_ylim = ax.get_ylim()
fig.canvas.draw()
renderer  = fig.canvas.get_renderer()

fig_w, fig_h = fig.get_size_inches()
ax_bbox = ax.get_position()
x_range = ax_xlim[1] - ax_xlim[0]
y_range = ax_ylim[1] - ax_ylim[0]
ax_w_in = fig_w * ax_bbox.width
ax_h_in = fig_h * ax_bbox.height

# Labels start at their branchpoint. Tiny nudges (max 0.5 * label_h upward,
# max 0.5 * label_w horizontal) are allowed only to resolve overlaps.
# Binary-search for the largest font size where nudged positions don't overlap.

MAX_NUDGE_Y = 1.5   # in units of label_h
MAX_NUDGE_X = 0.5   # in units of label_w

def label_dims(fs):
    lh = (fs / 72.0 * 1.5) * 3 * (y_range / ax_h_in)
    lw = (fs / 72.0 * 0.62) * 10 * (x_range / ax_w_in)
    return lw, lh

def boxes_overlap(ax1, ay1, ax2, ay2, lw, lh, margin=0.0):
    return (abs(ax1 - ax2) < lw * (1 + margin) and
            abs((ay1 + lh/2) - (ay2 + lh/2)) < lh * (1 + margin))

def nudge_positions(anchors, lw, lh):
    """Resolve overlaps with minimal nudges. Returns list of (x, y) placed positions."""
    # Start at anchor positions
    pos = list(anchors)
    # Iterate: for each overlapping pair nudge the higher-indexed one slightly
    for _ in range(200):
        moved = False
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                xi, yi = pos[i]
                xj, yj = pos[j]
                if boxes_overlap(xi, yi, xj, yj, lw, lh):
                    # Nudge j: try small upward step first, then horizontal
                    step_y = lh * 0.12
                    step_x = lw * 0.12
                    ax_i, ay_i = anchors[i]
                    ax_j, ay_j = anchors[j]
                    best = None
                    best_dist = float('inf')
                    # Try candidate nudges for j
                    for dy in np.arange(0, MAX_NUDGE_Y * lh + step_y, step_y):
                        for dx in np.arange(-MAX_NUDGE_X * lw,
                                             MAX_NUDGE_X * lw + step_x, step_x):
                            cx, cy = ax_j + dx, ay_j + dy
                            # Check against all already-placed labels
                            ok = all(
                                not boxes_overlap(cx, cy, pos[k][0], pos[k][1], lw, lh)
                                for k in range(len(pos)) if k != j
                            )
                            d = np.hypot(cx - ax_j, cy - ay_j)
                            if ok and d < best_dist:
                                best_dist = d
                                best = (cx, cy)
                    if best is not None:
                        pos[j] = best
                        moved = True
                        break
            if moved:
                break
        if not moved:
            break
    return pos

def count_overlaps(pos, lw, lh):
    n = 0
    for i in range(len(pos)):
        for j in range(i + 1, len(pos)):
            if boxes_overlap(pos[i][0], pos[i][1], pos[j][0], pos[j][1], lw, lh):
                n += 1
    return n

# Anchor each label at its branchpoint
anchors = [(node_x[bp['node_id']], bp['height']) for bp in branchpoints_sorted]

# Fixed font size (binary-search optimum + 2pt)
font_size = args.font_size
lw, lh = label_dims(font_size)
best_pos = nudge_positions(anchors, lw, lh)
best_fs  = font_size
label_w, label_h = label_dims(font_size)
print(f"Optimal font size: {font_size:.2f}pt")

# Draw labels; draw a short connector line only where nudged
for bp, anchor, (xp, yp) in zip(branchpoints_sorted, anchors, best_pos):
    xa, ya = anchor
    pct_mods  = pct_mods_lookup.get(bp['label'], np.nan)
    pct_genes = pct_genes_lookup.get(bp['label'], np.nan)
    label_txt = f"{bp['label']}\n{pct_mods:.1f}% (m)\n{pct_genes:.1f}% (g)"

    dist = np.hypot(xp - xa, yp - ya)
    if dist > label_h * 0.05:
        ax.plot([xa, xp], [ya, yp], color=LINE_COLOR,
                linewidth=0.7, alpha=0.6, zorder=2)

    ax.text(
        xp, yp, label_txt,
        fontsize=font_size,
        ha='center', va='bottom',
        color=TEXT_COLOR,
        fontfamily=FONT_FAMILY,
        fontweight='bold',
        bbox=dict(
            boxstyle='round,pad=0.3',
            facecolor='none',
            alpha=1.0,
            edgecolor='none',
        ),
        zorder=4,
    )

plt.tight_layout()
plt.savefig(OUTPUT_FILE, bbox_inches='tight')
print(f"Saved {OUTPUT_FILE}")

svg_file = OUTPUT_FILE.replace('.pdf', '.svg')
plt.savefig(svg_file, bbox_inches='tight')
print(f"Saved {svg_file}")
