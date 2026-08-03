#!/usr/bin/env python
# Composite the vector labeled dendrogram onto the (dendrogram-less) Jorstad panel_C,
# using UNIFORM scale (no text distortion): align the leaf axis [data 0..10n] to the
# measured heatmap column band [x0,x1], and drop the leaf baseline flush onto the top
# of the barplot stack.
#
# Writes two files (same coordinate system, so they stay column-aligned):
#   panel_C_consensusMin_Jorstad.svg  - full: labeled dendrogram + barplots + heatmap
#   panel_B.svg                       - cropped at the heatmap top: dendrogram + barplots
#
# Usage: python composite_dendrogram_panelC.py <v5_dir>
import os
import sys
from svgutils.transform import fromfile, SVGFigure

V5 = sys.argv[1] if len(sys.argv) > 1 else \
    os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v5")

def read_kv(path):
    d = {}
    for line in open(path):
        k, v = line.strip().split("=")
        d[k] = float(v)
    return d

pc = read_kv(f"{V5}/jorstad_panel_coords.txt")    # W_px,H_px,x0,x1,ybar_top,ybody_top
dc = read_kv(f"{V5}/jorstad_dendro_coords.txt")   # svg_w,svg_h,sx0,sx1,sy_leaf
Wp, Hp = pc["W_px"], pc["H_px"]
x0, x1, ybar, ybody = pc["x0"], pc["x1"], pc["ybar_top"], pc["ybody_top"]

scale = (x1 - x0) / (dc["sx1"] - dc["sx0"])       # uniform: leaf-span -> column band
tx = x0 - scale * dc["sx0"]
leaf_h = scale * dc["sy_leaf"]
Doff = max(0.0, leaf_h - ybar)                    # push panel down so dendro top >= 0
ty = (Doff + ybar) - leaf_h                       # leaves land exactly on the barplot top

def build(canvas_h):
    panel = fromfile(f"{V5}/panel_C_consensusMin_Jorstad_noDend.svg")
    dend  = fromfile(f"{V5}/jorstad_dendro_labeled.svg")
    proot, droot = panel.getroot(), dend.getroot()
    proot.moveto(0, Doff)
    droot.root.set("transform", f"translate({tx},{ty}) scale({scale})")
    fig = SVGFigure()
    fig.set_size((f"{Wp}pt", f"{canvas_h}pt"))
    fig.root.set("viewBox", f"0 0 {Wp} {canvas_h}")   # viewBox => content below canvas_h is clipped
    fig.append([proot, droot])
    return fig

# Full panel C (dendrogram + barplots + heatmap)
build(Doff + Hp).save(f"{V5}/panel_C_consensusMin_Jorstad.svg")
# Panel B: crop at the heatmap-body top -> dendrogram + barplots only
build(Doff + ybody).save(f"{V5}/panel_B.svg")
print(f"Saved panel_C_consensusMin_Jorstad.svg (H={Doff+Hp:.0f}) and panel_B.svg (H={Doff+ybody:.0f})"
      f"  [scale={scale:.4f}, band {x0:.0f}-{x1:.0f} of {Wp:.0f}]")
