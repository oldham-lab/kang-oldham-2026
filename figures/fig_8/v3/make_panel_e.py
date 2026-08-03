#!/usr/bin/env python3
"""
Panel e (Fig. 8) — GSEA summary table, cleaned for v3.

Source = the rendered gt table from fig_8/v1 (panel_G_gsea_summary_table_dfc_fdr.html).
fig8.R bakes two cosmetics into the SCZ-associated-genes column that v3 drops:
  * a red highlight on one gene  (<span style="color:red">RB1CC1</span>)
  * reference superscripts        (<sup>N</sup>)
Re-running fig8.R to change these would be heavy and would overwrite v1, so we
post-process the existing HTML (identical table content, just the two cosmetics
removed) and re-render to PDF with headless Chrome — the same engine the pipeline
(gt -> webshot2) used. assemble_figure.py crops this PDF to the table for panel e.
"""
import os
import re, subprocess

V1   = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v1")
V3   = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v3")
SRC_HTML = f"{V1}/panel_G_gsea_summary_table_dfc_fdr.html"
OUT_HTML = f"{V3}/panel_e_gsea_table.html"
OUT_PDF  = f"{V3}/panel_e_gsea_table.pdf"

html = open(SRC_HTML, encoding="utf-8").read()
html, n_red = re.subn(r'<span style="color:red">(.*?)</span>', r'\1', html)  # drop red highlight
html, n_sup = re.subn(r'<sup>.*?</sup>', '', html)                           # drop ref superscripts
open(OUT_HTML, "w", encoding="utf-8").write(html)
print(f"cleaned: removed {n_red} red span(s), {n_sup} reference superscript(s)")

# Render with webshot2 (chromote) — the same engine gt::gtsave used — cropping to
# the <table> element so it keeps its natural landscape proportions. (Chrome's
# --print-to-pdf instead constrains to the page width, wrapping the wide GO-geneset
# column into a tall/portrait block that then under-fills panel e's box.)
r_render = f'''
library(webshot2)
webshot("file://{OUT_HTML}", "{OUT_PDF}", selector = "table", zoom = 2)
'''
subprocess.run(["Rscript", "-e", r_render], check=True,
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
print("wrote:", OUT_PDF)
