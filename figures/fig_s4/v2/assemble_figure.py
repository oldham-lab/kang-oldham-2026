#!/usr/bin/env python3
"""
Reusable figure assembly for Fig. S4 (SEAAD2024 DFC): take the representative
Kang_Figure_S4_v2.pptx layout and swap in the current panel outputs, writing back to the same
file. No manual steps. Modeled on fig_2/v4/assemble_figure.py.

Strategy: copy the pptx verbatim, replacing ONLY the embedded panel media (each
panel is stored as an SVG + a PNG raster fallback). Picture positions/sizes and
all text are inherited unchanged from the existing layout.

NB (differs from fig_2): no aspect-fit. Each main panel is a tall 24-donor
column, so the author deliberately sized the picture boxes tall (the panel does
not sit at its natural aspect). Aspect-fitting would shrink them and break the
layout, so the authored boxes are taken as-is. Panels are emitted by
fig_s4/v2.R at fixed dimensions, so the embedded aspect is stable across rebuilds.

Layout (slide 1), left -> right = subclass r2, subclass rmse, supertype r2,
supertype rmse:
  main row  /  sanity row
"""
import os
import zipfile, os, shutil, subprocess

BASE  = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s4/v2")
PPTX  = f"{BASE}/Kang_Figure_S4_v2.pptx"   # template and output (rewritten in place)
STAGE = "/tmp/figs4_assets"        # staging dir for placed assets (gitignored)

# template media (svg vector, png fallback) -> source panel SVG, relative to BASE.
# Order/mapping derived from the representative pptx (matched by content).
PANELS = [
    ("image10.svg", "image9.png",  "part1_sub.svg"),            # subclass r2
    ("image8.svg",  "image7.png",  "part1_sub_rmse.svg"),       # subclass rmse
    ("image6.svg",  "image5.png",  "part1_super.svg"),          # supertype r2
    ("image4.svg",  "image3.png",  "part1_super_rmse.svg"),     # supertype rmse
    ("image2.svg",  "image1.png",  "part1_sub_san.svg"),        # subclass r2 (sanity)
    ("image12.svg", "image11.png", "part1_sub_san_rmse.svg"),   # subclass rmse (sanity)
    ("image14.svg", "image13.png", "part1_super_san.svg"),      # supertype r2 (sanity)
    ("image16.svg", "image15.png", "part1_super_san_rmse.svg"), # supertype rmse (sanity)
]

def _sh(cmd): subprocess.run(cmd, shell=True, check=True, capture_output=True)

# Stage each panel SVG (flattened name) and rasterize a matching PNG fallback.
def stage_assets():
    os.makedirs(STAGE, exist_ok=True)
    for _, _, src in PANELS:
        key = src.replace("/", "__")
        shutil.copy(f"{BASE}/{src}", f"{STAGE}/{key}")
        _sh(f"rsvg-convert -d 300 -p 300 {STAGE}/{key} -o {STAGE}/{key[:-4]}.png")

stage_assets()

# template media path -> staged file to embed in its place
MEDIA_MAP = {}
for svg_media, png_media, src in PANELS:
    key = src.replace("/", "__")
    MEDIA_MAP[f"ppt/media/{svg_media}"] = f"{STAGE}/{key}"
    MEDIA_MAP[f"ppt/media/{png_media}"] = f"{STAGE}/{key[:-4]}.png"

# Rewrite the pptx in place: write a temp copy (replacing only the panel media),
# then atomically replace -- PPTX is both source and destination.
z = zipfile.ZipFile(PPTX)
tmp = PPTX + ".tmp"
with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
    for item in z.infolist():
        repl = MEDIA_MAP.get(item.filename)
        data = open(repl, "rb").read() if repl else z.read(item.filename)
        out.writestr(item, data)
z.close()
os.replace(tmp, PPTX)
print(f"wrote: {PPTX} ({len(PANELS)} panels swapped)")

# Also emit PDF + PNG of the final figure (LibreOffice -> pdf, then rasterise).
PDF = PPTX[:-5] + ".pdf"
PNG = PPTX[:-5] + ".png"
_sh(f"soffice --headless -env:UserInstallation=file:///tmp/figs4_soffice "
    f"--convert-to pdf --outdir {BASE} {PPTX}")
_sh(f"pdftoppm -png -r 300 -singlefile {PDF} {PNG[:-4]}")
print(f"wrote: {PDF}\nwrote: {PNG}")
