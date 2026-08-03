#!/usr/bin/env python3
"""
Reusable figure assembly for Fig. S5 (DFC + MTG): take the representative
Kang_Figure_S5_v2.pptx layout and swap in the current panel outputs, writing
back to the same file. No manual steps. Modeled on fig_s3/v3/assemble_figure.py.

Strategy: copy the pptx verbatim, replacing ONLY the embedded panel media (each
panel is stored as an SVG + a PNG raster fallback). Picture positions/sizes and
all text are inherited unchanged from the existing layout.

NB (as in fig_s3/fig_s4): no aspect-fit. Panels are emitted by fig_s5/v2/v2.R at
fixed dimensions and the authored picture boxes already match the panel aspect,
so the boxes are taken as-is.

Only the DFC and MTG regions appear in this figure; the dfc_10x/ and v1_10x/
panel sets in this folder are not part of the reference layout.

Layout (slide 1), left -> right = subclass r2, subclass rmse, supertype r2,
supertype rmse:
  DFC main row  /  DFC sanity row  /  MTG main row  /  MTG sanity row
"""
import os
import zipfile, os, shutil, subprocess

BASE  = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s5/v2")
PPTX  = f"{BASE}/Kang_Figure_S5_v2.pptx"   # template and output (rewritten in place)
STAGE = "/tmp/figs5_assets"                # staging dir for placed assets (gitignored)

# template media (svg vector, png fallback) -> source panel SVG, relative to BASE.
# Order/mapping derived from the representative pptx (matched by content).
PANELS = [
    ("image32.svg", "image31.png", "DFC/part1_sub.svg"),            # DFC subclass r2
    ("image30.svg", "image29.png", "DFC/part1_sub_rmse.svg"),       # DFC subclass rmse
    ("image28.svg", "image27.png", "DFC/part1_super.svg"),          # DFC supertype r2
    ("image26.svg", "image25.png", "DFC/part1_super_rmse.svg"),     # DFC supertype rmse
    ("image22.svg", "image21.png", "DFC/part1_sub_san.svg"),        # DFC subclass r2 (sanity)
    ("image24.svg", "image23.png", "DFC/part1_sub_san_rmse.svg"),   # DFC subclass rmse (sanity)
    ("image20.svg", "image19.png", "DFC/part1_super_san.svg"),      # DFC supertype r2 (sanity)
    ("image18.svg", "image17.png", "DFC/part1_super_san_rmse.svg"), # DFC supertype rmse (sanity)
    ("image8.svg",  "image7.png",  "MTG/part1_sub.svg"),            # MTG subclass r2
    ("image6.svg",  "image5.png",  "MTG/part1_sub_rmse.svg"),       # MTG subclass rmse
    ("image4.svg",  "image3.png",  "MTG/part1_super.svg"),          # MTG supertype r2
    ("image2.svg",  "image1.png",  "MTG/part1_super_rmse.svg"),     # MTG supertype rmse
    ("image10.svg", "image9.png",  "MTG/part1_sub_san.svg"),        # MTG subclass r2 (sanity)
    ("image12.svg", "image11.png", "MTG/part1_sub_san_rmse.svg"),   # MTG subclass rmse (sanity)
    ("image14.svg", "image13.png", "MTG/part1_super_san.svg"),      # MTG supertype r2 (sanity)
    ("image16.svg", "image15.png", "MTG/part1_super_san_rmse.svg"), # MTG supertype rmse (sanity)
]

def _sh(cmd): subprocess.run(cmd, shell=True, check=True, capture_output=True)

# Stage each panel SVG (flattened name, since filenames repeat across DFC/MTG) and
# rasterize a matching PNG fallback.
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
_sh(f"soffice --headless -env:UserInstallation=file:///tmp/figs5_soffice "
    f"--convert-to pdf --outdir {BASE} {PPTX}")
_sh(f"pdftoppm -png -r 300 -singlefile {PDF} {PNG[:-4]}")
print(f"wrote: {PDF}\nwrote: {PNG}")
