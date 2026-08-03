#!/usr/bin/env python3
"""
Reusable figure assembly for Fig. S3 (MTG + V1): take the representative
Kang_Figure_S3_v3.pptx layout and swap in the current panel outputs, writing back to the same
file. No manual steps. Modeled on fig_2/v4/assemble_figure.py.

Strategy: copy the pptx verbatim, replacing ONLY the embedded panel media (each
panel is stored as an SVG + a PNG raster fallback). Picture positions/sizes and
all text are inherited unchanged from the existing layout.

NB (differs from fig_2): no aspect-fit. The panels are emitted by fig_s3/v3.R at
fixed dimensions, so their aspect never drifts, and the authored picture boxes
are taken as-is. (fig_s4's tall donor-column boxes depend on this; fig_s3's
boxes already match the panel aspect, so it's moot here either way.)

Layout (slide 1), left -> right = subclass r2, subclass rmse, supertype r2,
supertype rmse:
  MTG main row   /  MTG sanity row  /  V1 main row  /  V1 sanity row
"""
import os
import zipfile, os, shutil, subprocess

BASE  = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s3/v3")
PPTX  = f"{BASE}/Kang_Figure_S3_v3.pptx"   # template and output (rewritten in place)
STAGE = "/tmp/figs3_assets"        # staging dir for placed assets (gitignored)

# template media (svg vector, png fallback) -> source panel SVG, relative to BASE.
# Order/mapping derived from the representative pptx (matched by content).
PANELS = [
    ("image16.svg", "image15.png", "MTG/part1_sub.svg"),            # MTG subclass r2
    ("image14.svg", "image13.png", "MTG/part1_sub_rmse.svg"),       # MTG subclass rmse
    ("image12.svg", "image11.png", "MTG/part1_super.svg"),          # MTG supertype r2
    ("image10.svg", "image9.png",  "MTG/part1_super_rmse.svg"),     # MTG supertype rmse
    ("image8.svg",  "image7.png",  "MTG/part1_sub_san.svg"),        # MTG subclass r2 (sanity)
    ("image6.svg",  "image5.png",  "MTG/part1_sub_san_rmse.svg"),   # MTG subclass rmse (sanity)
    ("image4.svg",  "image3.png",  "MTG/part1_super_san.svg"),      # MTG supertype r2 (sanity)
    ("image2.svg",  "image1.png",  "MTG/part1_super_san_rmse.svg"), # MTG supertype rmse (sanity)
    ("image32.svg", "image31.png", "V1/part1_sub.svg"),             # V1 subclass r2
    ("image30.svg", "image29.png", "V1/part1_sub_rmse.svg"),        # V1 subclass rmse
    ("image28.svg", "image27.png", "V1/part1_super.svg"),           # V1 supertype r2
    ("image26.svg", "image25.png", "V1/part1_super_rmse.svg"),      # V1 supertype rmse
    ("image24.svg", "image23.png", "V1/part1_sub_san.svg"),         # V1 subclass r2 (sanity)
    ("image22.svg", "image21.png", "V1/part1_sub_san_rmse.svg"),    # V1 subclass rmse (sanity)
    ("image20.svg", "image19.png", "V1/part1_super_san.svg"),       # V1 supertype r2 (sanity)
    ("image18.svg", "image17.png", "V1/part1_super_san_rmse.svg"),  # V1 supertype rmse (sanity)
]

def _sh(cmd): subprocess.run(cmd, shell=True, check=True, capture_output=True)

# Stage each panel SVG (flattened name, since filenames repeat across MTG/V1) and
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
_sh(f"soffice --headless -env:UserInstallation=file:///tmp/figs3_soffice "
    f"--convert-to pdf --outdir {BASE} {PPTX}")
_sh(f"pdftoppm -png -r 300 -singlefile {PDF} {PNG[:-4]}")
print(f"wrote: {PDF}\nwrote: {PNG}")
