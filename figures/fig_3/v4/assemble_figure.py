#!/usr/bin/env python3
"""
Figure assembly for Fig. 3 (v4), modeled on fig_2/v4/assemble_figure.py.

Source layout = the hand-laid-out fig_3/v3/v3.pptx (positions, sizes, panel
letters, text). This script copies that pptx verbatim, swaps ONLY the embedded
panel media for the current fig_3/v4 panel SVGs (each panel is stored as an SVG
plus a PNG raster fallback), aspect-fits each swapped picture's box so a
regenerated panel isn't stretched, and writes the result to a NEW v4 deliverable
(the v3 template is never modified).

Fig-3-specific edits on top of the fig_2 recipe:
  * Panel A (top schematic) picture is swapped for layout_final.svg.
  * The slide title is rewritten to the Fig. 3 caption and spans the slide width.
  * The "r = ..." correlation labels next to the bottom two panels (E, F) are
    removed -- panels E and F now carry their brackets + r labels baked into
    their SVGs.

Layout (slide 1), by panel letter:
  a (top banner) .. layout_final.svg        (was image12.svg/image11.png)
  b ............... panel_B_bc.svg
  c ............... panel_C.svg
  d ............... panel_D.svg
  e ............... panel_E_native_log_with_brackets.svg
  f ............... panel_F_REI_with_brackets.svg
Note: panel_B_seed.svg exists in v4/ but is not placed in this layout (panel B
is panel_B_bc only), matching the v3 reference.
"""
import os
import zipfile, os, re, shutil, subprocess
import xml.etree.ElementTree as ET

TEMPLATE = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_3/v3/v3.pptx")   # source layout (read-only)
OUT      = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_3/v4/Kang_Figure_3_v4.pptx")
V4       = os.path.dirname(OUT)        # fig_3/v4 = canonical panel sources
STAGE    = "/tmp/fig3_assets"          # staging dir for placed assets (gitignored)

TITLE_NEW = ("Fig. 3 | Covariation Projection Analysis (CoPA) reveals the "
             "cellular origins of bulk gene coexpression modules")

# template media (svg vector, png fallback) -> source SVG in V4.
PANELS = [
    ("image12.svg", "image11.png", "layout_final.svg"),                    # Panel A (schematic banner)
    ("image6.svg",  "image5.png",  "panel_B_bc.svg"),                      # Panel B
    ("image10.svg", "image9.png",  "panel_C.svg"),                         # Panel C
    ("image8.svg",  "image7.png",  "panel_D.svg"),                         # Panel D
    ("image4.svg",  "image3.png",  "panel_E_native_log_with_brackets.svg"),# Panel E (+ corr brackets)
    ("image2.svg",  "image1.png",  "panel_F_REI_with_brackets.svg"),       # Panel F (+ corr brackets)
]

def _sh(cmd): subprocess.run(cmd, shell=True, check=True, capture_output=True)

# Stage each panel SVG from the v4 outputs and rasterize a matching PNG fallback;
# also author a blank placeholder (white box, faint border) for panel A.
def stage_assets():
    os.makedirs(STAGE, exist_ok=True)
    for _, _, src in PANELS:
        key = os.path.splitext(src)[0]
        shutil.copy(f"{V4}/{src}", f"{STAGE}/{key}.svg")
        # Rasterize with Inkscape, not rsvg-convert: rsvg mis-renders the panel A
        # schematic (layout_final.svg) -- its colored correlation matrices and
        # legend swatches come out solid black -- whereas Inkscape renders it
        # faithfully (matching PowerPoint / manual LibreOffice export).
        _sh(f"inkscape {STAGE}/{key}.svg --export-type=png --export-filename={STAGE}/{key}.png "
            f"--export-dpi=300 --export-background=white --export-background-opacity=1")
        # Inkscape still writes an alpha channel; strip it so the embedded raster is
        # fully opaque. A transparent PNG becomes a PDF soft-mask that some viewers
        # render as a black box, so flatten onto white and drop alpha entirely.
        _sh(f"convert {STAGE}/{key}.png -background white -alpha remove -alpha off {STAGE}/{key}.png")

stage_assets()

# template media path -> staged file to embed in its place
MEDIA_MAP = {}
for svg_media, png_media, src in PANELS:
    key = os.path.splitext(src)[0]
    MEDIA_MAP[f"ppt/media/{svg_media}"] = f"{key}.svg"
    MEDIA_MAP[f"ppt/media/{png_media}"] = f"{key}.png"
# (Panel A media is left in the package but unreferenced -- its picture is deleted below.)

# svg media basename -> panel key; aspect (w/h) read from each staged SVG.
SVG_PANEL = {svg_media: os.path.splitext(src)[0] for svg_media, _, src in PANELS}
def _svg_aspect(path):
    s = open(path, encoding="utf-8", errors="ignore").read(2000)
    w = float(re.search(r'width=["\']([\d.]+)', s).group(1))
    h = float(re.search(r'height=["\']([\d.]+)', s).group(1))
    return w / h
ASPECT = {os.path.splitext(src)[0]: _svg_aspect(f"{STAGE}/{os.path.splitext(src)[0]}.svg")
          for _, _, src in PANELS}

NS = {"a":"http://schemas.openxmlformats.org/drawingml/2006/main",
      "p":"http://schemas.openxmlformats.org/presentationml/2006/main",
      "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "asvg":"http://schemas.microsoft.com/office/drawing/2016/SVG/main"}
for k,v in NS.items(): ET.register_namespace("" if k=="" else k, v)

z = zipfile.ZipFile(TEMPLATE)

# rId -> media basename (from slide rels)
rels = {}
for rel in ET.fromstring(z.read("ppt/slides/_rels/slide1.xml.rels")):
    rels[rel.get("Id")] = rel.get("Target").split("/")[-1]

slide = ET.fromstring(z.read("ppt/slides/slide1.xml"))
parent = {c: p for p in slide.iter() for c in list(p)}   # child -> parent, for element removal

# slide width (EMU), for spanning the title across the top
SLIDE_W = int(ET.fromstring(z.read("ppt/presentation.xml")).find("p:sldSz", NS).get("cx"))
MARGIN  = 91440  # 0.1 in

# (1) Rewrite the slide title and span it across the top (extend the box to the
# right margin; raise its height so the long caption can wrap to ~2 lines).
TITLE_PT = 1200   # title font size in centipoints (12 pt; template default was 16)
title_done = False
for sp in slide.iter("{%s}sp" % NS["p"]):
    ts = sp.findall(".//a:t", NS)
    if ts and "".join(t.text or "" for t in ts).strip().startswith("Fig. 3"):
        ts[0].text = TITLE_NEW
        for t in ts[1:]: t.text = ""
        # shrink the title font on every run/end-paragraph property
        for rpr in sp.findall(".//a:rPr", NS) + sp.findall(".//a:endParaRPr", NS):
            rpr.set("sz", str(TITLE_PT))
        # the template box is wrap="none" (single line, overflows) -> let it wrap
        bodyPr = sp.find(".//a:bodyPr", NS)
        if bodyPr is not None:
            bodyPr.set("wrap", "square")
        xfrm = sp.find(".//a:xfrm", NS)
        off, ext = xfrm.find("a:off", NS), xfrm.find("a:ext", NS)
        ext.set("cx", str(SLIDE_W - int(off.get("x")) - MARGIN))
        ext.set("cy", str(640080))   # ~0.7 in
        title_done = True
        break

# (2) Remove the "r = ..." correlation labels sitting next to the bottom panels (E, F).
removed_txt = 0
for sp in list(slide.iter("{%s}sp" % NS["p"])):
    txt = "".join(t.text or "" for t in sp.findall(".//a:t", NS)).strip()
    if re.match(r"^r\s*=", txt):
        parent[sp].remove(sp); removed_txt += 1

# (3) Aspect-fit each swapped picture box to its new asset, centered in its box.
fitted = []
for pic in slide.iter("{%s}pic" % NS["p"]):
    svgblip = pic.find(".//asvg:svgBlip", NS)
    if svgblip is None:
        continue
    rid = svgblip.get("{%s}embed" % NS["r"])
    panel = SVG_PANEL.get(rels.get(rid, ""))
    if panel is None:
        continue
    xfrm = pic.find(".//a:xfrm", NS)
    off, ext = xfrm.find("a:off", NS), xfrm.find("a:ext", NS)
    ox, oy = int(off.get("x")), int(off.get("y"))
    cx, cy = int(ext.get("cx")), int(ext.get("cy"))
    asp, box = ASPECT[panel], cx/cy
    if asp > box:            # too wide for box -> limit by width
        ncx, ncy = cx, round(cx/asp)
    else:                    # too tall for box -> limit by height
        ncy, ncx = cy, round(cy*asp)
    nox, noy = ox + (cx - ncx)//2, oy + (cy - ncy)//2
    off.set("x", str(nox)); off.set("y", str(noy))
    ext.set("cx", str(ncx)); ext.set("cy", str(ncy))
    fitted.append(f"{panel}: {cx}x{cy} -> {ncx}x{ncy} (centered)")

slide_xml = (b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\r\n'
             + ET.tostring(slide, encoding="unicode").encode("utf-8"))

# Write the new pptx: copy everything from the template, replace the panel media
# (incl. blanked A) and the edited slide XML. OUT != TEMPLATE, so the v3 source
# is untouched. Write to a temp file then atomically move into place.
os.makedirs(os.path.dirname(OUT), exist_ok=True)
tmp = OUT + ".tmp"
with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
    for item in z.infolist():
        name = item.filename
        if name in MEDIA_MAP:
            data = open(f"{STAGE}/{MEDIA_MAP[name]}", "rb").read()
        elif name == "ppt/slides/slide1.xml":
            data = slide_xml
        else:
            data = z.read(name)
        out.writestr(item, data)
z.close()
os.replace(tmp, OUT)

print("wrote:", OUT)
print("title set:", title_done, "(spanning slide width) ->", TITLE_NEW)
print("r= labels removed:", removed_txt)
print("aspect-fit per panel:")
for f in fitted: print("  ", f)

# Reproducible PDF + PNG deliverables. LibreOffice's headless SVG importer
# mis-renders the embedded panel SVGs (heatmaps/violins -> black boxes), so convert
# from a raster-only copy of the pptx: strip the <asvg:svgBlip> extension from each
# blip so the primary <a:blip> (the rsvg-rendered PNG fallback staged above) is used.
# The SVG-based OUT pptx is the PowerPoint deliverable and is left untouched.
_SVGEXT = re.compile(r'<a:extLst>\s*<a:ext uri="\{96DAC541[^}]*\}">\s*'
                     r'<asvg:svgBlip\b[^>]*/>\s*</a:ext>\s*</a:extLst>')
base   = os.path.splitext(OUT)[0]
raster = base + ".raster.pptx"
zin = zipfile.ZipFile(OUT)
with zipfile.ZipFile(raster, "w", zipfile.ZIP_DEFLATED) as zo:
    for it in zin.infolist():
        d = zin.read(it.filename)
        if re.match(r"ppt/slides/slide\d+\.xml$", it.filename):
            d = _SVGEXT.sub("", d.decode("utf-8")).encode("utf-8")
        zo.writestr(it, d)
zin.close()
prof = "file://" + os.path.join(STAGE, "lo_profile")
_sh(f"soffice -env:UserInstallation={prof} --headless --convert-to pdf "
    f"--outdir {os.path.dirname(OUT)} {raster}")
os.replace(base + ".raster.pdf", base + ".pdf")
_sh(f"pdftoppm -png -singlefile -r 300 {base}.pdf {base}")
os.remove(raster)
print("wrote:", base + ".pdf", "and", base + ".png")
