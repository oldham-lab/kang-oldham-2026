#!/usr/bin/env python3
"""
Figure assembly for Fig. 6 (v6), modeled on fig_4/v6/assemble_figure.py.

Source layout = the hand-laid-out fig_6/v5/v5.pptx (panel positions, sizes, panel
letters, slide title). This script copies that pptx verbatim, swaps ONLY the
embedded panel media for the current fig_6/v6 panel SVGs (each panel is stored as
an SVG plus a PNG raster fallback), aspect-fits each swapped picture's box so a
regenerated panel isn't stretched, and writes the result to a NEW v6 deliverable
(the v5 template is never modified).

The v5 layout has three panels:
  * a  (top banner, full width) -- the dCoPA schematic
  * left column                 -- one module snapshot (expr / GSEA / projection /
                                   boxplots), sub-lettered b/d/f/h
  * right column                -- a second module snapshot, sub-lettered c/e/g/i

For v6 the schematic is decopa_schematic_final.svg, and the two columns are the
fig6_specific_modules.R outputs for modules 83 (left) and 147 (right). The module
snapshots are ggsaved at 6x14 in, matching the v5 column boxes (3.12x7.28 in), so
the aspect-fit is effectively a no-op but is kept for robustness.
"""
import os
import zipfile, os, re, shutil, subprocess
import xml.etree.ElementTree as ET

TEMPLATE = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/v5/v5.pptx")  # source layout (read-only)
OUT      = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/v6/Kang_Figure_6_v6.pptx")
V6       = os.path.dirname(OUT)        # fig_6/v6 = canonical panel sources
STAGE    = "/tmp/fig6_assets"          # staging dir for placed assets (gitignored)

TITLE_NEW = ("Fig. 6 | Differential CoPA (dCoPA) identifies gene coexpression modules that "
             "are significantly and uniformly dysregulated in specific cell types in disease")

# template media (svg vector, png fallback) -> source SVG in V6.
PANELS = [
    ("image2.svg", "image1.png", "decopa_schematic_final.svg"),        # a: dCoPA schematic
    ("image4.svg", "image3.png", "83_MTC_allAD_Con_bulk_megaset.svg"),  # left column: module 83
    ("image6.svg", "image5.png", "147_MTC_allAD_Con_bulk_megaset.svg"), # right column: module 147
]

def _sh(cmd): subprocess.run(cmd, shell=True, check=True, capture_output=True)

# Stage each panel SVG from the v6 outputs and rasterize a matching PNG fallback.
def stage_assets():
    os.makedirs(STAGE, exist_ok=True)
    for _, _, src in PANELS:
        key = os.path.splitext(src)[0]
        shutil.copy(f"{V6}/{src}", f"{STAGE}/{key}.svg")
        # Rasterize with Inkscape, not rsvg-convert: rsvg mis-renders the panel A
        # schematic (decopa_schematic_final.svg) -- its colored boxplots and legend
        # swatches come out solid black -- whereas Inkscape renders it faithfully
        # (matching PowerPoint / manual LibreOffice export). See fig_3/v4.
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
      "asvg":"http://schemas.microsoft.com/office/drawing/2016/SVG/main",
      "mc":"http://schemas.openxmlformats.org/markup-compatibility/2006",
      "p14":"http://schemas.microsoft.com/office/powerpoint/2010/main",
      "p15":"http://schemas.microsoft.com/office/powerpoint/2012/main"}
for k,v in NS.items(): ET.register_namespace(k, v)

z = zipfile.ZipFile(TEMPLATE)

# rId -> media basename (from slide rels)
rels = {}
for rel in ET.fromstring(z.read("ppt/slides/_rels/slide1.xml.rels")):
    rels[rel.get("Id")] = rel.get("Target").split("/")[-1]

slide = ET.fromstring(z.read("ppt/slides/slide1.xml"))

# slide width (EMU), for spanning the title across the top
SLIDE_W = int(ET.fromstring(z.read("ppt/presentation.xml")).find("p:sldSz", NS).get("cx"))
MARGIN  = 91440  # 0.1 in

# (1) Rewrite the slide title and span it across the full top of the slide: pin the
# box to the left margin, stretch it to the right margin, shrink the font so the long
# caption fits, and force word-wrap (the v5 box is wrap="none", a single overflowing
# line). The box is top-anchored, so the extra height never pushes text onto panel a.
TITLE_PT = 1200   # 12 pt in centipoints (v5 template default was 16 pt)
title_done = False
for sp in slide.iter("{%s}sp" % NS["p"]):
    ts = sp.findall(".//a:t", NS)
    if ts and "".join(t.text or "" for t in ts).strip().startswith("Fig. 6"):
        ts[0].text = TITLE_NEW
        for t in ts[1:]: t.text = ""
        for rpr in sp.findall(".//a:rPr", NS) + sp.findall(".//a:endParaRPr", NS):
            rpr.set("sz", str(TITLE_PT))
        bodyPr = sp.find(".//a:bodyPr", NS)
        if bodyPr is not None:
            bodyPr.set("wrap", "square")
        xfrm = sp.find(".//a:xfrm", NS)
        off, ext = xfrm.find("a:off", NS), xfrm.find("a:ext", NS)
        off.set("x", str(MARGIN))
        ext.set("cx", str(SLIDE_W - 2 * MARGIN))
        ext.set("cy", str(640080))   # ~0.7 in (room for the caption to wrap)
        title_done = True
        break

# (2) Aspect-fit each swapped <p:pic> box to its new asset, centered in its box.
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

# (2b) The 2-line spanning title needs more vertical room than the v5 1-line title,
# so nudge panel A (schematic picture + its "a" letter) down to clear it. Panel A's
# bottom stays well above the two module columns (y=3.08 in), so nothing else moves.
PANEL_A_DY = round(0.16 * 914400)   # 0.16 in down
for pic in slide.iter("{%s}pic" % NS["p"]):
    svgblip = pic.find(".//asvg:svgBlip", NS)
    if svgblip is None:
        continue
    if SVG_PANEL.get(rels.get(svgblip.get("{%s}embed" % NS["r"]), "")) == "decopa_schematic_final":
        off = pic.find(".//a:xfrm/a:off", NS)
        off.set("y", str(int(off.get("y")) + PANEL_A_DY))
for sp in slide.iter("{%s}sp" % NS["p"]):
    if "".join(t.text or "" for t in sp.findall(".//a:t", NS)).strip() == "a":
        off = sp.find(".//a:xfrm/a:off", NS)
        if off is not None:
            off.set("y", str(int(off.get("y")) + PANEL_A_DY))

# (3) Panel-letter labels (a..i) must paint on top of the panel images: the PNG
# fallbacks are opaque white, so a letter sitting behind a picture box would be
# covered. Move each single-character label to the end of the shape tree (topmost).
spTree = slide.find(".//p:spTree", NS)
raised = 0
for sp in list(spTree.findall("p:sp", NS)):
    txt = "".join(t.text or "" for t in sp.findall(".//a:t", NS)).strip()
    if len(txt) == 1 and txt.isalpha():
        spTree.remove(sp); spTree.append(sp); raised += 1

slide_xml = (b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\r\n'
             + ET.tostring(slide, encoding="unicode").encode("utf-8"))

# Write the new pptx: copy everything from the template, replace the panel media
# and the edited slide XML. OUT != TEMPLATE, so the v5 source is untouched.
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
print("title set:", title_done, "->", TITLE_NEW)
print("panel-letter labels raised to top of z-order:", raised)
print("aspect-fit per panel:")
for f in fitted: print("  ", f)

# (4) Reproducible PDF + PNG deliverables. LibreOffice's headless SVG importer
# mis-renders embedded panel SVGs (black boxes / spurious frames), so convert from
# a raster-only copy of the pptx: strip the <asvg:svgBlip> extension from each blip
# so the primary <a:blip> (the rsvg-rendered PNG fallback) is used. The SVG-based
# OUT pptx is the PowerPoint deliverable and is left untouched.
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
