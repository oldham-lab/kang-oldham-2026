#!/usr/bin/env python3
"""
Figure assembly for Fig. 4 (v6), modeled on fig_3/v4/assemble_figure.py.

Source layout = the hand-laid-out fig_4/v5/v5.pptx (positions, sizes, panel
letters, "r = ..." labels, text). This script copies that pptx verbatim, swaps
ONLY the embedded panel media for the current fig_4/v6 panel SVGs (each panel is
stored as an SVG plus a PNG raster fallback), aspect-fits each swapped picture's
box so a regenerated panel isn't stretched, and writes the result to a NEW v6
deliverable (the v5 template is never modified).

Fig-4-specific notes vs. the fig_3 recipe:
  * There is no top schematic banner to remove (fig_3 deleted its panel A); all
    19 panels in the v5 layout are real outputs and are kept.
  * The "r = ..." correlation labels are kept -- the underlying correlation values
    did not change in v6, so the hand-placed labels remain valid.
  * The slide title is rewritten to the Fig. 4 caption and spanned across the full
    slide width at the top.

The 19 panels are 4 modules x 4 sub-panels (expression/bulk-cor/GSEA/projection)
plus the three summary panels Q (size distribution), R (seed-gene mean corr) and
S (projection correlation heatmap). The template media <-> source-panel mapping
below was recovered by content-hashing the embedded SVGs against fig_4/v5 and, for
the five panels that the v5 pptx embedded from an older export (the four projection
panels and S), by their column position on the slide. Panel S is stored as a shape
picture-fill (not a <p:pic>); its media is still swapped, and its 4:3 aspect is
unchanged so it needs no refit.
"""
import os
import zipfile, os, re, shutil, subprocess
import xml.etree.ElementTree as ET

TEMPLATE = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/v5/v5.pptx")  # source layout (read-only)
OUT      = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/v6/Kang_Figure_4_v6.pptx")
V6       = os.path.dirname(OUT)        # fig_4/v6 = canonical panel sources
STAGE    = "/tmp/fig4_assets"          # staging dir for placed assets (gitignored)

TITLE_NEW = ("Fig. 4 | Bulk gene coexpression modules in human neocortex reflect "
             "highly reproducible patterns of genomic activity among neocortical "
             "cell types")

# template media (svg vector, png fallback) -> source SVG in V6.
PANELS = [
    ("image14.svg", "image13.png", "panel_1_1.svg"),  # module 1, expression
    ("image12.svg", "image11.png", "panel_1_2.svg"),  # module 1, bulk cor
    ("image30.svg", "image29.png", "panel_1_3.svg"),  # module 1, GSEA
    ("image6.svg",  "image5.png",  "panel_1_4_with_brackets.svg"),  # module 1, projection + corr brackets
    ("image28.svg", "image27.png", "panel_2_1.svg"),  # module 2, expression
    ("image26.svg", "image25.png", "panel_2_2.svg"),  # module 2, bulk cor
    ("image32.svg", "image31.png", "panel_2_3.svg"),  # module 2, GSEA
    ("image4.svg",  "image3.png",  "panel_2_4_with_brackets.svg"),  # module 2, projection + corr brackets
    ("image24.svg", "image23.png", "panel_3_1.svg"),  # module 3, expression
    ("image22.svg", "image21.png", "panel_3_2.svg"),  # module 3, bulk cor
    ("image34.svg", "image33.png", "panel_3_3.svg"),  # module 3, GSEA
    ("image38.svg", "image37.png", "panel_3_4_with_brackets.svg"),  # module 3, projection + corr brackets
    ("image20.svg", "image19.png", "panel_4_1.svg"),  # module 4, expression
    ("image18.svg", "image17.png", "panel_4_2.svg"),  # module 4, bulk cor
    ("image36.svg", "image35.png", "panel_4_3.svg"),  # module 4, GSEA
    ("image2.svg",  "image1.png",  "panel_4_4_with_brackets.svg"),  # module 4, projection + corr brackets
    ("image10.svg", "image9.png",  "panel_Q.svg"),    # Q: module size distribution
    ("image8.svg",  "image7.png",  "panel_R.svg"),    # R: seed-gene mean corr
    ("image16.svg", "image15.png", "panel_S.svg"),    # S: projection corr heatmap
]

def _sh(cmd): subprocess.run(cmd, shell=True, check=True, capture_output=True)

# Stage each panel SVG from the v6 outputs and rasterize a matching PNG fallback.
def stage_assets():
    os.makedirs(STAGE, exist_ok=True)
    for _, _, src in PANELS:
        key = os.path.splitext(src)[0]
        shutil.copy(f"{V6}/{src}", f"{STAGE}/{key}.svg")
        # -b white: flatten on an opaque white background. A transparent fallback
        # (e.g. a panel SVG with no background rect) becomes a PDF soft-mask that
        # some viewers render as a black box; an opaque raster avoids that entirely.
        _sh(f"rsvg-convert -b white -d 300 -p 300 {STAGE}/{key}.svg -o {STAGE}/{key}.png")

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

# (1) Rewrite the slide title and span it across the full top of the slide:
# pin the box to the left margin, stretch it to the right margin, and raise its
# height so the long caption can wrap to a few lines.
title_done = False
for sp in slide.iter("{%s}sp" % NS["p"]):
    ts = sp.findall(".//a:t", NS)
    if ts and "".join(t.text or "" for t in ts).strip().startswith("Fig. 4"):
        ts[0].text = TITLE_NEW
        for t in ts[1:]: t.text = ""
        xfrm = sp.find(".//a:xfrm", NS)
        off, ext = xfrm.find("a:off", NS), xfrm.find("a:ext", NS)
        off.set("x", str(MARGIN))
        ext.set("cx", str(SLIDE_W - 2 * MARGIN))
        ext.set("cy", str(640080))   # ~0.7 in (room for ~2-3 wrapped lines)
        # The v5 template title box is wrap="none", so LibreOffice lays the long
        # caption on a single overflowing line (clipped) when exporting to PDF.
        # Force word-wrap so it breaks within the box width (PowerPoint already did).
        bodyPr = sp.find(".//a:bodyPr", NS)
        if bodyPr is not None:
            bodyPr.set("wrap", "square")
        title_done = True
        break

# (2) Remove the "r = ..." pairwise-correlation labels sitting to the right of the
# projection panels (m-p). The "Overall mean: r = ..." annotation by panel S does
# not start with "r =", so it is preserved.
parent = {c: p for p in slide.iter() for c in list(p)}   # child -> parent, for removal
removed_txt = 0
for sp in list(slide.iter("{%s}sp" % NS["p"])):
    txt = "".join(t.text or "" for t in sp.findall(".//a:t", NS)).strip()
    if re.match(r"^r\s*=", txt):
        parent[sp].remove(sp); removed_txt += 1

# (3) Aspect-fit each swapped <p:pic> box to its new asset, centered in its box.
# (Panel S is a shape blipFill, not a pic, and keeps its 4:3 aspect, so it is not
# refit here -- only its media bytes are swapped on write.)
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
    if panel.endswith("_4_with_brackets"):
        # m-p projection panels carry correlation brackets in a right margin, making
        # them wider than the placeholder box. Keep full box HEIGHT and grow width so
        # the bars stay full size and the brackets extend out (centered horizontally),
        # rather than letterboxing the whole thing shorter.
        ncy, ncx = cy, round(cy * asp)
    elif asp > box:          # too wide for box -> limit by width
        ncx, ncy = cx, round(cx/asp)
    else:                    # too tall for box -> limit by height
        ncy, ncx = cy, round(cy*asp)
    nox, noy = ox + (cx - ncx)//2, oy + (cy - ncy)//2
    off.set("x", str(nox)); off.set("y", str(noy))
    ext.set("cx", str(ncx)); ext.set("cy", str(ncy))
    fitted.append(f"{panel}: {cx}x{cy} -> {ncx}x{ncy} (centered)")

# (4) Lay out panel S. The panel_S.svg now carries the heatmap on the left and a
# right-hand column with the "Overall mean" label above the legend, so: (a) nudge the
# shape (image16 shape-fill) DOWN to clear the "Mean pairwise correlations of REI
# indices" title above it; (b) widen its box to the SVG's new 5.2:3 aspect, reclaiming
# the space the standalone "Overall mean" text used to occupy; and (c) delete that now-
# redundant template text box. (1 in = 914400 EMU.)
PANEL_S_DY     = round(0.06 * 914400)   # 0.06 in down (sits just under the title)
PANEL_S_ASPECT = 5.2 / 3                # matches ggsave(width = 5.2, height = 3)
parent = {c: p for p in slide.iter() for c in list(p)}   # refresh child -> parent
nudged = []
for sp in list(slide.iter("{%s}sp" % NS["p"])):
    b = sp.find(".//asvg:svgBlip", NS)
    if b is not None and rels.get(b.get("{%s}embed" % NS["r"])) == "image16.svg":
        xfrm = sp.find(".//a:xfrm", NS)
        off, ext = xfrm.find("a:off", NS), xfrm.find("a:ext", NS)
        off.set("y", str(int(off.get("y")) + PANEL_S_DY))
        ext.set("cx", str(round(int(ext.get("cy")) * PANEL_S_ASPECT)))
        nudged.append("panel_S +down, widened to %.3f aspect" % PANEL_S_ASPECT)
    txt = "".join(t.text or "" for t in sp.findall(".//a:t", NS)).strip()
    if txt.startswith("Overall mean"):
        parent[sp].remove(sp); nudged.append("overall-mean text removed (now in SVG)")

# (5) Panel-letter labels (a..p, q/r/s) must paint on top of the panel images.
# The GSEA panels (i-l) are now saved with an opaque white background, so wherever
# a picture box overlaps its letter the picture would cover it. Move each single-
# character label to the end of the shape tree (= topmost z-order).
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
# Write to a temp file then atomically move into place.
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
print("r= labels removed (right of m-p):", removed_txt)
print("panel-letter labels raised to top of z-order:", raised)
print("panel S declutter nudges:", nudged)
print("aspect-fit per panel:")
for f in fitted: print("  ", f)

# (6) Reproducible PDF + PNG deliverables. LibreOffice's headless SVG importer
# mis-renders the embedded panel SVGs (black boxes / spurious frames), so convert
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
