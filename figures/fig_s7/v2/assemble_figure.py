#!/usr/bin/env python3
"""
Figure assembly for Fig. S7 (v2), modeled on fig_4/v6/assemble_figure.py.

Source layout = the hand-laid-out Kang_Figure_S7_v2_template.pptx (positions,
sizes, title, subtitle). This script copies that pptx verbatim, swaps ONLY the
embedded panel media for the current fig_s7/v2 panel SVG (stored as an SVG plus a
PNG raster fallback), aspect-fits the swapped picture's box so a regenerated panel
isn't stretched, and writes the result to the Kang_Figure_S7_v2.pptx deliverable
(the *_template.pptx is never modified).

Fig-S7-specific notes vs. the fig_4 recipe:
  * S7 is a single panel (the per-subclass genome-wide expression barplots), so
    there are no panel-letter labels, no "r = ..." correlation labels, and no
    panel-S declutter step -- only the one <p:pic> is swapped and refit.
  * The title and subtitle text are already correct in the template, so they are
    left untouched (unlike fig_4, which rewrote its caption on the version bump).
"""
import os
import zipfile, os, re, shutil, subprocess
import xml.etree.ElementTree as ET

V2       = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v2")
TEMPLATE = f"{V2}/Kang_Figure_S7_v2_template.pptx"  # source layout (read-only)
OUT      = f"{V2}/Kang_Figure_S7_v2.pptx"           # deliverable
STAGE    = "/tmp/figs7_assets"                       # staging dir for placed assets

# template media (svg vector, png fallback) -> source SVG in V2.
PANELS = [
    ("image2.svg", "image1.png", "panel_1.svg"),  # per-subclass expression barplots
]

def _sh(cmd): subprocess.run(cmd, shell=True, check=True, capture_output=True)

# Stage each panel SVG from the v2 outputs and rasterize a matching PNG fallback.
def stage_assets():
    os.makedirs(STAGE, exist_ok=True)
    for _, _, src in PANELS:
        key = os.path.splitext(src)[0]
        shutil.copy(f"{V2}/{src}", f"{STAGE}/{key}.svg")
        # -b white: flatten on an opaque white background. A transparent fallback
        # becomes a PDF soft-mask that some viewers render as a black box; an opaque
        # raster avoids that entirely.
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

# Aspect-fit the swapped <p:pic> box to its new asset, centered in its box.
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
# and the edited slide XML. OUT != TEMPLATE, so the template source is untouched.
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
print("aspect-fit per panel:")
for f in fitted: print("  ", f)

# Reproducible PDF + PNG deliverables. LibreOffice's headless SVG importer
# mis-renders embedded panel SVGs (black boxes / spurious frames), so convert
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
