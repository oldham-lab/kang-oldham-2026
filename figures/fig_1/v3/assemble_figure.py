#!/usr/bin/env python3
"""
Reusable figure assembly: take the v2.pptx layout (positions, sizes, text labels)
and swap in the current panel outputs, producing a new .pptx. No manual steps.

Strategy: copy the template pptx verbatim, replace ONLY the embedded panel media
(each panel is stored as an SVG + a PNG raster fallback), and aspect-fit each
picture's shape box so the swapped panel isn't stretched. All text labels,
panel letters, and coordinates are inherited unchanged from the template.
"""
import os
import zipfile, os, re, shutil, subprocess
import xml.etree.ElementTree as ET

# PowerPoint's SVG importer ignores `fill-opacity` in the inline style and drops the fill, so
# svglite's semi-transparent shapes (geom_bar alpha=0.6 bars, UpSet empty-intersection dots at
# alpha 0.5) vanish in the .pptx. librsvg/LibreOffice honor it, so the raster fallbacks are fine.
# Fix: pre-blend each semi-transparent fill against the white panel background into an opaque hex
# and drop fill-opacity. Visually identical in correct renderers; renders in PowerPoint too.
_STYLE = re.compile(r"style='([^']*)'")
_FILLHEX = re.compile(r"fill:\s*#([0-9A-Fa-f]{6})")
_FILLOP  = re.compile(r"fill-opacity:\s*([0-9]*\.?[0-9]+)\s*;?")
def _flatten_fill_opacity(svg_text, bg=(255, 255, 255)):
    def repl(m):
        style = m.group(1)
        fo, fh = _FILLOP.search(style), _FILLHEX.search(style)
        if fo is None or fh is None:
            return m.group(0)
        a = float(fo.group(1))
        if a < 1.0:
            r, g, b = (int(fh.group(1)[i:i+2], 16) for i in (0, 2, 4))
            blended = "#%02X%02X%02X" % tuple(round(c*a + bgc*(1-a)) for c, bgc in ((r, bg[0]), (g, bg[1]), (b, bg[2])))
            style = _FILLHEX.sub("fill: " + blended, style, count=1)
        style = _FILLOP.sub("", style)                 # drop fill-opacity (blended in, or a>=1 no-op)
        return "style='" + re.sub(r"\s+", " ", style).strip() + "'"
    return _STYLE.sub(repl, svg_text)

TEMPLATE = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/v2/v2.pptx")
OUT      = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/v3/fig1_v3_assembled.pptx")
V3       = os.path.dirname(OUT)   # fig1_v3.R output folder = canonical panel sources
STAGE    = "/tmp/fig_assets"      # staging dir for the placed assets (gitignored)

# Stage every panel from the fig1_v3.R outputs (no make_panel_C.R).
def _sh(cmd): subprocess.run(cmd, shell=True, check=True, capture_output=True)
def stage_assets():
    os.makedirs(STAGE, exist_ok=True)
    shutil.copy(f"{V3}/panel_A.svg", f"{STAGE}/A.svg")              # Panel A  (fig1_v3.R)
    shutil.copy(f"{V3}/panel_B.svg", f"{STAGE}/B.svg")              # Panel B  (fig1_v3.R)
    shutil.copy(f"{V3}/panel_C.svg", f"{STAGE}/C.svg")             # Panel C UpSet (fig1_v3.R)
    shutil.copy(f"{V3}/panel_EF_combined.svg", f"{STAGE}/EF.svg")  # E+F (fig1_v3.R)
    # Panel D: crop the full-page table PDF to its ink bounding box, then -> SVG
    bb = subprocess.run(f"gs -q -dBATCH -dNOPAUSE -sDEVICE=bbox '{V3}/panel_D.pdf'",
                        shell=True, capture_output=True, text=True).stderr
    x0, y0, x1, y1 = [float(v) for v in re.search(
        r"HiResBoundingBox: ([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+)", bb).groups()]
    pad = 2.0; w = x1 - x0 + 2*pad; h = y1 - y0 + 2*pad
    _sh(f"gs -q -o {STAGE}/D_crop.pdf -sDEVICE=pdfwrite -dFIXEDMEDIA "
        f"-dDEVICEWIDTHPOINTS={w:.2f} -dDEVICEHEIGHTPOINTS={h:.2f} "
        f"-c '<</PageOffset [{-(x0-pad):.2f} {-(y0-pad):.2f}]>> setpagedevice' -f '{V3}/panel_D.pdf'")
    _sh(f"pdf2svg {STAGE}/D_crop.pdf {STAGE}/D.svg")
    for k in ("A", "B", "C", "D", "EF"):                            # flatten fill-opacity in place (PowerPoint fix)
        p = f"{STAGE}/{k}.svg"
        with open(p, encoding="utf-8") as fh: flat = _flatten_fill_opacity(fh.read())
        with open(p, "w", encoding="utf-8") as fh: fh.write(flat)
    for k in ("A", "B", "C", "D", "EF"):                            # PNG raster fallbacks
        _sh(f"rsvg-convert -d 300 -p 300 {STAGE}/{k}.svg -o {STAGE}/{k}.png")

stage_assets()
ASSETS = STAGE

# Which template media file gets which staged asset.
# (template pairs each picture as imageN.svg [vector] + imageM.png [fallback])
MEDIA_MAP = {
    "ppt/media/image4.svg": "A.svg",  "ppt/media/image3.png": "A.png",   # Panel A
    "ppt/media/image2.svg": "B.svg",  "ppt/media/image1.png": "B.png",   # Panel B
    "ppt/media/image6.svg": "C.svg",  "ppt/media/image5.png": "C.png",   # Panel C
    "ppt/media/image8.svg": "D.svg",  "ppt/media/image7.png": "D.png",   # Panel D
    "ppt/media/image10.svg": "EF.svg","ppt/media/image9.png": "EF.png",  # Panels E+F
}
# svg media basename -> panel key; aspect (w/h) read from each staged SVG
SVG_PANEL = {"image4.svg":"A","image2.svg":"B","image6.svg":"C","image8.svg":"D","image10.svg":"EF"}
def _svg_aspect(path):
    s = open(path, encoding="utf-8", errors="ignore").read(2000)
    w = float(re.search(r'width=["\']([\d.]+)', s).group(1))
    h = float(re.search(r'height=["\']([\d.]+)', s).group(1))
    return w / h
ASPECT = {k: _svg_aspect(f"{STAGE}/{k}.svg") for k in ("A","B","C","D","EF")}

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

# Edit slide XML: aspect-fit each picture box to its new asset (preserve top-left)
slide = ET.fromstring(z.read("ppt/slides/slide1.xml"))
fitted = []
# Nudge Panel C down so it clears its "Overlap of unique subclass markers" title.
PANEL_C_DOWN_EMU = 137160   # ~0.15 in (tunable)
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
    # Placement within the original box: most panels are centered so a shrunk panel
    # doesn't hug the top/left of its slot. Panel C is anchored at the original
    # top-left so it sits in the exact same spot as the original PowerPoint.
    if panel == "C":
        nox, noy = ox, oy + PANEL_C_DOWN_EMU   # nudge down so C clears its title
    else:
        nox, noy = ox + (cx - ncx)//2, oy + (cy - ncy)//2
    off.set("x", str(nox)); off.set("y", str(noy))
    ext.set("cx", str(ncx)); ext.set("cy", str(ncy))
    fitted.append(f"{panel}: {cx}x{cy} -> {ncx}x{ncy} ({'anchored@orig' if panel=='C' else 'centered'})")
    if panel == "C":
        cbox = (nox, noy, ncx, ncy)   # placed box of Panel C, for the y-axis title
    if panel == "EF":
        efbox = (nox, noy, ncx, ncy)  # placed box of E+F, for label alignment

# Panel C y-axis title: UpSetR's title is suppressed (its positioning is unfixable),
# so we add it here as a rotated text box in the empty space left of the bars.
EMU = 914400
CY_TITLE = {"text":"# of intersecting genes", "sz":500,
            "fx":0.315, "fy":0.26,     # center as fraction of C's placed box (tune); +fx = toward axis
            "w_in":0.70, "h_in":0.16}  # unrotated box (long side = text length)
cx0, cy0, cw, ch = cbox
cenx = cx0 + CY_TITLE["fx"]*cw
ceny = cy0 + CY_TITLE["fy"]*ch
bw, bh = int(CY_TITLE["w_in"]*EMU), int(CY_TITLE["h_in"]*EMU)
tox, toy = int(cenx - bw/2), int(ceny - bh/2)
sp_xml = (
 '<p:sp xmlns:p="{p}" xmlns:a="{a}">'
 '<p:nvSpPr><p:cNvPr id="90" name="C_yaxis_title"/><p:cNvSpPr txBox="1"/><p:nvPr/></p:nvSpPr>'
 '<p:spPr><a:xfrm rot="16200000"><a:off x="{x}" y="{y}"/><a:ext cx="{w}" cy="{h}"/></a:xfrm>'
 '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr>'
 '<p:txBody><a:bodyPr wrap="none" lIns="0" tIns="0" rIns="0" bIns="0" anchor="ctr"/>'
 '<a:p><a:pPr algn="ctr"/><a:r><a:rPr lang="en-US" sz="{sz}"/><a:t>{txt}</a:t></a:r></a:p>'
 '</p:txBody></p:sp>'
).format(p=NS["p"], a=NS["a"], x=tox, y=toy, w=bw, h=bh, sz=CY_TITLE["sz"], txt=CY_TITLE["text"])
spTree = slide.find(".//p:cSld/p:spTree", NS)
spTree.append(ET.fromstring(sp_xml))

# Title: set the full figure title and make it span the top of the slide. The template box
# was sized (0.86 in, wrap="none" + spAutoFit) for the short "Fig. 1 v2" placeholder, which
# leaves the long title hugging the left edge. Widen the box to the full slide width (minus
# small side margins), center the text, and switch to normal wrapping so it spans the top.
TITLE = "Fig. 1 | Cell type marker genes in human neocortex vary by choice of study and algorithm"
SLIDE_W = int(ET.fromstring(z.read("ppt/presentation.xml")).find("p:sldSz", NS).get("cx"))
TITLE_MARGIN = 91440  # 0.1 in side margins
for sp in slide.iter("{%s}sp" % NS["p"]):
    runs = sp.findall(".//a:t", NS)
    joined = "".join(t.text or "" for t in runs)
    if joined.strip().startswith("Fig"):
        runs[0].text = TITLE
        for t in runs[1:]:
            t.text = ""
        xfrm = sp.find(".//a:xfrm", NS)
        off, ext = xfrm.find("a:off", NS), xfrm.find("a:ext", NS)
        off.set("x", str(TITLE_MARGIN))
        ext.set("cx", str(SLIDE_W - 2 * TITLE_MARGIN))   # span full width
        bodyPr = sp.find(".//a:bodyPr", NS)
        bodyPr.set("wrap", "square")                      # allow wrap within the wide box
        af = bodyPr.find("a:spAutoFit", NS)
        if af is not None:                                # drop autofit so width is honored
            bodyPr.remove(af)
        for p in sp.findall(".//a:p", NS):                # left-justify each paragraph
            pPr = p.find("a:pPr", NS)
            if pPr is None:
                pPr = ET.Element("{%s}pPr" % NS["a"]); p.insert(0, pPr)
            pPr.set("algn", "l")
        print(f"title: '{joined}' -> '{TITLE}' (full-width, left-justified)")
        break

# Programmatic E/F label alignment: the cell-class Set1 strip spans exactly the heatmap
# body width, so its x-extent (per panel) gives each body's center. Map those centers into
# the placed EF picture box, then give F's panel letter + title the SAME body-relative
# offset that E's have (E aligns correctly in the template) so F matches. Re-measured each
# build, so it adapts if the combined-SVG layout changes.
def ef_body_centers(svg_path):
    from collections import Counter
    s = open(svg_path, encoding="utf-8", errors="ignore").read()
    VW = float(re.search(r"viewBox='0 0 ([\d.]+) ", s).group(1))
    S1 = ("#4DAF4A", "#377EB8", "#E41A1C")   # Cell-class Set1 strip colors
    rects = [(float(m[1]), float(m[2]), float(m[3]))
             for m in re.finditer(r"<rect x='([\d.]+)' y='([\d.]+)' width='([\d.]+)' "
                                  r"height='[\d.]+' style='[^']*fill: (#[0-9A-Fa-f]{6})", s)
             if m[4].upper() in S1]
    strip_y = Counter(round(r[1], 1) for r in rects).most_common(1)[0][0]   # strip = densest row
    strip = sorted(r for r in rects if abs(r[1] - strip_y) < 0.5)
    gi = max(range(len(strip) - 1), key=lambda i: strip[i+1][0] - (strip[i][0] + strip[i][2]))
    cen = lambda g: (min(r[0] for r in g) + max(r[0] + r[2] for r in g)) / 2 / VW
    return cen(strip[:gi+1]), cen(strip[gi+1:])   # (E, F) center fractions of SVG width

def _find_sp(pred):
    for sp in slide.iter("{%s}sp" % NS["p"]):
        if pred("".join(t.text or "" for t in sp.findall(".//a:t", NS)).strip()):
            return sp
    return None

# Gene counts written by fig*_v*.R ({dataset: nrow of the metacell matrix}); used to rewrite
# "n = <count> genes" in the E/F titles so they always match the data / region.
def read_gene_counts(path):
    d = {}
    try:
        for line in open(path):
            if "=" in line: k, v = line.strip().split("=", 1); d[k] = int(v)
    except FileNotFoundError:
        pass
    return d
GENE_COUNTS = read_gene_counts(f"{V3}/panel_EF_gene_counts.txt")

try:
    E_cen_f, F_cen_f = ef_body_centers(f"{STAGE}/EF.svg")
    efx, _, efcx, _ = efbox
    E_cen, F_cen = efx + E_cen_f * efcx, efx + F_cen_f * efcx   # body centers in slide EMU
    # NB: Panels A/B titles also contain "Jorstad"/"Gabitto" — match the E/F heatmap titles
    # specifically by "Pairwise correlations" + the dataset name.
    shapes = {"e": _find_sp(lambda t: t == "e"), "f": _find_sp(lambda t: t == "f"),
              "eT": _find_sp(lambda t: "Pairwise" in t and "Jorstad" in t),
              "fT": _find_sp(lambda t: "Pairwise" in t and "Gabitto" in t)}
    # Rewrite "n = <count> genes" in each E/F title from the sidecar (region-correct gene counts).
    for key, ds in (("eT", "jorstad"), ("fT", "gabitto")):
        if shapes[key] is not None and ds in GENE_COUNTS:
            for t in shapes[key].findall(".//a:t", NS):
                if t.text and re.search(r"n\s*=\s*[\d,]+\s*genes", t.text):
                    t.text = re.sub(r"n\s*=\s*[\d,]+\s*genes", f"n = {GENE_COUNTS[ds]:,} genes", t.text)
                    print(f"gene count [{ds}] -> {GENE_COUNTS[ds]:,}")
    if all(shapes.values()):
        ox = lambda sp: sp.find(".//a:off", NS)
        letter_off = int(ox(shapes["e"]).get("x")) - E_cen   # E's letter offset from its body
        title_off  = int(ox(shapes["eT"]).get("x")) - E_cen  # E's title offset from its body
        f_x0, fT_x0 = int(ox(shapes["f"]).get("x")), int(ox(shapes["fT"]).get("x"))
        ox(shapes["f"]).set("x",  str(int(F_cen + letter_off)))
        ox(shapes["fT"]).set("x", str(int(F_cen + title_off)))
        print(f"aligned F (E/F body fracs {E_cen_f:.3f}/{F_cen_f:.3f}): "
              f"letter {f_x0}->{int(F_cen+letter_off)}, title {fT_x0}->{int(F_cen+title_off)} EMU")
    else:
        print("WARN: could not find all e/f letter+title shapes; skipped alignment",
              {k: v is not None for k, v in shapes.items()})
except Exception as e:
    print("WARN: E/F label alignment skipped:", e)

slide_xml = (b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\r\n'
             + ET.tostring(slide, encoding="unicode").encode("utf-8"))

# Write the new pptx: copy everything, replace the 10 media files + slide XML
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as out:
    for item in z.infolist():
        name = item.filename
        if name in MEDIA_MAP:
            data = open(f"{ASSETS}/{MEDIA_MAP[name]}", "rb").read()
        elif name == "ppt/slides/slide1.xml":
            data = slide_xml
        else:
            data = z.read(name)
        out.writestr(item, data)

print("wrote:", OUT)
print("aspect-fit per panel:")
for f in fitted: print("  ", f)

# Reproducible PDF + PNG deliverables (modeled on fig_3/v4/assemble_figure.py).
# LibreOffice's headless SVG importer mis-renders the embedded panel SVGs, so convert
# from a raster-only copy of the pptx: strip the <asvg:svgBlip> extension from each blip
# so the primary <a:blip> (the rsvg-rendered PNG fallback staged above) is used instead.
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
