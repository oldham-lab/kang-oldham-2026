#!/usr/bin/env python3
"""Fig 5 panel A, step 2/2: draw the REI-method schematic from panelA_cells.csv.

Reads the per-cell colours produced by build_panelA_cells.R and emits a
self-contained SVG: 8 heatmap grids laid out in the original three-column
schematic (projection -> pairwise correlation -> consensus) for Jorstad 2023
(red) and Gabitto 2024 (blue), plus the crossing/min/output arrows and labels.

Geometry (mm) was measured from the hand-built panel_A.svg so the output lines
up with the original; colorbar legends are intentionally omitted.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

# ---- document ------------------------------------------------------------
W, H = 169.60965, 87.396172

# ---- panel rectangles (mm): x, y = top-left; w, h = extent ---------------
# back/front pairs (cor matrices) drawn back-first so the front overlaps it.
PANELS = {
    "A1": dict(x=24.7,  y=17.1, w=25.0, h=25.0),   # Jorstad projection (red 5x4)
    "A2": dict(x=24.5,  y=62.3, w=25.0, h=25.0),   # Gabitto projection (blue 5x4)
    "A3": dict(x=69.8,  y=19.9, w=18.3, h=18.3),   # cor, center-top, back  (red)
    "A4": dict(x=72.3,  y=22.7, w=18.3, h=18.3),   # cor, center-top, front (blue)
    "A5": dict(x=69.9,  y=64.0, w=19.3, h=19.3),   # cor, center-bot, back  (red)
    "A6": dict(x=72.0,  y=66.3, w=19.3, h=19.3),   # cor, center-bot, front (blue)
    "A7": dict(x=112.2, y=16.7, w=25.0, h=25.0),   # consensus top    (grey)
    "A8": dict(x=112.2, y=62.2, w=25.0, h=25.0),   # consensus bottom (grey)
}
DRAW_ORDER = ["A1", "A2", "A3", "A5", "A4", "A6", "A7", "A8"]  # backs before fronts

# ---- type & colour -------------------------------------------------------
FONT = "Arial, Helvetica, sans-serif"
RED, BLUE = "#461210", "#333167"
SZ_TITLE, SZ_ROW, SZ_AXIS, SZ_MIN, SZ_BCD = 4.2, 4.4, 3.5, 4.0, 4.6


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def load_cells(path):
    cells = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            cells[r["panel"]].append(r)
    return cells


def grid(panel, rows):
    """Return SVG for one heatmap grid (cells + light gridlines + border)."""
    p = PANELS[panel]
    nr, nc = int(rows[0]["nrow"]), int(rows[0]["ncol"])
    cw, ch = p["w"] / nc, p["h"] / nr
    out = [f'<g>']
    for r in rows:
        i, j = int(r["row"]), int(r["col"])
        x = p["x"] + (j - 1) * cw
        y = p["y"] + (i - 1) * ch
        out.append(
            f'<rect x="{x:.3f}" y="{y:.3f}" width="{cw:.3f}" height="{ch:.3f}" '
            f'fill="{r["hex"]}" stroke="#ffffff" stroke-width="0.12"/>'
        )
    out.append(
        f'<rect x="{p["x"]:.3f}" y="{p["y"]:.3f}" width="{p["w"]:.3f}" '
        f'height="{p["h"]:.3f}" fill="none" stroke="#000000" stroke-width="0.35"/>'
    )
    out.append("</g>")
    return "\n".join(out)


def text(x, y, s, size, *, anchor="middle", italic=False, color="#000000", rot=None):
    style = (f'font-family:{FONT};font-size:{size}px;fill:{color};'
             f'text-anchor:{anchor};'
             f'font-style:{"italic" if italic else "normal"}')
    tr = f' transform="rotate(-90 {x:.3f} {y:.3f})"' if rot else ""
    return f'<text x="{x:.3f}" y="{y:.3f}" style="{style}"{tr}>{esc(s)}</text>'


def cx(panel):
    p = PANELS[panel]; return p["x"] + p["w"] / 2


def axis_labels(panel, top, left, *, left_panel=None):
    """Top label centred above panel; left label rotated to the left.

    left_panel lets a stacked pair anchor its rotated label to the back
    (left-most) matrix and centre the top label over the whole stack.
    """
    p = PANELS[panel]
    lp = PANELS[left_panel] if left_panel else p
    top_y = min(p["y"], lp["y"]) - 1.4
    top_x = (lp["x"] + p["x"] + p["w"]) / 2 if left_panel else cx(panel)
    out = [text(top_x, top_y, top, SZ_AXIS)]
    out.append(text(lp["x"] - 0.8, lp["y"] + lp["h"] / 2, left, SZ_AXIS, rot=True))
    return "\n".join(out)


def arrow(x1, y1, x2, y2, *, dotted=True, marker="dot"):
    dash = ' stroke-dasharray="1,1"' if dotted else ""
    w = 0.5 if dotted else 1.4
    return (f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="#000000" stroke-width="{w}"{dash} marker-end="url(#{marker})"/>')


def build(cells):
    body = []

    # arrowhead markers
    body.append(
        '<defs>'
        '<marker id="dot" markerWidth="6" markerHeight="6" refX="4.5" refY="3" '
        'orient="auto"><path d="M0,0 L5,3 L0,6 Z" fill="#000000"/></marker>'
        '<marker id="solid" markerWidth="5" markerHeight="5" refX="3.5" refY="2.5" '
        'orient="auto"><path d="M0,0 L4.5,2.5 L0,5 Z" fill="#000000"/></marker>'
        '</defs>'
    )

    # column titles (two lines, centred over each column)
    for label, xc in (("Projection", 37), ("Pairwise cor", 80), ("Consensus", 124.5)):
        body.append(text(xc, 5.4, label, SZ_TITLE))
        body.append(text(xc, 10.3, "matrices", SZ_TITLE))

    # row labels (italic, coloured), two lines
    body.append(text(4.5, 27.0, "Jorstad", SZ_ROW, anchor="start", italic=True, color=RED))
    body.append(text(4.5, 32.0, "2023",    SZ_ROW, anchor="start", italic=True, color=RED))
    body.append(text(4.5, 71.0, "Gabitto", SZ_ROW, anchor="start", italic=True, color=BLUE))
    body.append(text(4.5, 76.0, "2024",    SZ_ROW, anchor="start", italic=True, color=BLUE))

    # crossing dotted arrows: each projection feeds both cor stacks
    a1 = PANELS["A1"]; a2 = PANELS["A2"]
    r1y, r2y = a1["y"] + a1["h"] / 2, a2["y"] + a2["h"] / 2
    xr = a1["x"] + a1["w"]
    body.append(arrow(xr, r1y, 68.5, 27.0))       # Jorstad -> cor top
    body.append(arrow(xr, r1y + 1, 68.5, 70.0))   # Jorstad -> cor bottom (cross)
    body.append(arrow(xr, r2y - 1, 68.5, 28.5))   # Gabitto -> cor top (cross)
    body.append(arrow(xr, r2y, 68.5, 72.0))       # Gabitto -> cor bottom

    # min arrows cor -> consensus, with "min" label
    body.append(arrow(92.5, 29.0, 111.0, 29.0))
    body.append(text(101.0, 27.5, "min", SZ_MIN, italic=True))
    body.append(arrow(93.0, 74.5, 111.0, 74.5))
    body.append(text(102.0, 73.0, "min", SZ_MIN, italic=True))

    # solid output arrows consensus -> downstream panels
    body.append(arrow(138.0, 29.2, 147.0, 29.2, dotted=False, marker="solid"))
    body.append(text(149.4, 31.2, "(b)", SZ_BCD, anchor="start"))
    body.append(arrow(138.0, 74.7, 147.0, 74.7, dotted=False, marker="solid"))
    body.append(text(148.5, 76.5, "(c-d rows)", SZ_BCD, anchor="start"))

    # heatmap grids (backs before fronts)
    for panel in DRAW_ORDER:
        body.append(grid(panel, cells[panel]))

    # axis labels (match the original schematic's labelling)
    body.append(axis_labels("A1", "Cell types", "Modules (REI)"))
    body.append(axis_labels("A2", "Cell types", "Modules (REI)"))
    body.append(axis_labels("A4", "Modules", "Modules", left_panel="A3"))      # center-top stack
    body.append(axis_labels("A6", "Cell types", "Cell types", left_panel="A5"))  # center-bottom
    body.append(axis_labels("A7", "Modules", "Modules"))
    body.append(axis_labels("A8", "Cell types", "Cell types"))

    return body


def main():
    here = Path(__file__).resolve().parent
    cells = load_cells(here / "panelA_cells.csv")
    out = here / (sys.argv[1] if len(sys.argv) > 1 else "panel_A_generated.svg")
    svg = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{H}mm" '
        f'viewBox="0 0 {W} {H}" version="1.1">',
        *build(cells),
        "</svg>",
    ]
    out.write_text("\n".join(svg))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
