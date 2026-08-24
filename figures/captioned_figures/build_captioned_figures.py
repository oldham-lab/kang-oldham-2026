#!/usr/bin/env python3
"""Build standalone captioned versions of the 19 finalized manuscript figures.

Each output is ONE letter page: the figure with its baked-in title band cropped
off at the top, then the figure number + title in bold and the full legend
beneath it -- the layout used for the figures embedded in the dissertation
(docs/thesis), but as a self-contained file per figure.

  figures  Code_for_figures/Kang_Oldham_Figures/Kang_Figure_<N>.{pdf,png}
           (the finalized staging set that also feeds the zip bundle, the Box
           upload and the thesis build)
  legends  the FIGURE LEGENDS and EXTENDED DATA FIGURE LEGENDS sections of
           LEGEND_PDF, copied verbatim -- see SOURCE_NOTES.md for why the
           legends come from that PDF rather than the paper_draft docx, and
           for what differs between the two

Captions keep JOURNAL numbering ("Fig. 6 | ..."), matching the manuscript and
the Kang_Figure_<N> filenames -- not the thesis's UCSF "Figure 2.6:" scheme.

The figure is placed from the source PDF, so it stays vector; only S3/S5 need a
painted-over region (see LEFTLABEL). Since everything must fit on one page,
tall figures are scaled down to leave room for their legend.

Usage:  python3 build_captioned_figures.py [--out DIR]
"""
import argparse
import os
import re
import subprocess

import fitz
import numpy as np
from PIL import Image, ImageDraw

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAGE = os.path.join(BASE, "Kang_Oldham_Figures")
LEGEND_PDF = "/home/gugene/docs/Kang&Oldham_08.18.26.pdf"
FONT_DIR = "/usr/share/fonts/truetype/liberation2"

FIGS = ['1', '2', '3', '4', '5', '6', '7', '8',
        'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11']

# ---- page geometry (Letter, 1 in margins), all in points ----
PAGE_W, PAGE_H = 612.0, 792.0
MARGIN = 72.0
TEXT_W = PAGE_W - 2 * MARGIN          # 468
TEXT_H = PAGE_H - 2 * MARGIN          # 648
IMG_AFTER = 8.0                       # gap between the figure and its caption
INSET = 1.0                           # insert_htmlbox insets its text by ~1 pt


# ---------------------------------------------------------------- legends ----
# The legends are lifted verbatim out of the manuscript PDF, styling included:
# bold, superscript reference numbers, the log2 subscripts, and the monospaced
# tool names (edgeR / DESeq2) the manuscript sets in Courier.
SEC_START = 'FIGURE LEGENDS'
SEC_MID = 'EXTENDED DATA FIGURE LEGENDS'          # second half of the same run
SEC_END = 'EXTENDED DATA TABLE LEGENDS'           # stop here

# Line numbers sit at x0 in [28, 35] and page numbers below y=730, both left of
# / outside the 72 pt text block. Filtering them by POSITION is essential: they
# are set in Aptos, but so are a few characters of real legend text (the ")."
# closing the CoPA Cabana URL in Fig. 7 / Fig. 8 / Fig. S10), so filtering by
# font silently truncates those three legends.
TEXT_LEFT = 60.0
TEXT_BOTTOM = 730.0

# The source bolds a legend's opening panel marker along with the title, so it
# has to be handed back to the legend (same rule as the thesis build).
PANEL_TAIL_RE = re.compile(r'(?<=\.)(\s+[a-z](?:\s*[-–]\s*[a-z])?\)?\s*)$')
FIG_RE = re.compile(r'^Fig\.\s+(S?\d+)\s*\|')

# The ONLY departure from copying the legends verbatim. Fig. S1's title in the
# source PDF still carries a literal "(ref)." placeholder -- it was deleted from
# the paper_draft docx on 2026-08-17 but not from this PDF. Dropping it here is
# a stopgap: fixing the source PDF would make this entry unnecessary, and the
# build reports whether it actually fired.
TITLE_FIXES = {'S1': (' (ref).', '')}


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def legend_runs(pdf):
    """Styled runs of both figure-legend sections, in reading order.

    Each run is (text, bold, mono, script) with script in (None, 'sup', 'sub').
    """
    doc = fitz.open(pdf)
    runs = []
    started = False
    for page in doc:
        for block in page.get_text('dict')['blocks']:
            if block['type'] != 0:
                continue
            for line in block['lines']:
                spans = [s for s in line['spans']
                         if s['bbox'][0] >= TEXT_LEFT and s['bbox'][3] <= TEXT_BOTTOM]
                if not spans:
                    continue
                heading = ''.join(s['text'] for s in spans).strip()
                if heading == SEC_START:
                    started = True
                    continue
                if heading == SEC_END:
                    doc.close()
                    return runs
                if heading == SEC_MID or not started:
                    continue
                # super- vs subscript is a baseline shift, not just a size drop:
                # the log2 subscripts are the same 7.92 pt as the citations.
                body = max(s['size'] for s in spans)
                bases = [s['origin'][1] for s in spans if abs(s['size'] - body) < 0.5]
                base = max(set(bases), key=bases.count)
                for s in spans:
                    script = None
                    if s['size'] < body - 0.5:
                        dy = s['origin'][1] - base
                        script = 'sup' if dy < -0.5 else ('sub' if dy > 0.5 else None)
                    runs.append((s['text'],
                                 bool(s['flags'] & 16) or 'Bold' in s['font'],
                                 'Courier' in s['font'],
                                 script))
    doc.close()
    raise RuntimeError("%r heading never reached in %s" % (SEC_END, pdf))


def runs_to_html(runs, bold=True):
    """Render runs as HTML. bold=False where the paragraph is already bold."""
    out = []
    for text, is_bold, is_mono, script in runs:
        if not text:
            continue
        s = esc(text)
        if script:
            s = '<%s>%s</%s>' % (script, s, script)
        if is_mono:
            s = '<span class="mono">%s</span>' % s
        if is_bold and bold:
            s = '<b>%s</b>' % s
        out.append(s)
    return ''.join(out)


def split_runs(runs, offset):
    """Split a run list at a character offset into the concatenated text."""
    head, tail, pos = [], [], 0
    for run in runs:
        text = run[0]
        end = pos + len(text)
        if end <= offset:
            head.append(run)
        elif pos >= offset:
            tail.append(run)
        else:                                     # the split falls inside a run
            head.append((text[:offset - pos],) + run[1:])
            tail.append((text[offset - pos:],) + run[1:])
        pos = end
    return head, tail


def trim_runs(runs):
    """Drop leading/trailing whitespace, including whitespace held inside a
    styled run (the panel marker is bolded together with the space before it)."""
    runs = [r for r in runs if r[0]]
    while runs and not runs[0][0].strip():
        runs = runs[1:]
    while runs and not runs[-1][0].strip():
        runs = runs[:-1]
    if not runs:
        return []
    runs[0] = (runs[0][0].lstrip(),) + runs[0][1:]
    runs[-1] = (runs[-1][0].rstrip(),) + runs[-1][1:]
    return runs


def extract_legends(pdf):
    """{fid: (title_html, legend_html)}."""
    out, current = {}, None
    for run in legend_runs(pdf):
        match = FIG_RE.match(run[0].lstrip())
        if run[1] and match:
            current = match.group(1)
            out[current] = []
        if current:
            out[current].append(run)

    legends = {}
    for fid, runs in out.items():
        title, rest, done = [], [], False
        for run in runs:
            if done:
                rest.append(run)
            elif run[1] or not run[0].strip():    # bold, or spacing between bold runs
                title.append(run)
            else:
                done = True
                rest.append(run)
        tail = PANEL_TAIL_RE.search(''.join(r[0] for r in title))
        if tail:
            title, marker = split_runs(title, tail.start())
            rest = marker + rest
        title_html = runs_to_html(trim_runs(title), bold=False)
        if fid in TITLE_FIXES:
            old, new = TITLE_FIXES[fid]
            fixed = title_html.replace(old, new)
            print("Fig. %s title: %s" % (fid, "removed %r" % old if fixed != title_html
                                         else "%r NOT FOUND -- fixed upstream?" % old))
            title_html = fixed
        legends[fid] = (title_html, runs_to_html(trim_runs(rest)))
    return legends


# ------------------------------------------------------------------ crop ----
# Mirrors docs/thesis/pipeline/prepare_figures.py, which crops the same 19
# staging figures for the dissertation -- keep the two in sync. The difference
# is that this returns the box in PNG pixels instead of a cropped image, so it
# can be mapped onto the vector PDF.
FORCE_TOP = {'3': 185, '5': 120, '6': 170, 'S8': 135, 'S10': 206}
LEFTLABEL = {'S3', 'S5'}              # big "Fig. SN" label to white out, top-left
LEFTLABEL_BOX = (0, 0, 420, 125)      # px, in the 1836x2376 staging PNGs


def segments(mask):
    segs, start = [], None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        elif not v and start is not None:
            segs.append((start, i - 1))
            start = None
    if start is not None:
        segs.append((start, len(mask) - 1))
    return segs


def seg_width(content, s, e):
    band = content[s:e + 1]
    cols = band.any(0)
    if not cols.any():
        return 0.0
    width = content.shape[1]
    wl = np.argmax(cols)
    wr = width - 1 - np.argmax(cols[::-1])
    return (wr - wl + 1) / width


def content_box(im, remove_title=True, force_top=None):
    """(left, top, right, bottom) in pixels, title band and margins removed."""
    a = np.asarray(im.convert("RGB"))
    gray = a.mean(2)
    content = gray < 245
    H, wid = gray.shape
    rsegs = segments(content.any(1))
    if not rsegs:
        return (0, 0, wid, H)
    body = rsegs
    if force_top is not None:
        body = [(s, e) for s, e in rsegs if e >= force_top] or rsegs
        csegs = segments(content[force_top:].any(0))
        return (max(0, csegs[0][0] - 15), force_top,
                min(wid, csegs[-1][1] + 15), min(H, body[-1][1] + 15))
    if remove_title and len(rsegs) >= 2:
        MAXH, MAXGAP, MINW, TOPZONE = 0.035 * H, 0.014 * H, 0.20, 0.12 * H
        def is_titleline(s, e):
            return (e - s) < MAXH and seg_width(content, s, e) > MINW
        s0, e0 = rsegs[0]
        title_end = None
        if s0 < TOPZONE and is_titleline(s0, e0):
            k = 0
            while True:
                title_end = rsegs[k][1]
                if k + 1 >= len(rsegs):
                    break
                sn, en = rsegs[k + 1]
                if sn - rsegs[k][1] > MAXGAP or not is_titleline(sn, en):
                    break
                k += 1
        if title_end is not None and title_end < rsegs[-1][1]:
            nxt = next((s for s, e in rsegs if s > title_end), title_end + 1)
            cut = (title_end + nxt) // 2
            body = [(s, e) for s, e in rsegs if s >= cut] or rsegs
    csegs = segments(content.any(0))
    pad = 15
    return (max(0, csegs[0][0] - pad), max(0, body[0][0] - pad),
            min(wid, csegs[-1][1] + pad), min(H, body[-1][1] + pad))


def crop_rect(fid):
    """Content box of Kang_Figure_<fid> in SOURCE PDF points."""
    im = Image.open(os.path.join(STAGE, "Kang_Figure_%s.png" % fid))
    if fid in LEFTLABEL:
        im = im.convert("RGB")
        ImageDraw.Draw(im).rectangle(list(LEFTLABEL_BOX), fill="white")
        box = content_box(im, remove_title=False)
    else:
        box = content_box(im, force_top=FORCE_TOP.get(fid))
    src = fitz.open(os.path.join(STAGE, "Kang_Figure_%s.pdf" % fid))
    page = src[0].rect
    sx, sy = im.size[0] / page.width, im.size[1] / page.height
    src.close()
    return fitz.Rect(box[0] / sx, box[1] / sy, box[2] / sx, box[3] / sy)


# ------------------------------------------------------------- typesetting ---
CSS = """
@font-face {font-family: arial; src: url(arial-r.ttf);}
@font-face {font-family: arial; font-weight: bold; src: url(arial-b.ttf);}
@font-face {font-family: arial; font-style: italic; src: url(arial-i.ttf);}
@font-face {font-family: arial; font-weight: bold; font-style: italic; src: url(arial-bi.ttf);}
@font-face {font-family: mono; src: url(mono-r.ttf);}
* {font-family: arial; font-size: 12px; line-height: 1.15;}
p {margin: 0;}
p.title {font-weight: bold; text-align: left; margin-bottom: 8px;}
p.legend {text-align: justify; text-indent: 36px;}
span.mono {font-family: mono;}
sup, sub {font-size: 8px;}
"""


def archive():
    a = fitz.Archive()
    fonts = [("LiberationSans-Regular", "arial-r"), ("LiberationSans-Bold", "arial-b"),
             ("LiberationSans-Italic", "arial-i"), ("LiberationSans-BoldItalic", "arial-bi"),
             ("LiberationMono-Regular", "mono-r")]
    for src, name in fonts:
        with open("%s/%s.ttf" % (FONT_DIR, src), "rb") as fh:
            a.add(fh.read(), "%s.ttf" % name)
    return a


def caption_html(title, legend):
    return '<p class="title">%s</p><p class="legend">%s</p>' % (title, legend)


def caption_height(html, arch):
    """Height the caption needs at the text width, in points.

    Measured by trial insertion rather than with fitz.Story: Story's filled
    rect over-reports by two or three lines, which would shrink every figure
    for space the caption never uses.
    """
    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=PAGE_H)
    rect = fitz.Rect(MARGIN - INSET, MARGIN, PAGE_W - MARGIN + INSET, PAGE_H - MARGIN)
    spare = page.insert_htmlbox(rect, html, archive=arch, css=CSS)
    doc.close()
    if spare[0] < 0:
        raise RuntimeError("caption does not fit on a page at all")
    return rect.height - spare[0]


def build(fid, title, legend, arch, out_dir, fit_to_page=True):
    """Render one captioned page.

    fit_to_page  Letter page throughout; a figure too tall to sit above its
                 legend is scaled down until it fits (the dissertation's look).
    else         The figure always spans the full 6.5 in text width and the
                 PAGE grows downwards instead, so nothing is ever shrunk to
                 make room for text. Page width stays 8.5 in.
    """
    html = caption_html(title, legend)
    cap_h = caption_height(html, arch)

    clip = crop_rect(fid)
    if fit_to_page:
        avail_h = TEXT_H - IMG_AFTER - cap_h
        if avail_h <= 0:
            raise RuntimeError("Fig. %s: legend alone fills the page" % fid)
        scale = min(TEXT_W / clip.width, avail_h / clip.height)
        page_h = PAGE_H
    else:
        scale = TEXT_W / clip.width
        page_h = MARGIN + clip.height * scale + IMG_AFTER + cap_h + MARGIN
    fw, fh = clip.width * scale, clip.height * scale
    fx = MARGIN + (TEXT_W - fw) / 2                     # centered, as in the thesis
    fig_rect = fitz.Rect(fx, MARGIN, fx + fw, MARGIN + fh)

    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=page_h)
    src = fitz.open(os.path.join(STAGE, "Kang_Figure_%s.pdf" % fid))
    page.show_pdf_page(fig_rect, src, 0, clip=clip)
    if fid in LEFTLABEL:
        # the label shares rows with a subtitle that must be kept, so it cannot
        # be cropped away -- paint it out at its placed position instead
        im_w = Image.open(os.path.join(STAGE, "Kang_Figure_%s.png" % fid)).size[0]
        s = src[0].rect.width / im_w                    # px -> source pt
        lx0, ly0, lx1, ly1 = [v * s for v in LEFTLABEL_BOX]
        page.draw_rect(fitz.Rect(fig_rect.x0 + (lx0 - clip.x0) * scale,
                                 fig_rect.y0 + (ly0 - clip.y0) * scale,
                                 fig_rect.x0 + (lx1 - clip.x0) * scale,
                                 fig_rect.y0 + (ly1 - clip.y0) * scale),
                       color=None, fill=(1, 1, 1))
    src.close()

    cap_rect = fitz.Rect(MARGIN - INSET, fig_rect.y1 + IMG_AFTER,
                         PAGE_W - MARGIN + INSET, page_h - MARGIN)
    spare = page.insert_htmlbox(cap_rect, html, archive=arch, css=CSS)
    if spare[0] < 0:
        raise RuntimeError("Fig. %s: caption overflowed the page by %.1f pt"
                           % (fid, -spare[0]))

    pdf = os.path.join(out_dir, "Kang_Figure_%s_with_legend.pdf" % fid)
    doc.save(pdf, garbage=4, deflate=True)
    doc.close()
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile", pdf, pdf[:-4]],
                   check=True)
    return fig_rect, cap_h, page_h


# Two renditions of the same 19 pages: dir name -> fit_to_page.
MODES = [("Fit_to_page", True), ("Full_size", False)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(BASE, "Kang_Oldham_Figures_with_legends"))
    args = ap.parse_args()

    legends = extract_legends(LEGEND_PDF)
    missing = [f for f in FIGS if f not in legends]
    assert not missing, "legends not found: %s" % missing

    arch = archive()
    for name, fit in MODES:
        out_dir = os.path.join(args.out, name)
        os.makedirs(out_dir, exist_ok=True)
        print("\n=== %s ===" % name)
        for fid in FIGS:
            title, legend = legends[fid]
            rect, cap_h, page_h = build(fid, title, legend, arch, out_dir, fit)
            print("Fig %-4s figure %5.2f x %5.2f in   caption %5.2f in   page %5.2f in"
                  % (fid, rect.width / 72, rect.height / 72, cap_h / 72, page_h / 72))
        print("%d figures -> %s" % (len(FIGS), out_dir))


if __name__ == "__main__":
    main()
