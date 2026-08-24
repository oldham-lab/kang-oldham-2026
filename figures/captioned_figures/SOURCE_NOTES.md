# Where the captions come from

`build_captioned_figures.py` copies each figure's number, title and legend **verbatim** from

```
/home/gugene/docs/Kang&Oldham_08.18.26.pdf
```

specifically its `FIGURE LEGENDS` section (Fig. 1–8) and its `EXTENDED DATA FIGURE LEGENDS`
section (Fig. S1–S11), reading up to the `EXTENDED DATA TABLE LEGENDS` heading. All 19 legends
are present there, and nothing is paraphrased, re-wrapped or re-punctuated — the text, the
superscript reference numbers, the `log2` subscripts, the bold panel markers and the monospaced
tool names (`edgeR`, `DESeq2`) are reproduced as the manuscript sets them.

This replaced an earlier version that read the legends out of
`/home/gugene/docs/paper_draft/Kang&Oldham_final_review_kang_edit.docx`. **The two do not
agree** — see "What changed" below.

## Extraction details, and one trap

- **Line and page numbers are filtered by POSITION, not by font.** The manuscript is
  line-numbered; those numbers sit at `x0` ≈ 28–35, well left of the 72 pt text block, and the
  page number sits below `y` = 730. Everything kept is `x0 >= 60` and `y1 <= 730`.

  The obvious filter — drop spans set in `Aptos`, the line-number font — is **wrong**. A few
  characters of real legend text are also set in Aptos: the `).` that closes the CoPA Cabana URL
  at the end of Fig. 7, Fig. 8 and Fig. S10. Filtering by font truncates those three legends
  mid-URL and the loss is easy to miss, because what remains still reads like a finished sentence.

- **Superscript vs subscript is a baseline shift, not a size drop.** Both are set at 7.92 pt
  against 12 pt body text, so size alone cannot tell `studies^29,30` from `log_2`. Each run is
  compared against the dominant baseline of its own line.

- **Spans are concatenated verbatim.** The PDF emits explicit space characters (often as their
  own spans), so no gap-based space insertion is needed, and line joins rely on the trailing space
  the source already carries. That is what makes hyphenated line breaks (`e-` / `f)` → `e-f)`,
  `studies39-` / `44` → `studies39-44`) come out right.

- The title is the run of leading **bold** text; a bolded panel marker at its end (`a-b)`) is
  handed back to the legend, matching the dissertation build's rule.

## What changed vs. the paper_draft docx

14 of 19 legends are character-identical. The five that differ:

| Figure | Difference |
|---|---|
| Fig. 3 | citation renumbered, `studies38-44` → `studies39-44` |
| Fig. 7 | **adds** "Interactive dot plots (a-b) are available on the CoPA Cabana web site (https://oldhamlab.shinyapps.io/copacabana/)." |
| Fig. 8 | same added sentence; citation `cohort41` → `cohort42` |
| Fig. S10 | same added sentence |
| Fig. S1 | title **gains** a literal `(ref).` placeholder |

So this PDF is ahead of the docx in prose, but behind it on one point: the `(ref).` placeholder in
Fig. S1's title was deleted from the docx on 2026-08-17 and is still present here.

**`(ref).` is the one and only departure from verbatim copying.** It is stripped by `TITLE_FIXES`
in the build script, which prints whether the substitution actually fired — so if the placeholder
is ever removed from the source PDF, the build says `NOT FOUND` and the entry can be deleted.
Everything else in every legend is reproduced exactly as the manuscript sets it.

## Two renditions

Both are built from the same captions in one run, into sibling folders:

| Folder | Page | Figure |
|---|---|---|
| `Fit_to_page/` | always 8.5 × 11 in | scaled down when needed so figure + legend share one letter page (the dissertation's look). Fig. 6 lands at 4.02 in wide, Fig. S3/S5 at ~4.3 in. |
| `Full_size/` | 8.5 in wide, height grows to fit (5.78 – 14.42 in) | always the full 6.5 in text width; never shrunk to make room for text |

They upload to `box:Kang_Oldham_Figures/With_legends_fit_to_page/` and
`…/With_legends_full_size/` respectively.

## Figure images and legends disagree on three titles

The title baked into the top of the figure image is cropped away and replaced by the legend title,
so each output page is self-consistent. But the two sources do not say the same thing for three
figures, which is worth knowing when reconciling the manuscript:

| Figure | Title baked into the image | Title in the legends section |
|---|---|---|
| Fig. 1 | …marker genes in human **neocortex** vary… | …marker genes in human **frontal cortex** vary… |
| Fig. S2 | **Unique subclass** marker genes in human **MTG** vary… | **Cell type** marker genes in human **temporal cortex** vary… |
| Fig. S4 | Genome-wide expression variation in **human DFC (Gabitto et al.)**… | Genome-wide expression variation in **pseudobulked snRNA-seq data from human frontal cortex**… |

The captioned outputs use the legend-section wording in all three cases.

## Verification

The build is checked by re-extracting the text of each rendered page below its caption title and
comparing it, whitespace-normalized, against the legend text pulled from the source PDF. All 19
match exactly; no page retains a duplicate title band; every legend that contains a superscript
renders one.
