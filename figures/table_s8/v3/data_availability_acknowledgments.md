# SUPERSEDED — do not edit or use

This file was an early, Table-S8-only draft of the Data Availability / Acknowledgements language,
written when the required statements were still unknown (it carried `[INSERT]` and `[LIKE THIS]`
placeholders throughout).

**It has been superseded as of 2026-07-29.** The canonical documents are:

| Document | Contents |
|---|---|
| `docs/paper_draft/acknowledgements_draft.md` | The retrieved per-resource statements, their authoritative sources, and the record of every fix applied |
| `docs/paper_draft/acknowledgements_open_items.md` | Remaining action list, with the supporting evidence for each finding |

The finished text itself now lives in the manuscript:
`docs/paper_draft/Kang&Oldham_final_review_kang_edit.docx` (ACKNOWLEDGMENTS and DATA AVAILABILITY).

## Why this file should not be consulted

Every placeholder in it has since been filled from the authoritative source, and two of its
assumptions turned out to be wrong in ways that could mislead:

- It grouped **BrainSeq under PsychENCODE** and flagged that as its least-certain assumption. The
  grouping is in fact **correct** — Synapse tags every BrainSeq fastq `consortium = PEC`,
  `study = LIBD_szControl` — but BrainSeq *additionally* requires two specific citations (Jaffe et al.
  2018; *BrainSeq: A Human Brain Genomics Consortium*, Neuron 88:1078–1083, 2015), which this draft
  does not mention.
- It assumed the required statements were all behind access walls. Most were not: the AD Knowledge
  Portal Studies table (`syn17083367`), the PsychENCODE shared statement (`syn24240356`), the CMC
  Data Terms of Use wiki (`syn2759792/wiki/197282`) and the AnVIL citation policy are all public.

It also predates several material findings — five stale accessions, the CommonMind and PsychENCODE
migrations to the NIMH Data Archive, and the resolution of reviewer comments C153/C154/C155.

Retained rather than deleted so the earlier reasoning stays on the record.
