# actually fixed gene column entrez id issue

library(qs)
library(openxlsx)
library(tidyverse)

# ============================================================
# INPUT / OUTPUT PATHS  —  edit these before running
# ============================================================

# edgeR results (.qs files — named lists of per-celltype data frames)
EDGER_JORSTAD_DFC  <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/DE/edgeR/jorstad/edgeR_jorstad_DFC_subclass_matchedGenes.qs")
EDGER_JORSTAD_MTG  <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/DE/edgeR/jorstad/edgeR_jorstad_MTG_subclass_matchedGenes.qs")
EDGER_SEAAD_DFC    <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/DE/edgeR/seaad/edgeR_seaad_DFC_subclass_matchedGenes.qs")
EDGER_SEAAD_MTG    <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/DE/edgeR/seaad/edgeR_seaad_MTG_subclass_matchedGenes.qs")

# DESeq2 results (.qs files — named lists of per-celltype data frames)
DESEQ2_JORSTAD_DFC <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/DE/DESeq2/jorstad/DESeq2_jorstad_DFC_subclass_matchedGenes.qs")
DESEQ2_JORSTAD_MTG <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/DE/DESeq2/jorstad/DESeq2_jorstad_MTG_subclass_matchedGenes.qs")
DESEQ2_SEAAD_DFC   <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/DE/DESeq2/seaad/DESeq2_seaad_DFC_subclass_matchedGenes.qs")
DESEQ2_SEAAD_MTG   <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/DE/DESeq2/seaad/DESeq2_seaad_MTG_subclass_matchedGenes.qs")

# Output directory
version <- "v2"
OUT_DIR <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s2_3_4_5/"), version)

# Adjusted p-value threshold for significance filtering
PADJ_THRESH <- 0.05

# Maximum significant figures for decimal columns
SIG_FIGS <- 4


# ============================================================
# Column legend definitions
# ============================================================
EDGER_LEGEND <- data.frame(
  Column      = c("celltype", "Gene", "logFC", "logCPM", "LR", "PValue", "FDR"),
  Description = c(
    "Cell subclass for which DE was tested (one-vs-all contrast)",
    "Gene symbol",
    "Log2 fold change: positive = higher in this celltype vs all others",
    "Log2 counts per million (average expression across all samples)",
    "Likelihood ratio statistic from the generalized linear model test",
    "Nominal p-value from the likelihood ratio test",
    "False discovery rate-adjusted p-value (Benjamini-Hochberg); rows are filtered to FDR < threshold"
  ),
  stringsAsFactors = FALSE
)

DESEQ2_LEGEND <- data.frame(
  Column      = c("celltype", "Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
  Description = c(
    "Cell subclass for which DE was tested (one-vs-all contrast)",
    "Gene symbol",
    "Mean normalized count across all samples",
    "Log2 fold change: positive = higher in this celltype vs all others",
    "Standard error of the log2 fold change estimate",
    "Wald or LRT statistic",
    "Nominal p-value",
    "Adjusted p-value (Benjamini-Hochberg); rows are filtered to padj < threshold"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# Helper functions
# ============================================================

round_sig_df <- function(df, n = SIG_FIGS) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) signif(col, n) else col
  })
  df
}

# Compile significant results from a named list of per-celltype data frames.
# - For edgeR: gene symbols are in the "genes" column; rename it to "Gene".
# - For DESeq2: gene symbols are in rownames; promote them to a "Gene" column.
compile_de_results <- function(reslist, padj_col) {
  bind_rows(lapply(names(reslist), function(ct) {
    df <- reslist[[ct]]

    # edgeR: "genes" column already contains gene symbols — just rename it
    if ("genes" %in% colnames(df)) {
      colnames(df)[colnames(df) == "genes"] <- "Gene"
    }

    # DESeq2: gene symbols are in rownames — promote to a "Gene" column
    if (!"Gene" %in% colnames(df)) {
      df <- tibble::rownames_to_column(df, var = "Gene")
    }

    df <- df[!is.na(df[[padj_col]]) & df[[padj_col]] < PADJ_THRESH, ]
    df$celltype <- ct

    # Final column order: celltype, Gene, then all remaining columns
    other_cols <- setdiff(colnames(df), c("celltype", "Gene"))
    df[, c("celltype", "Gene", other_cols)]
  }))
}


# ============================================================
# Formatting helpers
# ============================================================

format_sheet <- function(wb, sheet, df) {
  n_rows <- nrow(df)
  n_cols <- ncol(df)

  header_style <- createStyle(
    fontName       = "Arial",
    fontSize       = 11,
    textDecoration = "bold",
    halign         = "center",
    valign         = "center",
    fgFill         = "#D9E1F2",
    border         = "Bottom",
    borderColour   = "#4472C4",
    borderStyle    = "medium"
  )
  body_style <- createStyle(
    fontName = "Arial",
    fontSize = 11,
    halign   = "center",
    valign   = "center"
  )

  addStyle(wb, sheet, header_style,
           rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  if (n_rows > 0)
    addStyle(wb, sheet, body_style,
             rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)

  setColWidths(wb, sheet, cols = 1:n_cols, widths = "auto")
  freezePane(wb, sheet, firstRow = TRUE)
}

format_legend_sheet <- function(wb, sheet, df) {
  n_cols <- ncol(df)
  n_rows <- nrow(df)

  header_style <- createStyle(
    fontName       = "Arial",
    fontSize       = 11,
    textDecoration = "bold",
    halign         = "center",
    valign         = "center",
    fgFill         = "#D9E1F2",
    border         = "Bottom",
    borderColour   = "#4472C4",
    borderStyle    = "medium"
  )
  col_style <- createStyle(
    fontName       = "Arial",
    fontSize       = 11,
    textDecoration = "bold",
    halign         = "left",
    valign         = "center",
    wrapText       = TRUE
  )
  body_style <- createStyle(
    fontName = "Arial",
    fontSize = 11,
    halign   = "left",
    valign   = "center",
    wrapText = TRUE
  )

  addStyle(wb, sheet, header_style,
           rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  addStyle(wb, sheet, col_style,
           rows = 2:(n_rows + 1), cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet, body_style,
           rows = 2:(n_rows + 1), cols = 2, gridExpand = TRUE)

  setColWidths(wb, sheet, cols = 1, widths = 20)
  setColWidths(wb, sheet, cols = 2, widths = 70)
  freezePane(wb, sheet, firstRow = TRUE)
}


# ============================================================
# Workbook builder
# ============================================================

build_workbook <- function(edger_path, deseq2_path, label) {
  message("Building: ", label)

  edger_res  <- qread(edger_path)
  deseq2_res <- qread(deseq2_path)

  edger_df  <- compile_de_results(edger_res,  padj_col = "FDR")  |> round_sig_df()
  deseq2_df <- compile_de_results(deseq2_res, padj_col = "padj") |> round_sig_df()

  wb <- createWorkbook()
  modifyBaseFont(wb, fontName = "Arial", fontSize = 11)

  # edgeR tab
  addWorksheet(wb, "edgeR")
  writeData(wb, "edgeR", edger_df, rowNames = FALSE)
  format_sheet(wb, "edgeR", edger_df)

  # DESeq2 tab
  addWorksheet(wb, "DESeq2")
  writeData(wb, "DESeq2", deseq2_df, rowNames = FALSE)
  format_sheet(wb, "DESeq2", deseq2_df)

  # Legend tab
  addWorksheet(wb, "Legend")

  writeData(wb, "Legend", data.frame(x = "edgeR columns"),
            rowNames = FALSE, colNames = FALSE, startRow = 1)
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                       fgFill = "#4472C4", fontColour = "#FFFFFF"),
           rows = 1, cols = 1:2, gridExpand = TRUE)

  writeData(wb, "Legend", EDGER_LEGEND, rowNames = FALSE, startRow = 2)
  format_legend_sheet(wb, "Legend", EDGER_LEGEND)

  deseq2_start <- nrow(EDGER_LEGEND) + 4
  writeData(wb, "Legend", data.frame(x = "DESeq2 columns"),
            rowNames = FALSE, colNames = FALSE, startRow = deseq2_start)
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                       fgFill = "#4472C4", fontColour = "#FFFFFF"),
           rows = deseq2_start, cols = 1:2, gridExpand = TRUE)

  writeData(wb, "Legend", DESEQ2_LEGEND, rowNames = FALSE, startRow = deseq2_start + 1)
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                       halign = "left", valign = "center", wrapText = TRUE),
           rows = (deseq2_start + 1):(deseq2_start + nrow(DESEQ2_LEGEND)),
           cols = 1, gridExpand = TRUE)
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 11,
                       halign = "left", valign = "center", wrapText = TRUE),
           rows = (deseq2_start + 1):(deseq2_start + nrow(DESEQ2_LEGEND)),
           cols = 2, gridExpand = TRUE)

  note_row <- deseq2_start + nrow(DESEQ2_LEGEND) + 2
  writeData(wb, "Legend",
            data.frame(x = paste0("Note: All results are filtered to adjusted p-value < ", PADJ_THRESH, ".")),
            rowNames = FALSE, colNames = FALSE, startRow = note_row)
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 10, textDecoration = "italic"),
           rows = note_row, cols = 1:2, gridExpand = TRUE)

  out_path <- file.path(OUT_DIR, paste0("DE_results_", label, ".xlsx"))
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message("Saved: ", out_path)
}


# ============================================================
# Build all four workbooks
# ============================================================
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

build_workbook(EDGER_JORSTAD_DFC, DESEQ2_JORSTAD_DFC, label = "jorstad_DFC")
build_workbook(EDGER_JORSTAD_MTG, DESEQ2_JORSTAD_MTG, label = "jorstad_MTG")
build_workbook(EDGER_SEAAD_DFC,   DESEQ2_SEAAD_DFC,   label = "seaad_DFC")
build_workbook(EDGER_SEAAD_MTG,   DESEQ2_SEAAD_MTG,   label = "seaad_MTG")

message("All workbooks complete.")
