library(qs)
library(dplyr)
library(openxlsx)

# ============================================================
# Paths
# ============================================================
de_dir   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/DE_ADvsCon")
OUT_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s11/v1/Kang_Table_S11_v1.xlsx")

files <- list(
  list(path = file.path(de_dir, "MIT_DE_results_DFC/mit_dfc_DESeq2_ADvsCon_by_celltype.qs"),  dataset = "Liu_DFC",     method = "DESeq2"),
  list(path = file.path(de_dir, "MIT_DE_results_DFC/mit_dfc_edgeR_ADvsCon_by_celltype.RDS"),  dataset = "Liu_DFC",     method = "edgeR"),
  list(path = file.path(de_dir, "MIT_DE_results_MTG/mit_mtg_DESeq2_ADvsCon_by_celltype.qs"),  dataset = "Liu_MTG",     method = "DESeq2"),
  list(path = file.path(de_dir, "MIT_DE_results_MTG/mit_mtg_edgeR_ADvsCon_by_celltype.RDS"),  dataset = "Liu_MTG",     method = "edgeR"),
  list(path = file.path(de_dir, "SEA_DE_results_DFC/sea_dfc_DESeq2_ADvsCon_by_celltype.qs"),  dataset = "Gabitto_DFC", method = "DESeq2"),
  list(path = file.path(de_dir, "SEA_DE_results_DFC/sea_dfc_edgeR_ADvsCon_by_celltype.RDS"),  dataset = "Gabitto_DFC", method = "edgeR"),
  list(path = file.path(de_dir, "SEA_DE_results_MTG/sea_mtg_DESeq2_ADvsCon_by_celltype.qs"),  dataset = "Gabitto_MTG", method = "DESeq2"),
  list(path = file.path(de_dir, "SEA_DE_results_MTG/sea_mtg_edgeR_ADvsCon_by_celltype.RDS"),  dataset = "Gabitto_MTG", method = "edgeR")
)

read_file <- function(path) {
  if (grepl("\\.qs$", path)) qread(path) else readRDS(path)
}

process_method <- function(method_files) {
  bind_rows(lapply(method_files, function(f) {
    obj <- read_file(f$path)
    df <- bind_rows(lapply(names(obj), function(ct) {
      obj[[ct]]$celltype <- ct
      obj[[ct]]
    }))
    df$dataset <- f$dataset
    df$method  <- f$method
    df
  }))
}

filter_and_reorder <- function(df) {
  df <- df[!is.na(df$FDR) & df$FDR < 0.05, ]
  df <- rename(df, gene = genes)
  rest <- setdiff(colnames(df), c("dataset", "method", "celltype", "gene"))
  df[, c("dataset", "method", "celltype", "gene", rest)]
}

deseq2 <- filter_and_reorder(process_method(files[sapply(files, \(f) f$method == "DESeq2")]))
edger  <- filter_and_reorder(process_method(files[sapply(files, \(f) f$method == "edgeR")]))

by_method <- list(edgeR = edger, DESeq2 = deseq2)

# ============================================================
# Shared styles
# ============================================================
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
  halign   = "left",
  valign   = "center"
)

# ============================================================
# Build workbook
# ============================================================
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 11)

write_tab <- function(wb, tab_name, df) {
  addWorksheet(wb, tab_name)
  writeData(wb, tab_name, df, rowNames = FALSE)

  n_rows <- nrow(df)
  n_cols <- ncol(df)

  addStyle(wb, tab_name, header_style,
           rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  if (n_rows > 0) {
    addStyle(wb, tab_name, body_style,
             rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)
  }

  # Size each column to the widest of its header / values (+ padding) so nothing clips
  col_widths <- vapply(seq_len(n_cols), function(j) {
    max(nchar(names(df)[j]), nchar(as.character(df[[j]])), na.rm = TRUE)
  }, numeric(1))
  setColWidths(wb, tab_name, cols = 1:n_cols, widths = col_widths + 2)
  freezePane(wb, tab_name, firstRow = TRUE)

  message("  Added tab: ", tab_name, "  (", n_rows, " rows x ", n_cols, " cols)")
}

# One tab per dataset x method x region, in the requested order. dataset/method
# are constant within a tab (encoded in the tab name), so drop them from the table.
# Excel caps sheet names at 31 chars, so the pipe separators carry no padding
# spaces ("Gabitto | AD vs. CTRL DEseq2 | MTG" would be 34 chars and error).
tab_spec <- list(
  list(method = "edgeR",  dataset = "Gabitto_MTG", title = "Gabitto|AD vs. CTRL edgeR|MTG"),
  list(method = "DESeq2", dataset = "Gabitto_MTG", title = "Gabitto|AD vs. CTRL DEseq2|MTG"),
  list(method = "edgeR",  dataset = "Liu_MTG",     title = "Liu|AD vs. CTRL edgeR|MTG"),
  list(method = "DESeq2", dataset = "Liu_MTG",     title = "Liu|AD vs. CTRL DEseq2|MTG"),
  list(method = "edgeR",  dataset = "Gabitto_DFC", title = "Gabitto|AD vs. CTRL edgeR|DFC"),
  list(method = "DESeq2", dataset = "Gabitto_DFC", title = "Gabitto|AD vs. CTRL DEseq2|DFC"),
  list(method = "edgeR",  dataset = "Liu_DFC",     title = "Liu|AD vs. CTRL edgeR|DFC"),
  list(method = "DESeq2", dataset = "Liu_DFC",     title = "Liu|AD vs. CTRL DEseq2|DFC")
)

for (t in tab_spec) {
  sub <- by_method[[t$method]]
  sub <- sub[sub$dataset == t$dataset, ]
  sub <- sub[, setdiff(colnames(sub), c("dataset", "method")), drop = FALSE]
  write_tab(wb, t$title, sub)
}

# ============================================================
# Legend tab
# ============================================================
LEGEND <- data.frame(
  Column = c("celltype", "gene", "logFC", "logCPM", "LR", "baseMean",
             "lfcSE", "stat", "PValue", "FDR"),
  Description = c(
    "Single-nucleus RNA-seq cell type (subclass) tested",
    "Gene symbol",
    "Log2 fold change of expression in AD relative to control",
    "Average log2 counts per million across samples (edgeR)",
    "Likelihood-ratio test statistic (edgeR)",
    "Mean of normalized counts across all samples (DESeq2)",
    "Standard error of the log2 fold change estimate (DESeq2)",
    "Wald test statistic (DESeq2)",
    "Raw p-value for the AD-vs-control test",
    "Benjamini-Hochberg false discovery rate; only genes with FDR < 0.05 are reported"
  ),
  stringsAsFactors = FALSE
)

addWorksheet(wb, "Legend")

TITLE_OFFSET <- 2  # title row + blank spacer row before legend content

# Title row
writeData(wb, "Legend",
          data.frame(x = "Table S11 | Differentially expressed genes (AD vs control) by subclass"),
          rowNames = FALSE, colNames = FALSE, startRow = 1)
mergeCells(wb, "Legend", cols = 1:2, rows = 1)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold",
                     halign = "left", valign = "center", wrapText = TRUE),
         rows = 1, cols = 1:2, gridExpand = TRUE)

writeData(wb, "Legend",
          data.frame(x = paste0(
            "Each data tab holds significant DE genes (FDR < 0.05) for one dataset, ",
            "differential-expression method, and brain region, named ",
            "'Dataset | AD vs. CTRL Method | Region'. ",
            "Datasets: Liu (Liu et al.) and Gabitto (Gabitto et al. / SEA-AD). ",
            "Methods: edgeR and DESeq2. ",
            "Regions: DFC (dorsolateral prefrontal cortex) and MTG (middle temporal gyrus). ",
            "All tests compare AD against control (CTRL).")),
          rowNames = FALSE, colNames = FALSE, startRow = 1 + TITLE_OFFSET)
mergeCells(wb, "Legend", cols = 1:2, rows = 1 + TITLE_OFFSET)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11,
                     halign = "left", valign = "center", wrapText = TRUE),
         rows = 1 + TITLE_OFFSET, cols = 1:2, gridExpand = TRUE)

COL_OFFSET <- TITLE_OFFSET + 2  # leave a spacer row after the description note

writeData(wb, "Legend", data.frame(x = "Column descriptions"),
          rowNames = FALSE, colNames = FALSE, startRow = 1 + COL_OFFSET)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     fontColour = "#FFFFFF", fgFill = "#4472C4",
                     halign = "left", valign = "center"),
         rows = 1 + COL_OFFSET, cols = 1:2, gridExpand = TRUE)

writeData(wb, "Legend", LEGEND, rowNames = FALSE, startRow = 2 + COL_OFFSET)

addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     fgFill = "#D9E1F2", halign = "left", valign = "center",
                     border = "Bottom", borderColour = "#4472C4", borderStyle = "medium"),
         rows = 2 + COL_OFFSET, cols = 1:2, gridExpand = TRUE)

addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     halign = "left", valign = "center"),
         rows = (3 + COL_OFFSET):(nrow(LEGEND) + 2 + COL_OFFSET), cols = 1, gridExpand = TRUE)

addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11,
                     halign = "left", valign = "center", wrapText = TRUE),
         rows = (3 + COL_OFFSET):(nrow(LEGEND) + 2 + COL_OFFSET), cols = 2, gridExpand = TRUE)

setColWidths(wb, "Legend", cols = 1, widths = max(nchar(LEGEND$Column)) + 2)
setColWidths(wb, "Legend", cols = 2, widths = 90)
freezePane(wb, "Legend", firstRow = TRUE)
message("  Added tab: Legend")

# Legend tab first
leg_idx <- which(names(wb) == "Legend")
worksheetOrder(wb) <- c(leg_idx, setdiff(seq_along(names(wb)), leg_idx))

# ============================================================
# Save
# ============================================================
saveWorkbook(wb, OUT_PATH, overwrite = TRUE)
message("Saved: ", OUT_PATH)
