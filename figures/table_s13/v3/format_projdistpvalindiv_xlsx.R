library(data.table)
library(openxlsx)

# ============================================================
# INPUT / OUTPUT PATHS  —  edit these before running
# ============================================================

# Directory containing all *_projdistpvalindiv.csv files produced by the runcode scripts.
INPUT_DIR <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s13/v2/")

OUT_PATH  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s13/v3/Kang_Table_S13_v3.xlsx")

# A cell is highlighted only if it is 100%/0% AND its (FDR) p-value is below this.
SIG_THRESHOLD <- 0.05


# ============================================================
# Legend definition
# ============================================================
LEGEND <- data.frame(
  Column      = c("(row index)", "(cell type columns)"),
  Description = c(
    "Module index",
    "FDR-corrected (Benjamini-Hochberg) p-value for the euclidean distance between case and control mean expression projections for that module and cell type, computed against a null distribution of 10,000 random gene sets of matching size. FDR correction was applied across all p-values within each table."
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# Helper: derive a short tab name from filename
# Excel sheet names are limited to 31 characters.
# ============================================================
make_tab_name <- function(filename) {
  tab <- sub("_projdistpvalindiv_FDR\\.csv$", "", filename)

  # Shorten comparison labels
  tab <- gsub("brainSCOPE_", "",  tab)
  tab <- gsub("SZBDMulti-Seq", "SZBD",  tab)


  tab <- gsub("allAD_vs_Con",      "allAD.vs.Con",  tab)
  tab <- gsub("Schizophrenia_vs_control", "SCZ.vs.Con", tab)

  tab <- gsub("MIT_MTC", "Liu_MTG", tab)
  tab <- gsub("MIT_PFC", "Liu_DFC", tab)
  tab <- gsub("SEAAD_MTC", "Gab_MTG", tab)
  tab <- gsub("SEAAD_PFC", "Gab_DFC", tab)


  # Shorten mod type
  tab <- gsub("_bulk_megaset", "_CTRLmods", tab)
  tab <- gsub("_rosmap", "_ADmods", tab)
  tab <- gsub("_brainseq_scz", "_SCZmods", tab)


  # Truncate to 31 chars if needed (Excel limit)
  if (nchar(tab) > 31) tab <- substr(tab, 1, 31)
  return(tab)
}


# ============================================================
# Collect all CSVs
# ============================================================
files <- list.files(INPUT_DIR, pattern = "_projdistpvalindiv_FDR\\.csv$", full.names = TRUE)
files <- files[c(9:12, 5:8, 2, 1, 4, 3)]

if (length(files) == 0) stop("No projdistpvalindiv CSVs found in: ", INPUT_DIR)

csv_entries <- lapply(files, \(f) list(tab_name = make_tab_name(basename(f)), path = f))

message("Found ", length(csv_entries), " CSV file(s).")



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
  halign   = "center",
  valign   = "center"
)
pval_style <- createStyle(
  fontName = "Arial",
  fontSize = 11,
  halign   = "center",
  valign   = "center",
  numFmt   = "0.00E+00;0.00E+00;0"
)
# Sign-percentage highlight fills, stacked on top of pval_style:
#   100% of module genes higher in the less-severe (control) group -> green
#     0% (all lower in the less-severe group)                      -> red
hl_100 <- createStyle(fgFill = "#C6EFCE")
hl_0   <- createStyle(fgFill = "#FFC7CE")


# ============================================================
# Build workbook
# ============================================================
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 11)

# --- Data tabs ---
for (entry in csv_entries) {
  tab  <- entry$tab_name
  path <- entry$path

  df <- fread(path, data.table = FALSE)

  # First column is the module index (row names written by write.csv)
  colnames(df)[1] <- "Module"

  addWorksheet(wb, tab)
  writeData(wb, tab, df, rowNames = FALSE)

  n_rows <- nrow(df)
  n_cols <- ncol(df)

  addStyle(wb, tab, header_style,
           rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  addStyle(wb, tab, body_style,
           rows = 2:(n_rows + 1), cols = 1, gridExpand = TRUE)
  addStyle(wb, tab, pval_style,
           rows = 2:(n_rows + 1), cols = 2:n_cols, gridExpand = TRUE)

  # Size each column to the widest of its header / values (+ padding) so nothing clips
  col_widths <- vapply(seq_len(n_cols), function(j) {
    max(nchar(names(df)[j]), nchar(as.character(df[[j]])), na.rm = TRUE)
  }, numeric(1))
  setColWidths(wb, tab, cols = 1:n_cols, widths = col_widths + 2)
  freezePane(wb, tab, firstRow = TRUE)

  # ── Overlay sign-percentage highlights ──────────────────────────────────
  # Matching *_pctposdiff.csv (same module rows, same cell-type columns):
  # fill cells where 100% (green) or 0% (red) of module genes are higher in
  # control. Ties/NA are neither, so they stay plain.
  pct_path <- file.path(INPUT_DIR, "pct_positive_signs",
                        paste0(sub("_projdistpvalindiv_FDR\\.csv$", "", basename(path)), "_pctposdiff.csv"))
  if (file.exists(pct_path)) {
    pct <- read.csv(pct_path, row.names = 1, check.names = FALSE)
    if (nrow(pct) == n_rows) {
      cell_types <- colnames(df)[-1]            # df col 1 is the "Module" index
      for (k in seq_along(cell_types)) {
        vals    <- pct[[cell_types[k]]]
        pvals   <- df[[cell_types[k]]]          # FDR p-values shown in this column
        ws_col  <- k + 1                        # +1 to skip the Module column
        rows100 <- which(vals == 100 & pvals < SIG_THRESHOLD)
        rows0   <- which(vals == 0   & pvals < SIG_THRESHOLD)
        if (length(rows100)) addStyle(wb, tab, hl_100, rows = rows100 + 1, cols = ws_col, gridExpand = TRUE, stack = TRUE)
        if (length(rows0))   addStyle(wb, tab, hl_0,   rows = rows0   + 1, cols = ws_col, gridExpand = TRUE, stack = TRUE)
      }
    } else {
      message("  ! pct rows (", nrow(pct), ") != tab rows (", n_rows, ") for ", tab, " - highlight skipped")
    }
  } else {
    message("  ! no pct file for ", tab, " - highlight skipped")
  }

  message("  Added tab: ", tab, "  (", n_rows, " rows x ", n_cols, " cols)")
}

# --- Legend tab ---
addWorksheet(wb, "Legend")

TITLE_OFFSET <- 2  # title row + blank spacer row before legend content

# Title row
writeData(wb, "Legend",
          data.frame(x = "Table S13 | Significant dCoPA modules by cell type"),
          rowNames = FALSE, colNames = FALSE, startRow = 1)
mergeCells(wb, "Legend", cols = 1:2, rows = 1)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold",
                     halign = "left", valign = "center", wrapText = TRUE),
         rows = 1, cols = 1:2, gridExpand = TRUE)

# Section header
writeData(wb, "Legend", data.frame(x = "Column descriptions"),
          rowNames = FALSE, colNames = FALSE, startRow = 1 + TITLE_OFFSET)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     fontColour = "#FFFFFF", fgFill = "#4472C4",
                     halign = "left", valign = "center"),
         rows = 1 + TITLE_OFFSET, cols = 1:2, gridExpand = TRUE)

# Legend table
writeData(wb, "Legend", LEGEND, rowNames = FALSE, startRow = 2 + TITLE_OFFSET)

addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     fgFill = "#D9E1F2", halign = "left", valign = "center",
                     border = "Bottom", borderColour = "#4472C4", borderStyle = "medium"),
         rows = 2 + TITLE_OFFSET, cols = 1:2, gridExpand = TRUE)

addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     halign = "left", valign = "center"),
         rows = (3 + TITLE_OFFSET):(nrow(LEGEND) + 2 + TITLE_OFFSET), cols = 1, gridExpand = TRUE)

addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11,
                     halign = "left", valign = "center", wrapText = TRUE),
         rows = (3 + TITLE_OFFSET):(nrow(LEGEND) + 2 + TITLE_OFFSET), cols = 2, gridExpand = TRUE)

setColWidths(wb, "Legend", cols = 1, widths = 20)
setColWidths(wb, "Legend", cols = 2, widths = 70)
freezePane(wb, "Legend", firstRow = TRUE)

# --- Cell highlight colour code ---
hc_start <- nrow(LEGEND) + 4 + TITLE_OFFSET   # leave a blank row after the column-description table

writeData(wb, "Legend", data.frame(x = "Cell highlight colour code"),
          rowNames = FALSE, colNames = FALSE, startRow = hc_start)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     fontColour = "#FFFFFF", fgFill = "#4472C4",
                     halign = "left", valign = "center"),
         rows = hc_start, cols = 1:2, gridExpand = TRUE)

color_code <- data.frame(
  pct  = c("", ""),
  desc = c("All module genes lower in pathological samples and FDR p < 0.05",
           "All module genes higher in pathological samples and FDR p < 0.05"),
  stringsAsFactors = FALSE
)
writeData(wb, "Legend", color_code, rowNames = FALSE, colNames = FALSE, startRow = hc_start + 1)

# Colour swatch + label in column 1, matching the in-table highlight fills
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     halign = "center", valign = "center", fgFill = "#C6EFCE",
                     border = "TopBottomLeftRight", borderColour = "#BFBFBF"),
         rows = hc_start + 1, cols = 1, gridExpand = TRUE)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     halign = "center", valign = "center", fgFill = "#FFC7CE",
                     border = "TopBottomLeftRight", borderColour = "#BFBFBF"),
         rows = hc_start + 2, cols = 1, gridExpand = TRUE)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, halign = "left", valign = "center", wrapText = TRUE),
         rows = (hc_start + 1):(hc_start + 2), cols = 2, gridExpand = TRUE)


# Legend tab first
leg_idx <- which(names(wb) == "Legend")
worksheetOrder(wb) <- c(leg_idx, setdiff(seq_along(names(wb)), leg_idx))


# ============================================================
# Save
# ============================================================
saveWorkbook(wb, OUT_PATH, overwrite = TRUE)
message("Saved: ", OUT_PATH)
