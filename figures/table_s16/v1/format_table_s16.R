library(data.table)
library(openxlsx)

# ============================================================
# INPUT / OUTPUT PATHS  —  edit these before running
# ============================================================
INPUT_CSV <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s16/v1/table_s16.csv")
OUT_PATH  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s16/v1/Kang_Table_S16_v1.xlsx")

TITLE <- paste0("Table S16 | Assignments of all genes to all modules identified ",
                "by CoPA in frontal cortex from patients with schizophrenia")


# ============================================================
# Legend definition
# ============================================================
LEGEND <- data.frame(
  Column      = c("Gene", "ensembl_id", "Module # (seed)", "Module # (topmodposbc)"),
  Description = c(
    "Gene symbol",
    "Ensembl ID",
    "Module seed genes",
    "Modules defined by all unique genes with positive kME values that are significant after Bonferroni correction"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# Load data
# ============================================================
df <- fread(INPUT_CSV, data.table = FALSE)

# Rename module columns for display
names(df)[names(df) == "seed"]        <- "Module # (seed)"
names(df)[names(df) == "topmodposbc"] <- "Module # (topmodposbc)"


# ============================================================
# Build workbook
# ============================================================
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 11)

# --- Data tab ---
addWorksheet(wb, "table_s16")
writeData(wb, "table_s16", df, rowNames = FALSE)

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

addStyle(wb, "table_s16", header_style,
         rows = 1, cols = 1:n_cols, gridExpand = TRUE)
addStyle(wb, "table_s16", body_style,
         rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)

# Size each column to the widest of its header / values (+ padding) so nothing clips
col_widths <- vapply(seq_len(n_cols), function(j) {
  max(nchar(names(df)[j]), nchar(as.character(df[[j]])), na.rm = TRUE)
}, numeric(1))
setColWidths(wb, "table_s16", cols = 1:n_cols, widths = col_widths + 2)
freezePane(wb, "table_s16", firstRow = TRUE)


# --- Legend tab ---
addWorksheet(wb, "Legend")

TITLE_OFFSET <- 2  # title row + blank spacer row before legend content

# Title row
writeData(wb, "Legend", data.frame(x = TITLE),
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

setColWidths(wb, "Legend", cols = 1, widths = 24)
setColWidths(wb, "Legend", cols = 2, widths = 70)
freezePane(wb, "Legend", firstRow = TRUE)

# Put Legend first, table_s16 second
worksheetOrder(wb) <- c(2, 1)


# ============================================================
# Save
# ============================================================
saveWorkbook(wb, OUT_PATH, overwrite = TRUE)
message("Saved: ", OUT_PATH)
