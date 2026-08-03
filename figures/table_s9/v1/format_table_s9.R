library(data.table)
library(openxlsx)

# ============================================================
# INPUT / OUTPUT PATHS  —  edit these before running
# ============================================================
INPUT_CSV <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s9/v1/table_s9.csv")
OUT_PATH  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s9/v1/table_s9.xlsx")


# ============================================================
# Legend definition
# ============================================================
LEGEND <- data.frame(
  Column      = c("Gene", "ensembl_id", "seed", "topmodposbc"),
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


# ============================================================
# Build workbook
# ============================================================
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 11)

# --- Data tab ---
addWorksheet(wb, "table_s9")
writeData(wb, "table_s9", df, rowNames = FALSE)

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

addStyle(wb, "table_s9", header_style,
         rows = 1, cols = 1:n_cols, gridExpand = TRUE)
addStyle(wb, "table_s9", body_style,
         rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)

setColWidths(wb, "table_s9", cols = 1:n_cols, widths = "auto")
freezePane(wb, "table_s9", firstRow = TRUE)


# --- Legend tab ---
addWorksheet(wb, "Legend")

# Section header
writeData(wb, "Legend", data.frame(x = "Column descriptions"),
          rowNames = FALSE, colNames = FALSE, startRow = 1)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     fontColour = "#FFFFFF", fgFill = "#4472C4",
                     halign = "left", valign = "center"),
         rows = 1, cols = 1:2, gridExpand = TRUE)

# Legend table
writeData(wb, "Legend", LEGEND, rowNames = FALSE, startRow = 2)

addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     fgFill = "#D9E1F2", halign = "left", valign = "center",
                     border = "Bottom", borderColour = "#4472C4", borderStyle = "medium"),
         rows = 2, cols = 1:2, gridExpand = TRUE)

addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     halign = "left", valign = "center"),
         rows = 3:(nrow(LEGEND) + 2), cols = 1, gridExpand = TRUE)

addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11,
                     halign = "left", valign = "center", wrapText = TRUE),
         rows = 3:(nrow(LEGEND) + 2), cols = 2, gridExpand = TRUE)

setColWidths(wb, "Legend", cols = 1, widths = 20)
setColWidths(wb, "Legend", cols = 2, widths = 70)
freezePane(wb, "Legend", firstRow = TRUE)


# ============================================================
# Save
# ============================================================
saveWorkbook(wb, OUT_PATH, overwrite = TRUE)
message("Saved: ", OUT_PATH)
