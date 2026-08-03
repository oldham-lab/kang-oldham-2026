library(data.table)
library(openxlsx)

# ============================================================
# INPUT / OUTPUT PATHS  —  edit these before running
# ============================================================

# Directory containing all *_projdistpvalindiv.csv files produced by the runcode scripts.
INPUT_DIR <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s13")

OUT_PATH  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s13/table_s13.xlsx")


# ============================================================
# Legend definition
# ============================================================
LEGEND <- data.frame(
  Column      = c("(row index)", "(cell type columns)"),
  Description = c(
    "Module index",
    "P-value for the euclidean distance between case and control mean expression projections for that module and cell type, computed against a null distribution of 10,000 random gene sets of matching size"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# Helper: derive a short tab name from filename
# Excel sheet names are limited to 31 characters.
# ============================================================
make_tab_name <- function(filename) {
  tab <- sub("_projdistpvalindiv\\.csv$", "", filename)

  # Shorten comparison labels
  tab <- gsub("allAD_vs_Con",      "allAD.Con",  tab)
  tab <- gsub("earlyAD_vs_Con",    "earAD.Con",  tab)
  tab <- gsub("lateAD_vs_earlyAD", "late.earAD", tab)
  tab <- gsub("APOE44_vs_APOE33",  "AP44.33",    tab)

  # Shorten mod type
  tab <- gsub("_bulk_megaset", "_bulk", tab)

  # Truncate to 31 chars if needed (Excel limit)
  if (nchar(tab) > 31) tab <- substr(tab, 1, 31)
  return(tab)
}


# ============================================================
# Collect all CSVs
# ============================================================
files <- list.files(INPUT_DIR, pattern = "_projdistpvalindiv\\.csv$", full.names = TRUE)

if (length(files) == 0) stop("No projdistpvalindiv CSVs found in: ", INPUT_DIR)

csv_entries <- lapply(sort(files), \(f) list(tab_name = make_tab_name(basename(f)), path = f))

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
           rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)

  setColWidths(wb, tab, cols = 1:n_cols, widths = "auto")
  freezePane(wb, tab, firstRow = TRUE)

  message("  Added tab: ", tab, "  (", n_rows, " rows x ", n_cols, " cols)")
}

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
