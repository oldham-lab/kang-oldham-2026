library(openxlsx)

# ============================================================
# v7 = v6 base, with the Legend tab rebuilt in house style.
#
# v6 (Kang_Table_S1_v6.xlsx) is taken verbatim as the base: its two data
# tabs (Table S1a / Table S1b) are preserved exactly. Only the Legend tab is
# replaced -- the flat v6 legend is rebuilt with the inlined, house-style
# formatting introduced in table_s1/v2 (title banner, blank spacer row, blue
# section headers, lavender column-name row, freeze pane, Legend tab first).
# ============================================================
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}
SAVE_DIR <- get_script_dir()

BASE_XLSX <- normalizePath(file.path(SAVE_DIR, "..", "v6", "Kang_Table_S1_v6.xlsx"))
OUT_PATH  <- file.path(SAVE_DIR, "Kang_Table_S1_v7.xlsx")

TITLE <- "Table S1 | Summary of all single-nucleus RNA-seq datasets analyzed in this study"

# ============================================================
# Legend definition (inlined; v6 carried these in the workbook itself)
# ============================================================
REGION_LEGEND <- data.frame(
  Abbreviation = c("DFC", "V1", "MTG"),
  Description  = c("Dorsolateral prefrontal cortex",
                   "Primary visual cortex",
                   "Medial temporal gyrus"),
  stringsAsFactors = FALSE
)
PLATFORM_LEGEND <- data.frame(
  Abbreviation = c("Cv3", "SSv4"),
  Description  = c("10x Chromium v3", "SMART-Seq v4"),
  stringsAsFactors = FALSE
)

# House-style legend section: blue section header, lavender column-name row,
# wrapped body. Returns the next free row, leaving a one-row gap. (Matches the
# write_legend_section helper used by table_s8 / table_s6_7 / table_s1 v2.)
write_legend_section <- function(wb, sheet, title, legend_df, start_row) {
  n_cols <- ncol(legend_df)

  section_style <- createStyle(
    fontName = "Arial", fontSize = 11, textDecoration = "bold",
    fontColour = "#FFFFFF", fgFill = "#4472C4",
    halign = "left", valign = "center"
  )
  mergeCells(wb, sheet, cols = 1:n_cols, rows = start_row)
  writeData(wb, sheet, x = title, startRow = start_row, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet, section_style, rows = start_row, cols = 1:n_cols, gridExpand = TRUE)

  subheader_style <- createStyle(
    fontName = "Arial", fontSize = 11, textDecoration = "bold",
    fgFill = "#D9E1F2", halign = "left", valign = "center",
    border = "Bottom", borderColour = "#4472C4", borderStyle = "medium"
  )
  writeData(wb, sheet, x = legend_df, startRow = start_row + 1, startCol = 1, colNames = TRUE)
  addStyle(wb, sheet, subheader_style,
           rows = start_row + 1, cols = 1:n_cols, gridExpand = TRUE)

  body_style <- createStyle(
    fontName = "Arial", fontSize = 11,
    halign = "left", valign = "center", wrapText = TRUE
  )
  body_rows <- (start_row + 2):(start_row + 1 + nrow(legend_df))
  addStyle(wb, sheet, body_style, rows = body_rows, cols = 1:n_cols, gridExpand = TRUE)

  start_row + 1 + nrow(legend_df) + 2
}

format_legend_sheet <- function(wb, sheet) {
  TITLE_OFFSET <- 2  # title row + blank spacer row before legend content

  # Title row
  writeData(wb, sheet, data.frame(x = TITLE),
            rowNames = FALSE, colNames = FALSE, startRow = 1)
  mergeCells(wb, sheet, cols = 1:2, rows = 1)
  addStyle(wb, sheet,
           createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold",
                       halign = "left", valign = "center", wrapText = TRUE),
           rows = 1, cols = 1:2, gridExpand = TRUE)

  next_row <- write_legend_section(wb, sheet, "Region",
                                   REGION_LEGEND, start_row = 1 + TITLE_OFFSET)
  next_row <- write_legend_section(wb, sheet, "Platform",
                                   PLATFORM_LEGEND, start_row = next_row)

  setColWidths(wb, sheet, cols = 1, widths = 24)
  setColWidths(wb, sheet, cols = 2, widths = 70)
  freezePane(wb, sheet, firstRow = TRUE)
}

# ============================================================
# Load v6, swap the Legend tab, save v7
# ============================================================
wb <- loadWorkbook(BASE_XLSX)

removeWorksheet(wb, "Legend")
addWorksheet(wb, "Legend")
format_legend_sheet(wb, "Legend")

# Legend tab first
leg_idx <- which(names(wb) == "Legend")
worksheetOrder(wb) <- c(leg_idx, setdiff(seq_along(names(wb)), leg_idx))

saveWorkbook(wb, OUT_PATH, overwrite = TRUE)
message("Saved: ", OUT_PATH)
