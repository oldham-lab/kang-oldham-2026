library(tidyverse)
library(data.table)
library(openxlsx)

# Table S8: Dataset summary — bulk RNA-seq cohorts included in the consensus analysis.
#
# Input:  table_s8/v3/v3.csv
# Output: table_s8/v3/Kang_Table_S8_v3.xlsx
#         One "Data" tab with all rows, plus a Legend tab.
#
# v3 changes vs v2: the former "Data_repository" column (which held accession IDs /
# URLs) is renamed "Accession_ID_or_URL", and a new human-readable "Data_Repository"
# column is added before it (mirrors the split used in table_s1).

# ============================================================
# INPUT / OUTPUT PATHS
# ============================================================
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s8/v3/"))
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

input_csv <- file.path(save_dir, "v3.csv")


# ============================================================
# openxlsx formatting helpers  (mirrored from table_s6_7)
# ============================================================
format_data_sheet <- function(wb, sheet, df) {
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

  # Size each column to the widest of its header / values (+ padding) so nothing clips
  col_widths <- vapply(seq_len(n_cols), function(j) {
    max(nchar(names(df)[j]), nchar(as.character(df[[j]])), na.rm = TRUE)
  }, numeric(1))
  setColWidths(wb, sheet, cols = 1:n_cols, widths = col_widths + 2)
  freezePane(wb, sheet, firstRow = TRUE)
}

write_legend_section <- function(wb, sheet, title, legend_df, start_row) {
  n_cols <- ncol(legend_df)

  section_style <- createStyle(
    fontName       = "Arial",
    fontSize       = 11,
    textDecoration = "bold",
    fontColour     = "#FFFFFF",
    fgFill         = "#4472C4",
    halign         = "left",
    valign         = "center"
  )
  mergeCells(wb, sheet, cols = 1:n_cols, rows = start_row)
  writeData(wb, sheet, x = title, startRow = start_row, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet, section_style, rows = start_row, cols = 1:n_cols, gridExpand = TRUE)

  subheader_style <- createStyle(
    fontName       = "Arial",
    fontSize       = 11,
    textDecoration = "bold",
    fgFill         = "#D9E1F2",
    halign         = "left",
    valign         = "center",
    border         = "Bottom",
    borderColour   = "#4472C4",
    borderStyle    = "medium"
  )
  writeData(wb, sheet, x = legend_df, startRow = start_row + 1, startCol = 1, colNames = TRUE)
  addStyle(wb, sheet, subheader_style,
           rows = start_row + 1, cols = 1:n_cols, gridExpand = TRUE)

  body_style <- createStyle(
    fontName = "Arial",
    fontSize = 11,
    halign   = "left",
    valign   = "center",
    wrapText = TRUE
  )
  body_rows <- (start_row + 2):(start_row + 1 + nrow(legend_df))
  addStyle(wb, sheet, body_style, rows = body_rows, cols = 1:n_cols, gridExpand = TRUE)

  start_row + 1 + nrow(legend_df) + 2
}


# ============================================================
# Read data
# ============================================================
dat <- fread(input_csv, data.table = FALSE)


# ============================================================
# Build workbook
# ============================================================
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 11)

# --- Data sheet ---
addWorksheet(wb, "Data")
writeData(wb, "Data", dat, rowNames = FALSE)
format_data_sheet(wb, "Data", dat)

# --- Legend tab ---
addWorksheet(wb, "Legend")

TITLE_OFFSET <- 2  # title row + blank spacer row before legend content

# Title row
writeData(wb, "Legend",
          data.frame(x = "Table S8 | Bulk RNA-seq datasets analyzed in this study"),
          rowNames = FALSE, colNames = FALSE, startRow = 1)
mergeCells(wb, "Legend", cols = 1:2, rows = 1)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold",
                     halign = "left", valign = "center", wrapText = TRUE),
         rows = 1, cols = 1:2, gridExpand = TRUE)

col_legend <- data.frame(
  Column = c(
    "Dataset",
    "First_author_last_name",
    "Last_author_last_name",
    "Journal",
    "Year",
    "PMID",
    "Data_Repository",
    "Accession_ID_or_URL",
    "Platform",
    "No_of_samples_filtered"
  ),
  Description = c(
    "Short identifier for the dataset / cohort",
    "Last name of the first author of the primary publication",
    "Last name of the last (senior) author of the primary publication",
    "Journal in which the dataset was published",
    "Year of publication",
    "PubMed identifier of the primary publication",
    "Repository / data portal hosting the data",
    "Repository accession ID or URL where raw or processed data are available",
    "Sequencing platform used to generate the bulk RNA-seq data",
    "Number of samples retained after quality-control filtering"
  ),
  stringsAsFactors = FALSE
)

next_row <- write_legend_section(wb, "Legend",
                                 title     = "Column descriptions",
                                 legend_df = col_legend,
                                 start_row = 1 + TITLE_OFFSET)

# Format Legend sheet
setColWidths(wb, "Legend", cols = 1, widths = 28)
setColWidths(wb, "Legend", cols = 2, widths = 80)
freezePane(wb, "Legend", firstRow = TRUE)

# Put Legend first, Data second
worksheetOrder(wb) <- c(2, 1)

# ============================================================
# Save
# ============================================================
out_path <- file.path(save_dir, "Kang_Table_S8_v3.xlsx")
saveWorkbook(wb, out_path, overwrite = TRUE)
message("Saved: ", out_path)
