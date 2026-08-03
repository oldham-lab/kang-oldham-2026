library(openxlsx)
library(data.table)

# ============================================================
# Paths
# ============================================================
# Resolve the directory containing this script (works via Rscript, source(), and RStudio)
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}
SAVE_DIR <- get_script_dir()

TABLE_S1_CSV  <- file.path(SAVE_DIR, "table_s1.csv")
TABLE_S1B_CSV <- file.path(SAVE_DIR, "table_s1b.csv")
LEGEND_CSV    <- file.path(SAVE_DIR, "legend.csv")
OUT_PATH      <- file.path(SAVE_DIR, "table_s1.xlsx")

# ============================================================
# Load data
# ============================================================
s1  <- as.data.frame(fread(TABLE_S1_CSV))
s1b <- as.data.frame(fread(TABLE_S1B_CSV))

# Restore original column names for table_s1b
colnames(s1b) <- c("Dataset", "Platform",
                   "# Cases", "# Controls",
                   "# Nuclei (Cases)", "# Nuclei (Controls)",
                   "Median # UMIs per Nucleus",
                   "Pubmed ID",
                   "Data repository", "Accession ID")

# Read legend exactly as CSV (two columns, no header).
# Force sep="," — fread's auto-detection otherwise picks space because the
# title row contains many spaces, which mangles the legend into many columns.
leg <- as.data.frame(fread(LEGEND_CSV, header = FALSE, fill = TRUE, sep = ","))
leg[is.na(leg)] <- ""

# ============================================================
# Style helpers (matching export_DE_to_xlsx_v2.2.R)
# ============================================================
fit_col_widths <- function(df, padding = 2) {
  sapply(seq_along(df), function(i) {
    max(nchar(c(colnames(df)[i], as.character(df[[i]]))), na.rm = TRUE) + padding
  })
}

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

  setColWidths(wb, sheet, cols = 1:n_cols, widths = fit_col_widths(df))
  freezePane(wb, sheet, firstRow = TRUE)
}

format_legend_sheet <- function(wb, sheet, df) {
  n_rows <- nrow(df)
  n_cols <- ncol(df)

  body_style <- createStyle(
    fontName = "Arial",
    fontSize = 12,
    halign   = "left",
    valign   = "center"
  )
  header_style <- createStyle(
    fontName       = "Arial",
    fontSize       = 12,
    textDecoration = "bold",
    halign         = "left",
    valign         = "center"
  )

  writeData(wb, sheet, df, rowNames = FALSE, colNames = FALSE)

  # Section-header rows are bold (title, blank separators, and the
  # "Region"/"Platform"/"Dataset" sub-headers); everything else is plain.
  header_rows <- which(df[[1]] %in% c("", "Region", "Platform") |
                         seq_len(n_rows) == 1)
  body_rows <- setdiff(seq_len(n_rows), header_rows)

  if (length(body_rows))
    addStyle(wb, sheet, body_style,
             rows = body_rows, cols = 1:n_cols, gridExpand = TRUE)
  if (length(header_rows))
    addStyle(wb, sheet, header_style,
             rows = header_rows, cols = 1:n_cols, gridExpand = TRUE)

  setColWidths(wb, sheet, cols = 1:n_cols, widths = c(22, 32))
}

# ============================================================
# Build workbook
# ============================================================
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 11)

addWorksheet(wb, "Legend")
format_legend_sheet(wb, "Legend", leg)

addWorksheet(wb, "Table S1a")
writeData(wb, "Table S1a", s1, rowNames = FALSE)
format_sheet(wb, "Table S1a", s1)

addWorksheet(wb, "Table S1b")
writeData(wb, "Table S1b", s1b, rowNames = FALSE)
format_sheet(wb, "Table S1b", s1b)

saveWorkbook(wb, OUT_PATH, overwrite = TRUE)
message("Saved: ", OUT_PATH)
