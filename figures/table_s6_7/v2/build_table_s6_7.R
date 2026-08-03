# Rebuild Tables S6 & S7 (.xlsx) from the backing CSVs — no model rerun.
#
# The per-gene modeling results live in s6.csv (Jorstad) and s7.csv (SEAAD/Gabitto),
# with columns: Gene, Donor, Platform, Region, Celltype, mean_expr_pcntile, adj_r2, rmse.
# This script only re-renders the workbooks from those CSVs, so it needs no raw
# snRNA-seq data. To regenerate the CSVs themselves from scratch, use
# table_6_7_v2.1.R (which reruns the pseudobulk modeling).
#
# v2 note: the SEAAD (S7) Region is DLPFC (Brodmann A9 = dorsolateral prefrontal
# cortex). The v1 workbook mislabeled it "MTG"; that is corrected here.

suppressMessages({library(openxlsx); library(dplyr); library(data.table)})

save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s6_7/v2")
SIG_FIGS <- 4

round_sig_df <- function(df, n = SIG_FIGS) {
  df[] <- lapply(df, function(col) if (is.numeric(col)) signif(col, n) else col)
  df
}

make_tab_name <- function(donor, platform, region, celltype) {
  substr(paste(donor, platform, region, celltype, sep = " | "), 1, 31)
}

format_data_sheet <- function(wb, sheet, df) {
  n_rows <- nrow(df); n_cols <- ncol(df)
  header_style <- createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                              halign = "center", valign = "center", fgFill = "#D9E1F2",
                              border = "Bottom", borderColour = "#4472C4", borderStyle = "medium")
  body_style <- createStyle(fontName = "Arial", fontSize = 11, halign = "center", valign = "center")
  addStyle(wb, sheet, header_style, rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  if (n_rows > 0) addStyle(wb, sheet, body_style, rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)
  col_widths <- vapply(seq_len(n_cols), function(j)
    max(nchar(names(df)[j]), nchar(as.character(df[[j]])), na.rm = TRUE), numeric(1))
  setColWidths(wb, sheet, cols = 1:n_cols, widths = col_widths + 2)
  freezePane(wb, sheet, firstRow = TRUE)
}

write_legend_section <- function(wb, sheet, title, legend_df, start_row) {
  n_cols <- ncol(legend_df)
  section_style <- createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                               fontColour = "#FFFFFF", fgFill = "#4472C4", halign = "left", valign = "center")
  mergeCells(wb, sheet, cols = 1:n_cols, rows = start_row)
  writeData(wb, sheet, x = title, startRow = start_row, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet, section_style, rows = start_row, cols = 1:n_cols, gridExpand = TRUE)
  subheader_style <- createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                                 fgFill = "#D9E1F2", halign = "left", valign = "center",
                                 border = "Bottom", borderColour = "#4472C4", borderStyle = "medium")
  writeData(wb, sheet, x = legend_df, startRow = start_row + 1, startCol = 1, colNames = TRUE)
  addStyle(wb, sheet, subheader_style, rows = start_row + 1, cols = 1:n_cols, gridExpand = TRUE)
  body_style <- createStyle(fontName = "Arial", fontSize = 11, halign = "left", valign = "center", wrapText = TRUE)
  body_rows <- (start_row + 2):(start_row + 1 + nrow(legend_df))
  addStyle(wb, sheet, body_style, rows = body_rows, cols = 1:n_cols, gridExpand = TRUE)
  start_row + 1 + nrow(legend_df) + 2
}

build_workbook <- function(data, filename, dataset_label, title) {
  wb <- createWorkbook()
  modifyBaseFont(wb, fontName = "Arial", fontSize = 11)
  combos <- data |> distinct(Donor, Platform, Region, Celltype) |>
    arrange(Region, Platform, Donor, Celltype)
  for (k in seq_len(nrow(combos))) {
    combo  <- combos[k, ]
    tab_df <- data |>
      dplyr::filter(Donor == combo$Donor, Platform == combo$Platform,
                    Region == combo$Region, Celltype == combo$Celltype) |>
      round_sig_df()
    tab_name <- make_tab_name(combo$Donor, combo$Platform, combo$Region, combo$Celltype)
    addWorksheet(wb, tab_name)
    writeData(wb, tab_name, tab_df, rowNames = FALSE)
    format_data_sheet(wb, tab_name, tab_df)
  }
  addWorksheet(wb, "Legend")
  TITLE_OFFSET <- 2
  writeData(wb, "Legend", data.frame(x = title), rowNames = FALSE, colNames = FALSE, startRow = 1)
  mergeCells(wb, "Legend", cols = 1:2, rows = 1)
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold",
                       halign = "left", valign = "center", wrapText = TRUE),
           rows = 1, cols = 1:2, gridExpand = TRUE)
  col_legend <- data.frame(
    Column = c("Gene", "Donor", "Platform", "Region", "Celltype",
               "mean_expr_pcntile", "adj_r2", "rmse"),
    Description = c(
      "HGNC gene symbol", "Donor ID",
      "Sequencing platform used to generate the data (Cv3 = 10x Chromium v3; SSv4 = SMART-seq v4)",
      "Brain region from which nuclei were collected",
      "Level of cell type annotation used in the model (subclass or supertype)",
      "Mean expression percentile of the gene across all cells; higher values indicate more highly expressed genes",
      "Adjusted R-squared from the gene expression model; reflects the proportion of variance in expression explained by cell type composition",
      "Root mean square error from the gene expression model; reflects the average prediction error in expression units"),
    stringsAsFactors = FALSE)
  overview <- data |> distinct(Donor, Platform, Region, Celltype) |>
    arrange(Region, Platform, Donor, Celltype)
  next_row <- write_legend_section(wb, "Legend", "Column descriptions", col_legend, 1 + TITLE_OFFSET)
  next_row <- write_legend_section(wb, "Legend",
                                   paste0("Combinations present — ", dataset_label), overview, next_row)
  writeData(wb, "Legend",
            x = "Note: Each tab is labelled as Donor | Platform | Region | Celltype. Tab names are truncated to 31 characters (Excel limit).",
            startRow = next_row, startCol = 1, colNames = FALSE)
  addStyle(wb, "Legend", createStyle(fontName = "Arial", fontSize = 10, textDecoration = "italic"),
           rows = next_row, cols = 1:2, gridExpand = TRUE)
  setColWidths(wb, "Legend", cols = 1, widths = 25)
  setColWidths(wb, "Legend", cols = 2, widths = 80)
  freezePane(wb, "Legend", firstRow = TRUE)
  leg_idx <- which(names(wb) == "Legend")
  worksheetOrder(wb) <- c(leg_idx, setdiff(seq_along(names(wb)), leg_idx))
  saveWorkbook(wb, file.path(save_dir, filename), overwrite = TRUE)
  message("Saved: ", file.path(save_dir, filename))
}

s6 <- fread(file.path(save_dir, "s6.csv"), data.table = FALSE)
s7 <- fread(file.path(save_dir, "s7.csv"), data.table = FALSE)

build_workbook(s6, "s6_jorstad.xlsx", "Jorstad et al.",
               "Table S6 | Pseudobulk modeling of gene expression as a function of cell-type abundance (Jorstad et al.)")
build_workbook(s7, "s7_seaad2024.xlsx", "SEAAD 2024 (Gabitto et al.)",
               "Table S7: Pseudobulk modeling of gene expression as a function of cell-type abundance (Gabitto et al.)")
