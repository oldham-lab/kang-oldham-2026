library(data.table)
library(openxlsx)

# ============================================================
# INPUT / OUTPUT PATHS
# ============================================================
IN_DIR   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s10/v1")
OUT_PATH <- file.path(IN_DIR, "Kang_Table_S10_v1.xlsx")

TITLE <- "Table S10 | REI values for all modules and meta-modules in human neocortex."

# Tabs: name shown on the sheet -> CSV file + dataset descriptor for the legend
TABS <- list(
  # --- Per-module REI tabs (Figure 3) ---
  list(sheet = "REI (Jorstad 2023)",  csv = "REI_Jorstad_2023.csv",
       desc = "Per-module REI values from Jorstad et al. 2023 (normal human dorsolateral prefrontal cortex)"),
  list(sheet = "REI (Gabitto 2024)",  csv = "REI_Gabitto_2024.csv",
       desc = "Per-module REI values from Gabitto et al. 2024 (normal human dorsolateral prefrontal cortex)"),
  list(sheet = "REI (Liu 2025)",      csv = "REI_Liu_2025.csv",
       desc = "Per-module REI values from Liu et al. 2025 (normal human dorsolateral prefrontal cortex)"),
  list(sheet = "REI (Morabito 2021)", csv = "REI_Morabito_2021.csv",
       desc = "Per-module REI values from Morabito et al. 2021 (normal human dorsolateral prefrontal cortex)"),
  list(sheet = "REI Average (4 datasets)", csv = "REI_Average.csv",
       desc = "Element-wise mean REI at each module x cell type cell, averaged across the four per-module dataset tabs above (Jorstad, Gabitto, Liu, Morabito; not the meta-module tabs). Cells absent in a dataset are skipped; a cell missing in all four datasets is left blank"),
  # --- Meta-module eigengene tabs, DFC (Figure 5 panel C) ---
  list(sheet = "Eigenmod (Jorstad DFC)", csv = "MetaMod_Eigenmod_Jorstad_DFC.csv",
       desc = "Meta-module eigengene values from Jorstad et al. 2023, dorsolateral prefrontal cortex (DFC); data underlying the Figure 5 panel C heatmap"),
  list(sheet = "Eigenmod (Gabitto DFC)", csv = "MetaMod_Eigenmod_Gabitto_DFC.csv",
       desc = "Meta-module eigengene values from Gabitto et al. 2024, dorsolateral prefrontal cortex (DFC), meta-modules ordered to match Jorstad; data underlying the Figure 5 panel C heatmap"),
  # --- Meta-module eigengene tabs, MTG (Figure S8 panel C) ---
  list(sheet = "Eigenmod (Jorstad MTG)", csv = "MetaMod_Eigenmod_Jorstad_MTG.csv",
       desc = "Meta-module eigengene values from Jorstad et al. 2023, middle temporal gyrus (MTG); data underlying the Figure S8 panel C heatmap"),
  list(sheet = "Eigenmod (Gabitto MTG)", csv = "MetaMod_Eigenmod_Gabitto_MTG.csv",
       desc = "Meta-module eigengene values from Gabitto et al. 2024, middle temporal gyrus (MTG), meta-modules ordered to match Jorstad; data underlying the Figure S8 panel C heatmap"),
  # --- Meta-module membership tabs, DFC (assignment of modules to meta-modules) ---
  list(sheet = "Membership (Jorstad DFC)", csv = "MetaMod_Membership_Jorstad_DFC.csv",
       desc = "Assignment of each module to its seed and topmodposbc meta-module in the Jorstad et al. 2023 DFC network (Figure 5 panel C)"),
  list(sheet = "Membership (Gabitto DFC)", csv = "MetaMod_Membership_Gabitto_DFC.csv",
       desc = "Assignment of each module to its seed and topmodposbc meta-module in the Gabitto et al. 2024 DFC network, meta-modules ordered to match Jorstad (Figure 5 panel C)"),
  # --- Meta-module membership tabs, MTG (assignment of modules to meta-modules) ---
  list(sheet = "Membership (Jorstad MTG)", csv = "MetaMod_Membership_Jorstad_MTG.csv",
       desc = "Assignment of each module to its seed and topmodposbc meta-module in the Jorstad et al. 2023 MTG network (Figure S8 panel C)"),
  list(sheet = "Membership (Gabitto MTG)", csv = "MetaMod_Membership_Gabitto_MTG.csv",
       desc = "Assignment of each module to its seed and topmodposbc meta-module in the Gabitto et al. 2024 MTG network, meta-modules ordered to match Jorstad (Figure S8 panel C)")
)

# ============================================================
# Legend definition
# ============================================================
LEGEND <- data.frame(
  Column = c("REI tabs", "Eigenmod tabs", "Membership tabs"),
  Description = c(
    paste0("Column 'Module #' is the module identifier (1..1016), matching the module numbering in ",
           "Figure 3. The remaining columns give the relative expression index (REI) of the module in ",
           "each cell type (mean module expression normalized to mean genome-wide expression in each ",
           "subclass) - the bar heights in the Figure 3 panel F REI barplots."),
    paste0("Column 'Meta-module #' is the meta-module identifier, matching the numbering in Figure 5 ",
           "panel C (DFC) and Figure S8 panel C (MTG); DFC has 103 meta-modules and MTG has 108, with ",
           "Gabitto ordered to match Jorstad. The remaining columns give the meta-module eigengene (PC1 ",
           "of the REIs of its constituent modules) in each cell type - the values in the Figure 5 / ",
           "Figure S8 panel C heatmaps."),
    paste0("Assignment of modules to meta-modules - the meta-module analogue of Table S9's assignment ",
           "of genes to modules. Column 'Module #' is the module identifier (1..1016, matching the REI ",
           "tabs and Figure 3). 'Meta-module # (seed)' is the meta-module the module seeds, and ",
           "'Meta-module # (topmodposbc)' the meta-module it joins with positive eigengene correlation ",
           "(kME) significant after Bonferroni correction; both use the Eigenmod-tab numbering (Figure 5 ",
           "panel C / Figure S8 panel C). A module unassigned to any meta-module is left blank.")
  ),
  stringsAsFactors = FALSE
)

TABLEGEND <- data.frame(
  Tab = vapply(TABS, function(t) t$sheet, character(1)),
  Description = vapply(TABS, function(t) t$desc, character(1)),
  stringsAsFactors = FALSE
)

# ============================================================
# Styles
# ============================================================
header_style <- createStyle(
  fontName = "Arial", fontSize = 11, textDecoration = "bold",
  halign = "center", valign = "center", fgFill = "#D9E1F2",
  border = "Bottom", borderColour = "#4472C4", borderStyle = "medium"
)
body_style <- createStyle(
  fontName = "Arial", fontSize = 11, halign = "center", valign = "center"
)

# ============================================================
# Build workbook
# ============================================================
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 11)

# --- Legend tab (first) ---
addWorksheet(wb, "Legend")
TITLE_OFFSET <- 2  # title row + blank spacer row

writeData(wb, "Legend", data.frame(x = TITLE), rowNames = FALSE, colNames = FALSE, startRow = 1)
mergeCells(wb, "Legend", cols = 1:2, rows = 1)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold",
                     halign = "left", valign = "center", wrapText = TRUE),
         rows = 1, cols = 1:2, gridExpand = TRUE)

section_header <- function(text, row){
  writeData(wb, "Legend", data.frame(x = text), rowNames = FALSE, colNames = FALSE, startRow = row)
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                       fontColour = "#FFFFFF", fgFill = "#4472C4",
                       halign = "left", valign = "center"),
           rows = row, cols = 1:2, gridExpand = TRUE)
}

write_legend_block <- function(df, start_row){
  writeData(wb, "Legend", df, rowNames = FALSE, startRow = start_row)
  # header row
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                       fgFill = "#D9E1F2", halign = "left", valign = "center",
                       border = "Bottom", borderColour = "#4472C4", borderStyle = "medium"),
           rows = start_row, cols = 1:2, gridExpand = TRUE)
  # body: col 1 bold, col 2 wrapped
  body_rows <- (start_row + 1):(start_row + nrow(df))
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                       halign = "left", valign = "center"),
           rows = body_rows, cols = 1, gridExpand = TRUE)
  addStyle(wb, "Legend",
           createStyle(fontName = "Arial", fontSize = 11,
                       halign = "left", valign = "center", wrapText = TRUE),
           rows = body_rows, cols = 2, gridExpand = TRUE)
  start_row + nrow(df)  # last used row
}

# Section 1: Tabs
tab_hdr_row <- 1 + TITLE_OFFSET
section_header("Tab descriptions", tab_hdr_row)
last_row <- write_legend_block(TABLEGEND, tab_hdr_row + 1)

# Section 2: Columns
col_hdr_row <- last_row + 2
section_header("Column descriptions", col_hdr_row)
invisible(write_legend_block(LEGEND, col_hdr_row + 1))

# Width col 1 to the widest label it must hold (section headers, both block headers, all entries)
legend_col1_text <- c("Tab descriptions", "Column descriptions",
                      "Tab", "Column", TABLEGEND$Tab, LEGEND$Column)
setColWidths(wb, "Legend", cols = 1, widths = max(nchar(legend_col1_text)) + 2)
setColWidths(wb, "Legend", cols = 2, widths = 90)
freezePane(wb, "Legend", firstRow = TRUE)

# --- Data tabs ---
for(t in TABS){
  df <- fread(file.path(IN_DIR, t$csv), data.table = FALSE, check.names = FALSE)
  addWorksheet(wb, t$sheet)
  writeData(wb, t$sheet, df, rowNames = FALSE)

  n_rows <- nrow(df); n_cols <- ncol(df)
  addStyle(wb, t$sheet, header_style, rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  addStyle(wb, t$sheet, body_style, rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)

  col_widths <- vapply(seq_len(n_cols), function(j) {
    max(nchar(names(df)[j]), nchar(format(df[[j]], trim = TRUE)), na.rm = TRUE)
  }, numeric(1))
  setColWidths(wb, t$sheet, cols = 1:n_cols, widths = col_widths + 2)
  freezePane(wb, t$sheet, firstRow = TRUE, firstCol = TRUE)
}

# ============================================================
# Save
# ============================================================
saveWorkbook(wb, OUT_PATH, overwrite = TRUE)
message("Saved: ", OUT_PATH)
