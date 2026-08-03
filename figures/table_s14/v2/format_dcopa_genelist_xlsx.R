library(data.table)
library(openxlsx)

# ============================================================
# Paths
# ============================================================
AD_CSV  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.1/panel_B_dcopa_genelist.csv")
SCZ_CSV <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v4/panel_B_dcopa_genelist.csv")
OUT_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s14/v2/Kang_Table_S14_v2.xlsx")


# ============================================================
# Load and clean
# ============================================================
clean_df <- function(path) {
  df <- fread(path, data.table = FALSE)
  df <- df[, -1]  # remove first column
  df[["Module type"]] <- gsub("^Bulk megaset$", "CTRL modules", df[["Module type"]])
  df[["Module type"]] <- gsub("^ROSMAP AD$",    "AD modules",   df[["Module type"]])
  df
}

df_ad  <- clean_df(AD_CSV)
df_scz <- clean_df(SCZ_CSV)


# ============================================================
# Cross-presence flag
# A gene is "shared" with the other disease if Region, Direction,
# and Celltype (subclass) all match (plus the gene symbol itself).
# ============================================================
make_key <- function(df) {
  paste(df[["Region"]], df[["Direction"]], df[["Celltype"]], df[["Gene"]],
        sep = "")
}
ad_keys  <- unique(make_key(df_ad))
scz_keys <- unique(make_key(df_scz))

df_ad[["Gene present in SCZ?"]] <- ifelse(make_key(df_ad)  %in% scz_keys, "Yes", "No")
df_scz[["Gene present in AD?"]] <- ifelse(make_key(df_scz) %in% ad_keys,  "Yes", "No")

# Disease-specific Direction labels (done AFTER cross-presence keys are computed, so the
# "more severe" wording still matches across tabs when building the shared-gene flag).
df_ad[["Direction"]]  <- gsub("in more severe", "in AD",  df_ad[["Direction"]])
df_scz[["Direction"]] <- gsub("in more severe", "in SCZ", df_scz[["Direction"]])


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
  addStyle(wb, tab_name, body_style,
           rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)

  # Size each column to the widest of its header / values (+ padding) so nothing clips
  col_widths <- vapply(seq_len(n_cols), function(j) {
    max(nchar(names(df)[j]), nchar(as.character(df[[j]])), na.rm = TRUE)
  }, numeric(1))
  setColWidths(wb, tab_name, cols = 1:n_cols, widths = col_widths + 2)
  freezePane(wb, tab_name, firstRow = TRUE)

  message("  Added tab: ", tab_name, "  (", n_rows, " rows x ", n_cols, " cols)")
}

write_tab(wb, "dCoPA genes (AD)",  df_ad)
write_tab(wb, "dCoPA genes (SCZ)", df_scz)

# --- Legend tab ---
LEGEND <- data.frame(
  Column = c("Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Module", "Gene",
             "Gene present in SCZ? / Gene present in AD?"),
  Description = c(
    "Source dataset for the co-expression modules",
    "Disease context of the co-expression modules",
    "Brain region of the co-expression modules",
    "Module type: CTRL modules = modules derived from control (bulk megaset) data; AD modules = modules derived from ROSMAP AD data",
    "Direction of differential module expression in disease (AD or SCZ) relative to controls",
    "Single-nucleus RNA-seq cell type in which the gene was identified as a dCoPA gene",
    "Index of the co-expression module the gene belongs to (numbered within each module type)",
    "Gene symbol",
    "Whether the gene is also a dCoPA gene in the other disease's tab (the AD tab reports presence in SCZ; the SCZ tab reports presence in AD). A gene counts as shared only if Region, Direction, and Celltype (subclass) all match between the two tabs."
  ),
  stringsAsFactors = FALSE
)

addWorksheet(wb, "Legend")

TITLE_OFFSET <- 2  # title row + blank spacer row before legend content

# Title row
writeData(wb, "Legend",
          data.frame(x = "Table S14 | Significant dCoPA genes by cell type"),
          rowNames = FALSE, colNames = FALSE, startRow = 1)
mergeCells(wb, "Legend", cols = 1:2, rows = 1)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold",
                     halign = "left", valign = "center", wrapText = TRUE),
         rows = 1, cols = 1:2, gridExpand = TRUE)

writeData(wb, "Legend", data.frame(x = "Column descriptions"),
          rowNames = FALSE, colNames = FALSE, startRow = 1 + TITLE_OFFSET)
addStyle(wb, "Legend",
         createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                     fontColour = "#FFFFFF", fgFill = "#4472C4",
                     halign = "left", valign = "center"),
         rows = 1 + TITLE_OFFSET, cols = 1:2, gridExpand = TRUE)

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

setColWidths(wb, "Legend", cols = 1, widths = max(nchar(LEGEND$Column)) + 2)
setColWidths(wb, "Legend", cols = 2, widths = 80)
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
