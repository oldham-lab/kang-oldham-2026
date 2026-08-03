library(data.table)
library(openxlsx)

# ============================================================
# INPUT / OUTPUT PATHS  —  edit these before running
# ============================================================
FIG7      <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7")
SUMMARY   <- list(
  DFC = file.path(FIG7, "v7.2/ad_db_summary_table_dfc.csv"),
  MTG = file.path(FIG7, "v7.2/ad_db_summary_table_mtg.csv")
)
OVERLAPS  <- list(
  DFC = file.path(FIG7, "v7.2/dfc_overlaps.csv"),
  MTG = file.path(FIG7, "v7.2/mtg_overlaps.csv")
)
# dCoPA gene list carrying Module / Module type / Direction / Celltype per gene.
PANEL_B   <- file.path(FIG7, "v7.1/panel_B_dcopa_genelist.csv")

OUT_PATH  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s15/v1/Kang_Table_S15_v1.xlsx")

MAX_WIDTH <- 90   # cap so long free-text columns (mechanism, references) stay readable


# ============================================================
# Cell-type display names (match figure 7)
# ============================================================
CT_RENAME <- c(
  "Lamp5"      = "LAMP5",
  "Lamp5 Lhx6" = "LAMP5 LHX6",
  "Pax6"       = "PAX6",
  "Pvalb"      = "PVALB",
  "Sncg"       = "SNCG",
  "Sst"        = "SST",
  "Vip"        = "VIP",
  "Sst Chodl"  = "SST CHODL"
)
rename_ct <- function(x) ifelse(x %in% names(CT_RENAME), CT_RENAME[x], x)

# dCoPA direction -> pathological framing
DIR_MAP <- c(
  "Lower in more severe"  = "Lower in pathological samples",
  "Higher in more severe" = "Higher in pathological samples"
)


# ============================================================
# Legend definition
# ============================================================
LEGEND <- data.frame(
  Column      = c("Gene", "Region", "Module index (CTRL)", "Module index (AD)",
                  "Associated subclass", "Direction", "AD_mechanism",
                  "References", "Pubmed_total_hits", "OpenTargets_AD_score"),
  Description = c(
    "Gene symbol",
    "Brain region (DFC = dorsolateral prefrontal cortex; MTG = middle temporal gyrus)",
    "dCoPA module index in the CTRL (bulk megaset) module set to which this gene belongs",
    "dCoPA module index in the AD (ROSMAP AD) module set to which this gene belongs",
    "Neuronal subclass(es) in which the gene is a significant and reproducible dCoPA gene (i.e. recovered across datasets)",
    "Direction of differential module expression in pathological samples relative to controls",
    "AD mechanism linking the gene to disease, summarised from the reviewed literature",
    "Verified PMIDs supporting the AD mechanism, or the curated database source for prefilter-database genes",
    "Total PubMed hits for the gene in an AD-specific query",
    "Open Targets Platform AD association score (disease MONDO_0004975); blank if the gene has no Open Targets AD score"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# Load dCoPA annotation and build per-gene module / subclass / direction
# ============================================================
pb <- fread(PANEL_B, data.table = FALSE)

# For each region, return a data.frame: Gene, CTRL_module, AD_module, Subclass, Direction
annotate_region <- function(region) {
  ov  <- fread(OVERLAPS[[region]], data.table = FALSE)        # reproducible (Celltype, Gene) pairs
  sub <- pb[pb$Region == region, ]
  # restrict to reproducible (gene, celltype) pairs only
  sub <- sub[paste(sub$Gene, sub$Celltype) %in% paste(ov$Gene, ov$Celltype), ]

  genes <- sort(unique(ov$Gene))
  out <- lapply(genes, function(g) {
    r  <- sub[sub$Gene == g, ]
    cm <- unique(r$Module[r$`Module type` == "Bulk megaset"])
    am <- unique(r$Module[r$`Module type` == "ROSMAP AD"])
    dr <- unique(r$Direction)
    cts <- sort(unique(rename_ct(ov$Celltype[ov$Gene == g])))
    data.frame(
      Gene        = g,
      CTRL_module = if (length(cm)) paste(sort(cm), collapse = ", ") else NA_character_,
      AD_module   = if (length(am)) paste(sort(am), collapse = ", ") else NA_character_,
      Subclass    = paste(cts, collapse = ", "),
      Direction   = if (length(dr)) paste(DIR_MAP[dr], collapse = ", ") else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}


# ============================================================
# Build one merged, reordered data.frame per region
# ============================================================
build_region_df <- function(region) {
  s   <- fread(SUMMARY[[region]], data.table = FALSE)
  ann <- annotate_region(region)

  # Keep only reproducible dCoPA genes (present in the Gabitto+Liu overlap set).
  # Summary-table genes absent from the overlap set are not reproducible and are dropped.
  dropped <- setdiff(s$Gene, ann$Gene)
  if (length(dropped))
    message("  ", region, ": dropping ", length(dropped),
            " non-reproducible gene(s): ", paste(sort(dropped), collapse = ", "))
  s <- s[s$Gene %in% ann$Gene, ]

  s   <- merge(s, ann, by = "Gene", all.x = TRUE, sort = FALSE)

  # Space out the PMID separator in References for readability ("a|b" -> "a | b")
  s$References <- gsub("|", " | ", s$References, fixed = TRUE)

  # Column order: identity -> dCoPA context -> literature review
  ordered <- c("Gene", "Region", "CTRL_module", "AD_module", "Subclass", "Direction",
               "AD_mechanism", "References", "Pubmed_total_hits",
               "OpenTargets_AD_score")
  s <- s[, ordered]
  # Friendly display headers for the new columns
  names(s)[names(s) == "CTRL_module"] <- "Module index (CTRL)"
  names(s)[names(s) == "AD_module"]   <- "Module index (AD)"
  names(s)[names(s) == "Subclass"]    <- "Associated subclass"
  s
}

regions <- c("MTG", "DFC")
region_dfs <- lapply(regions, build_region_df)
names(region_dfs) <- regions


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
  borderStyle    = "medium",
  wrapText       = TRUE
)
body_style <- createStyle(
  fontName = "Arial",
  fontSize = 11,
  halign   = "left",
  valign   = "top",
  wrapText = TRUE
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

  addStyle(wb, tab_name, header_style, rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  addStyle(wb, tab_name, body_style, rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)

  # Content-aware widths, capped so long free-text columns stay readable
  col_widths <- vapply(seq_len(n_cols), function(j) {
    max(nchar(names(df)[j]), nchar(as.character(df[[j]])), na.rm = TRUE)
  }, numeric(1))
  setColWidths(wb, tab_name, cols = 1:n_cols, widths = pmin(col_widths + 2, MAX_WIDTH))
  freezePane(wb, tab_name, firstRow = TRUE)

  message("  Added tab: ", tab_name, "  (", n_rows, " rows x ", n_cols, " cols)")
}

for (rg in regions) write_tab(wb, rg, region_dfs[[rg]])


# --- Legend tab ---
addWorksheet(wb, "Legend")

TITLE_OFFSET <- 2  # title row + blank spacer row before legend content

# Title row
writeData(wb, "Legend",
          data.frame(x = "Table S15 | AI-powered literature review for significant and reproducible dCoPA genes"),
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

setColWidths(wb, "Legend", cols = 1, widths = 24)
setColWidths(wb, "Legend", cols = 2, widths = 90)
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
