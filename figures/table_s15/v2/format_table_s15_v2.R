library(data.table)
library(openxlsx)

# ============================================================
# Table S15 v2 — AI-powered literature review for reproducible dCoPA genes.
# v2 adds the schizophrenia (SCZ, fig_8) gene database alongside the Alzheimer's
# (AD, fig_7) DFC + MTG databases. Sheets: Legend, AD (MTG), AD (DFC), SCZ (DFC).
# Adapted from table_s15/v1/format_table_s15.R.
# ============================================================
FIG7 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7")
FIG8 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8")
OUT_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s15/v2/Kang_Table_S15_v2.xlsx")
MAX_WIDTH <- 90

CT_RENAME <- c("Lamp5"="LAMP5","Lamp5 Lhx6"="LAMP5 LHX6","Pax6"="PAX6","Pvalb"="PVALB",
               "Sncg"="SNCG","Sst"="SST","Vip"="VIP","Sst Chodl"="SST CHODL")
rename_ct <- function(x) ifelse(x %in% names(CT_RENAME), CT_RENAME[x], x)
DIR_MAP <- c("Lower in more severe"="Lower in pathological samples",
             "Higher in more severe"="Higher in pathological samples")

# ---- Generic dCoPA annotation: per gene -> CTRL module, disease module, subclass, direction
annotate_region <- function(panel_b, overlaps, region, dis_modtype) {
  pb <- fread(panel_b, data.table = FALSE)
  ov <- fread(overlaps, data.table = FALSE)
  sub <- pb[pb$Region == region, ]
  sub <- sub[paste(sub$Gene, sub$Celltype) %in% paste(ov$Gene, ov$Celltype), ]
  genes <- sort(unique(ov$Gene))
  do.call(rbind, lapply(genes, function(g) {
    r  <- sub[sub$Gene == g, ]
    cm <- unique(r$Module[r$`Module type` == "Bulk megaset"])
    dm <- unique(r$Module[r$`Module type` == dis_modtype])
    dr <- unique(r$Direction)
    cts <- sort(unique(rename_ct(ov$Celltype[ov$Gene == g])))
    data.frame(Gene = g,
               CTRL_module = if (length(cm)) paste(sort(cm), collapse = ", ") else NA_character_,
               DIS_module  = if (length(dm)) paste(sort(dm), collapse = ", ") else NA_character_,
               Subclass    = paste(cts, collapse = ", "),
               Direction   = if (length(dr)) paste(DIR_MAP[dr], collapse = ", ") else NA_character_,
               stringsAsFactors = FALSE)
  }))
}

# ---- AD region sheet (unchanged columns from v1) ----
build_ad <- function(region) {
  s   <- fread(file.path(FIG7, "v7.2", paste0("ad_db_summary_table_", tolower(region), ".csv")), data.table = FALSE)
  ann <- annotate_region(file.path(FIG7, "v8/panel_B_dcopa_genelist.csv"),
                         file.path(FIG7, "v7.2", paste0(tolower(region), "_overlaps.csv")),
                         region, "ROSMAP AD")
  dropped <- setdiff(s$Gene, ann$Gene)
  if (length(dropped)) message("  AD ", region, ": dropping ", length(dropped), " non-reproducible: ",
                               paste(sort(dropped), collapse = ", "))
  s <- merge(s[s$Gene %in% ann$Gene, ], ann, by = "Gene", all.x = TRUE, sort = FALSE)
  s$References <- gsub("|", " | ", s$References, fixed = TRUE)
  s <- s[, c("Gene","Region","CTRL_module","DIS_module","Subclass","Direction",
             "AD_mechanism","References","Pubmed_total_hits","OpenTargets_AD_score")]
  names(s)[3:5] <- c("Module index (CTRL)","Module index (AD)","Associated subclass")
  s
}

# ---- SCZ DFC sheet (SCZ-specific columns; no Open Targets score) ----
build_scz <- function() {
  s   <- fread(file.path(FIG8, "v4/scz_db_summary_table_dfc.csv"), data.table = FALSE)
  ann <- annotate_region(file.path(FIG8, "v4/panel_B_dcopa_genelist.csv"),
                         file.path(FIG8, "v4/dfc_overlaps.csv"), "DFC", "SCZ mods")
  dropped <- setdiff(s$Gene, ann$Gene)
  if (length(dropped)) message("  SCZ DFC: dropping ", length(dropped), " non-reproducible: ",
                               paste(sort(dropped), collapse = ", "))
  s <- merge(s[s$Gene %in% ann$Gene, ], ann, by = "Gene", all.x = TRUE, sort = FALSE)
  s$References <- gsub("|", " | ", s$References, fixed = TRUE)
  s <- s[, c("Gene","Region","CTRL_module","DIS_module","Subclass","Direction",
             "SCZ_mechanism","References","Pubmed_total_hits")]
  names(s)[3:5] <- c("Module index (CTRL)","Module index (SCZ)","Associated subclass")
  s
}

sheets <- list("AD (MTG)" = build_ad("MTG"),
               "AD (DFC)" = build_ad("DFC"),
               "SCZ (DFC)" = build_scz())

# ---- Legend ----
LEGEND <- data.frame(
  Column = c("Gene","Region","Module index (CTRL)","Module index (AD) / (SCZ)",
             "Associated subclass","Direction","AD_mechanism / SCZ_mechanism",
             "References","Pubmed_total_hits","OpenTargets_AD_score"),
  Description = c(
    "Gene symbol",
    "Brain region (DFC = dorsolateral prefrontal cortex; MTG = middle temporal gyrus)",
    "dCoPA module index in the CTRL (bulk megaset) module set to which this gene belongs",
    "dCoPA module index in the disease module set (AD = ROSMAP AD modules; SCZ = brainseq SCZ modules)",
    "Neuronal subclass(es) in which the gene is a significant and reproducible dCoPA gene (recovered across datasets)",
    "Direction of differential module expression in pathological samples relative to controls",
    "Disease mechanism linking the gene, summarised from the reviewed literature (or the curated-database note for prefilter genes)",
    "Verified PMIDs supporting the mechanism, or the curated database source for prefilter-database genes",
    "Total PubMed hits for the gene in a disease-specific query (AD- or SCZ-specific per sheet)",
    "Open Targets Platform AD association score (MONDO_0004975); AD sheets only, blank if none. Not computed for the SCZ sheet."
  ),
  stringsAsFactors = FALSE)

# ---- Styles ----
header_style <- createStyle(fontName="Arial", fontSize=11, textDecoration="bold", halign="center",
                            valign="center", fgFill="#D9E1F2", border="Bottom", borderColour="#4472C4",
                            borderStyle="medium", wrapText=TRUE)
body_style <- createStyle(fontName="Arial", fontSize=11, halign="left", valign="top", wrapText=TRUE)

wb <- createWorkbook(); modifyBaseFont(wb, fontName="Arial", fontSize=11)
write_tab <- function(tab, df) {
  addWorksheet(wb, tab); writeData(wb, tab, df, rowNames=FALSE)
  nr <- nrow(df); nc <- ncol(df)
  addStyle(wb, tab, header_style, rows=1, cols=1:nc, gridExpand=TRUE)
  addStyle(wb, tab, body_style, rows=2:(nr+1), cols=1:nc, gridExpand=TRUE)
  w <- vapply(seq_len(nc), function(j) max(nchar(names(df)[j]), nchar(as.character(df[[j]])), na.rm=TRUE), numeric(1))
  setColWidths(wb, tab, cols=1:nc, widths=pmin(w+2, MAX_WIDTH)); freezePane(wb, tab, firstRow=TRUE)
  message("  Added tab: ", tab, "  (", nr, " rows x ", nc, " cols)")
}
for (nm in names(sheets)) write_tab(nm, sheets[[nm]])

# ---- Legend tab ----
addWorksheet(wb, "Legend"); TO <- 2
writeData(wb, "Legend", data.frame(x="Table S15 | AI-powered literature review for significant and reproducible dCoPA genes (AD and SCZ)"),
          rowNames=FALSE, colNames=FALSE, startRow=1)
mergeCells(wb, "Legend", cols=1:2, rows=1)
addStyle(wb, "Legend", createStyle(fontName="Arial", fontSize=14, textDecoration="bold", halign="left", valign="center", wrapText=TRUE),
         rows=1, cols=1:2, gridExpand=TRUE)
writeData(wb, "Legend", data.frame(x="Column descriptions"), rowNames=FALSE, colNames=FALSE, startRow=1+TO)
addStyle(wb, "Legend", createStyle(fontName="Arial", fontSize=11, textDecoration="bold", fontColour="#FFFFFF",
         fgFill="#4472C4", halign="left", valign="center"), rows=1+TO, cols=1:2, gridExpand=TRUE)
writeData(wb, "Legend", LEGEND, rowNames=FALSE, startRow=2+TO)
addStyle(wb, "Legend", createStyle(fontName="Arial", fontSize=11, textDecoration="bold", fgFill="#D9E1F2",
         halign="left", valign="center", border="Bottom", borderColour="#4472C4", borderStyle="medium"),
         rows=2+TO, cols=1:2, gridExpand=TRUE)
addStyle(wb, "Legend", createStyle(fontName="Arial", fontSize=11, textDecoration="bold", halign="left", valign="center"),
         rows=(3+TO):(nrow(LEGEND)+2+TO), cols=1, gridExpand=TRUE)
addStyle(wb, "Legend", createStyle(fontName="Arial", fontSize=11, halign="left", valign="center", wrapText=TRUE),
         rows=(3+TO):(nrow(LEGEND)+2+TO), cols=2, gridExpand=TRUE)
setColWidths(wb, "Legend", cols=1, widths=26); setColWidths(wb, "Legend", cols=2, widths=95)
freezePane(wb, "Legend", firstRow=TRUE); message("  Added tab: Legend")

leg <- which(names(wb) == "Legend")
worksheetOrder(wb) <- c(leg, setdiff(seq_along(names(wb)), leg))
saveWorkbook(wb, OUT_PATH, overwrite=TRUE)
message("Saved: ", OUT_PATH)
