library(dplyr)
if (!requireNamespace("gt", quietly = TRUE)) install.packages("gt")
library(gt)

# ── AD gene → reference number (ordered by strength/relevance) ────────────────
ad_gene_refs <- c(
  APOE=1,
  APP=2, PSEN1=2, PSEN2=2,
  TREM2=3,
  BIN1=4, CLU=4, CR1=4, PICALM=4,
  ABCA7=5,
  SORL1=6,
  ADAM10=7,
  FERMT2=8, PTK2B=8, INPP5D=8, MEF2C=8,
  MS4A6A=8, MS4A4A=8, EPHA1=8, CELF1=8, CD33=8, CASS4=8,
  ZCWPW1=8, SLC24A4=8, NME8=8,
  MAPT=9,
  GRN=10,
  PLCG2=11, ABI3=11,
  ACE=12
)

# ── Build one row per celltype ─────────────────────────────────────────────────
celltypes <- unique(dfc_overlaps$Celltype)
gsea_cols <- colnames(gsea_out)[4:ncol(gsea_out)]
# Map original celltype names to R-sanitised column names
ct_to_col <- setNames(make.names(celltypes), celltypes)

build_row <- function(ct) {
  genes   <- dfc_overlaps$Gene[dfc_overlaps$Celltype == ct]
  n_genes <- length(genes)

  # Top 5 gene sets for this celltype
  col_name <- ct_to_col[[ct]]
  if (col_name %in% gsea_cols) {
    pvals   <- gsea_out[[col_name]]
    top_idx <- order(pvals)[1:min(5, sum(!is.na(pvals)))]
    top_sets <- gsea_out$SetName[top_idx]
    top_pval <- pvals[top_idx]
    top_str  <- paste0(
      seq_along(top_sets), ". ", top_sets,
      " (p=", formatC(top_pval, format = "e", digits = 2), ")",
      collapse = "\n"
    )
  } else {
    top_str <- NA_character_
  }

  # AD genes: sort by ref number (most established first)
  genes_upper <- toupper(genes)
  is_ad       <- genes_upper %in% names(ad_gene_refs)
  ad_genes    <- genes[is_ad]

  if (length(ad_genes) > 0) {
    ref_nums   <- ad_gene_refs[toupper(ad_genes)]
    ord        <- order(ref_nums)
    ad_genes_s <- ad_genes[ord]
    ref_nums_s <- ref_nums[ord]
    ad_str <- paste0(
      ad_genes_s, "<sup>", ref_nums_s, "</sup>",
      collapse = ", "
    )
  } else {
    ad_str <- "-"
  }

  data.frame(
    Celltype     = ct,
    N_genes      = n_genes,
    Top_genesets = top_str,
    AD_genes     = ad_str,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(celltypes, build_row)) |>
  arrange(Celltype)

# ── Format and save ────────────────────────────────────────────────────────────
out_path <- file.path(save_dir, "panel_G_gsea_summary_table.html")

gt_table <- summary_df |>
  gt() |>
  cols_label(
    Celltype     = "Cell type",
    N_genes      = "Genes (n)",
    Top_genesets = "Top 5 gene sets (by p-value)",
    AD_genes     = "AD-associated genes"
  ) |>
  fmt_markdown(columns = c(Top_genesets, AD_genes)) |>
  cols_width(
    Celltype     ~ px(120),
    N_genes      ~ px(70),
    Top_genesets ~ px(400),
    AD_genes     ~ px(250)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_footnote(
    footnote  = "Superscripts refer to bibliography: gsea_summary_bibliography.md",
    locations = cells_column_labels(columns = AD_genes)
  ) |>
  tab_options(
    table.font.size       = 11,
    data_row.padding      = px(6),
    column_labels.padding = px(8)
  )

gtsave(gt_table, out_path)
message("Saved: ", out_path)
