library(qs)
library(dplyr)

de_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/DE")

files <- list(
  list(path = file.path(de_dir, "MIT_DE_results_DFC/mit_dfc_DESeq2_ADvsCon_by_celltype.qs"),  dataset = "Liu_DFC",     method = "DESeq2"),
  list(path = file.path(de_dir, "MIT_DE_results_DFC/mit_dfc_edgeR_ADvsCon_by_celltype.RDS"),  dataset = "Liu_DFC",     method = "edgeR"),
  list(path = file.path(de_dir, "MIT_DE_results_MTG/mit_mtg_DESeq2_ADvsCon_by_celltype.qs"),  dataset = "Liu_MTG",     method = "DESeq2"),
  list(path = file.path(de_dir, "MIT_DE_results_MTG/mit_mtg_edgeR_ADvsCon_by_celltype.RDS"),  dataset = "Liu_MTG",     method = "edgeR"),
  list(path = file.path(de_dir, "SEA_DE_results_DFC/sea_dfc_DESeq2_ADvsCon_by_celltype.qs"),  dataset = "Gabitto_DFC", method = "DESeq2"),
  list(path = file.path(de_dir, "SEA_DE_results_DFC/sea_dfc_edgeR_ADvsCon_by_celltype.RDS"),  dataset = "Gabitto_DFC", method = "edgeR"),
  list(path = file.path(de_dir, "SEA_DE_results_MTG/sea_mtg_DESeq2_ADvsCon_by_celltype.qs"),  dataset = "Gabitto_MTG", method = "DESeq2"),
  list(path = file.path(de_dir, "SEA_DE_results_MTG/sea_mtg_edgeR_ADvsCon_by_celltype.RDS"),  dataset = "Gabitto_MTG", method = "edgeR")
)

read_file <- function(path) {
  if (grepl("\\.qs$", path)) qread(path) else readRDS(path)
}

process_method <- function(method_files) {
  bind_rows(lapply(method_files, function(f) {
    obj <- read_file(f$path)
    df <- bind_rows(lapply(names(obj), function(ct) {
      obj[[ct]]$celltype <- ct
      obj[[ct]]
    }))
    df$dataset <- f$dataset
    df$method  <- f$method
    df
  }))
}

filter_and_reorder <- function(df) {
  df <- df[!is.na(df$FDR) & df$FDR < 0.05, ]
  df <- rename(df, gene = genes)
  rest <- setdiff(colnames(df), c("dataset", "method", "celltype", "gene"))
  df[, c("dataset", "method", "celltype", "gene", rest)]
}

deseq2 <- filter_and_reorder(process_method(files[sapply(files, \(f) f$method == "DESeq2")]))
edger   <- filter_and_reorder(process_method(files[sapply(files, \(f) f$method == "edgeR")]))

write.csv(deseq2, file = file.path(de_dir, "CTRLvsAD_DE_by_subclass_DESeq2.csv"), row.names = FALSE)
write.csv(edger,  file = file.path(de_dir, "CTRLvsAD_DE_by_subclass_edgeR.csv"),  row.names = FALSE)
cat("Saved", nrow(deseq2), "DESeq2 results to CTRLvsAD_DE_by_subclass_DESeq2.csv\n")
cat("Saved", nrow(edger),  "edgeR results  to CTRLvsAD_DE_by_subclass_edgeR.csv\n")
