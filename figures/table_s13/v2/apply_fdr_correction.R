library(tools)

csv_files <- list.files(
  path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s13/v2"),
  pattern = "\\.csv$",
  full.names = TRUE
)

for (f in csv_files) {
  df <- read.csv(f, row.names = 1, check.names = FALSE)

  # FDR correction across all p-values in the table
  pvals <- as.matrix(df)
  pvals_fdr <- matrix(
    p.adjust(pvals, method = "BH"),
    nrow = nrow(pvals),
    ncol = ncol(pvals),
    dimnames = dimnames(pvals)
  )

  out_path <- file.path(
    dirname(f),
    paste0(file_path_sans_ext(basename(f)), "_FDR.csv")
  )
  write.csv(as.data.frame(pvals_fdr), file = out_path)
  cat("Written:", basename(out_path), "\n")
}
