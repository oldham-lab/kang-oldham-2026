suppressMessages(library(openxlsx))
wb <- "Kang_Table_S10_v1.xlsx"
sheets <- getSheetNames(wb)
cat("Sheets:", paste(sheets, collapse=" | "), "\n\n")
for(s in sheets){
  d <- read.xlsx(wb, sheet = s, colNames = FALSE)
  cat(sprintf("  %-14s %d rows x %d cols\n", s, nrow(d), ncol(d)))
}
cat("\n--- Legend ---\n")
lg <- read.xlsx(wb, sheet="Legend", colNames=FALSE)
for(i in 1:nrow(lg)){
  row <- lg[i, !is.na(lg[i,]), drop=TRUE]
  if(length(row)) cat(i, ": ", paste(unlist(row), collapse=" || "), "\n", sep="")
}
cat("\n--- REI (Jorstad 2023): header + module 354 (first 6 cols) ---\n")
j <- read.xlsx(wb, sheet="REI (Jorstad 2023)")
print(colnames(j)[1:6]); print(j[354, 1:6])
cat("\n--- REI Average (4 datasets): module 354 (first 6 cols) ---\n")
a <- read.xlsx(wb, sheet="REI Average (4 datasets)")
print(a[354, 1:6])
