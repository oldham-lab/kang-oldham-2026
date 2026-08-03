library(pdftools)
library(magrittr)

# Convert PDFs to PNG files

#wd <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap")
wd <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/scz")
wdfiles <- list.files(wd, recursive = T, full.names = T) %>%
  grep("pdf", ., value = T)
wdfilesout <- gsub(".pdf", ".png", wdfiles)

for(i in seq_along(wdfiles)){
  pdf_convert(wdfiles[i], filenames = wdfilesout[i], dpi = 320, verbose = T)
}
