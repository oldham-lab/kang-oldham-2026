#!/usr/bin/env Rscript
# Fig 5 panel A, step 1/2: regenerate the 8 REI-schematic demo matrices and
# emit one tidy CSV of per-cell fill colours.
#
# This mirrors the data + colour definitions in fig_5_v3.R (lines 17-32). The
# matrices are illustrative random data (fixed seed), so the schematic is fully
# deterministic. Layout/labels/arrows are drawn by assemble_panel_A.py.
suppressPackageStartupMessages(library(circlize))

set.seed(23)
rei1 <- matrix(runif(20, 0, 1), nrow = 5, ncol = 4)
rei2 <- matrix(runif(20, 0, 1), nrow = 5, ncol = 4)

rei1c <- cor(rei1);    rei2c <- cor(rei2)      # celltype x celltype (4x4)
rei1m <- cor(t(rei1)); rei2m <- cor(t(rei2))   # module x module    (5x5)
reip1 <- pmin(rei1c, rei2c)                     # celltype consensus
reip2 <- pmin(rei1m, rei2m)                     # module consensus

h_col  <- colorRamp2(c(-1, 1), c("white", "grey"))     # consensus (grey)
h_col1 <- colorRamp2(c(0, 1),  c("white", "#420D09"))  # Jorstad   (red)
h_col2 <- colorRamp2(c(0, 1),  c("white", "#151B54"))  # Gabitto   (blue)

panels <- list(
  A1 = list(m = rei1,  col = h_col1), A2 = list(m = rei2,  col = h_col2),
  A3 = list(m = rei1c, col = h_col1), A4 = list(m = rei2c, col = h_col2),
  A5 = list(m = rei1m, col = h_col1), A6 = list(m = rei2m, col = h_col2),
  A7 = list(m = reip1, col = h_col),  A8 = list(m = reip2, col = h_col)
)

out <- do.call(rbind, lapply(names(panels), function(nm) {
  m <- panels[[nm]]$m; col <- panels[[nm]]$col
  nr <- nrow(m); nc <- ncol(m)
  g <- expand.grid(row = seq_len(nr), col = seq_len(nc))
  data.frame(panel = nm, row = g$row, col = g$col, nrow = nr, ncol = nc,
             hex = substr(col(m[cbind(g$row, g$col)]), 1, 7),  # drop alpha
             stringsAsFactors = FALSE)
}))

# write next to this script
args <- commandArgs(FALSE)
sd <- dirname(sub("^--file=", "", args[grep("^--file=", args)][1]))
if (length(sd) == 0 || is.na(sd) || sd == "") sd <- "."
fp <- file.path(sd, "panelA_cells.csv")
write.csv(out, fp, row.names = FALSE)
cat("wrote", fp, "(", nrow(out), "cells )\n")
