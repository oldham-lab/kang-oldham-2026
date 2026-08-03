library(tidyverse)

save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/schematic/")
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

################################
# Example projections (boxplots)
################################

distlistcon <- list(list(c(0, 1), c(3, 4.5), c(1, 2)), # passes filter
                    list(c(0, 1), c(2, 2.9), c(2, 3)), # not sig
                    list(c(1.3, 1.8), c(2, 2.6), c(1.8, 2.4)), # sig but not highest
                    list(c(0.1, 1.1), c(0.1, 0.8), c(2, 2.6)) # sig and highest but not consistent
                    )
distlistAD <- list(list(c(0.1, 1.1), c(2, 3), c(0.9, 1.9)),
                   list(c(0.1, 1.1), c(2, 2.6), c(1.9, 2.9)),
                   list(c(0.1, 1.1), c(2, 2.6), c(1.8, 2.4)),   
                   list(c(0.1, 1.1), c(0.1, 1), c(1.2, 2.4))
                   )
ylimvec <- c(5.3, 3.5, 3.5, 3.5)
siglist <- list(c(4.6, 1.7, 2.3), #y_position, xmin, xmax
                c(0, 0, 0),
                c(2.5, 0.7, 1.3),
                c(3, 2.7, 3.3)    
                )

for(z in 1:4){ # four plots
  b1 <- list()
  for(i in 1:3){ # three celltypes
    conmat <- data.frame("CT" = paste0("CT", i),
                        "Disease" = c("Con"),
                        "Value" = runif(20, distlistcon[[z]][[i]][1], distlistcon[[z]][[i]][2]))

    ADmat <- data.frame("CT" = paste0("CT", i),
                        "Disease" = c("AD"),
                        "Value" = runif(20, distlistAD[[z]][[i]][1], distlistAD[[z]][[i]][2]))

    b1[[i]] <- rbind(conmat, ADmat)
  }
  b1dat <- do.call(rbind, b1)

  p <- ggplot(b1dat, aes(x = CT, y = Value, fill = Disease)) +
    theme_classic() + 
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.6) +
    labs(x = "", y = "Expression") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.2, "cm"),
          legend.position = "none"
          #legend.position = "inside",
          #legend.position.inside = c(0.8, 0.2)
          ) +
    ylim(0, ylimvec[z]) +
    scale_fill_manual(values = c("blue", "red"))

  if(z != 2){
    p <- p +
      ggsignif::geom_signif(y_position = siglist[[z]][1], xmin = siglist[[z]][2], xmax = siglist[[z]][3], annotation = "***", tip_length = 0, vjust = 0.2) 
  }
  ggsave(p, file = file.path(save_dir, paste0("p", z, ".pdf")), width = 1.5, height = 1.5)
}

# Get legend 
library(grid)
library(gridExtra) 
p <- ggplot(b1dat, aes(x = CT, y = Value, fill = Disease)) +
    theme_classic() + 
    geom_boxplot(key_glyph = "rect", outlier.shape = NA, linewidth = 0.3, alpha = 0.6) +
    labs(x = "", y = "Expression") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.2, "cm")
          #legend.position = "inside",
          #legend.position.inside = c(0.8, 0.2)
          ) +
    ylim(0, ylimvec[z]) +
    scale_fill_manual(values = c("blue", "red"), labels = c("CTRL", "Disease"))
leg <- cowplot::get_legend(p)
pdf(file.path(save_dir, "legend.pdf"), width = 0.6, height = 0.4)
grid.newpage()
grid.draw(leg)
dev.off()


# # Alternative (facet_wrap)
# ball <- list()
# for(z in 1:4){ # four plots
#   b1 <- list()
#   for(i in 1:3){ # three celltypes
#     conmat <- data.frame("CT" = paste0("CT", i),
#                         "Disease" = c("Con"),
#                         "Value" = runif(20, distlistcon[[z]][[i]][1], distlistcon[[z]][[i]][2]))

#     ADmat <- data.frame("CT" = paste0("CT", i),
#                         "Disease" = c("AD"),
#                         "Value" = runif(20, distlistAD[[z]][[i]][1], distlistAD[[z]][[i]][2]))

#     b1[[i]] <- rbind(conmat, ADmat)
#   }
#   b1dat <- do.call(rbind, b1)
#   b1dat$mod <- paste0("Module ", z)
#   ball[[z]] <- b1dat
# }
# balldf <- do.call(rbind, ball)
# p <- ggplot(balldf, aes(x = CT, y = Value, fill = Disease)) +
#   theme_classic() + 
#   geom_boxplot(outlier.shape = NA) +
#   labs(x = "", y = "Expression") +
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 5),
#         legend.key.size = unit(0.2, "cm"),
#         legend.position = "none"
#         #legend.position = "inside",
#         #legend.position.inside = c(0.8, 0.2)
#         ) +
#   ylim(0, 3.5) +
#   facet_wrap(~mod, nrow = 4)
# ggsave(p, file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/schematic/pall_facetwrap.pdf")), width = 4, height = 5)

######################################
# Plots for individual genes (2 total)
######################################

indval <- list(c(1, 2, 1.6, 2.4, 2.3, 3.1),
               c(1.2, 1, 2.1, 2.9, 1.4, 2.5))

for(z in 1:2){
  genedf1 <- data.frame("Gene" = c(rep("Gene 1", 2), rep("Gene 2", 2), rep("Gene 3", 2)),
                        "Type" = rep(c("CTRL", "Disease"), 3),
                        "Values" = indval[[z]])
  genedf1$Type <- factor(genedf1$Type, levels = c("CTRL", "Disease"))
  p1 <- ggplot(genedf1, aes(x = Type, y = Values, group = Gene, color = Gene)) +
    theme_classic() + 
    geom_point(key_glyph = "point") + 
    geom_line(key_glyph = "point") +
    labs(x = "", y = "Expression") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.spacing.y = unit(2, "cm"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
    scale_color_manual(values = c("#E69F00", "#009E73", "#CC79A7")) +       
    guides(color = guide_legend(byrow = TRUE)) +
    ylim(0.5, NA)
  ggsave(p1, file = file.path(save_dir, paste0("g", z, ".pdf")), width = 2.1, height = 1.5)
}
