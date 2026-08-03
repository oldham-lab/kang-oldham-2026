library(ggh4x)

make_expr_line_plots <- function(expr,
                                 mods,
                                 these_mods,
                                 datkme,
                                 highlight_genes = NULL # gene names to colour red in the legend
                                 #seed = 23
                                 ){
  #set.seed(seed)
  num_of_genes <- 10 # number of genes to plot
  cols2 <- RColorBrewer::brewer.pal(10, "Spectral")

  # Calculate expression z-scores
  expr_t <- t(expr[,3:ncol(expr)])
  colnames(expr_t) <- expr[,2]
  expr_z <- apply(expr_t, 2, function(x) (x-mean(x))/sd(x))

  # Calculate # var explained
  r2_vec_list <- list()
  for(j in seq_along(mods)){
    expr_plot <- as.data.frame(expr_z[ ,colnames(expr_z) %in% mods[[j]]])
    pc1 <-  prcomp(expr_plot[apply(expr_plot,2,var)>0], scale=T)$x[,1]
    r2_vec <- apply(expr_plot[apply(expr_plot,2,var)>0],2,function(x) summary(lm(x~pc1))$r.squared)
    r2_vec_list[[j]] <- r2_vec
  }
  vevec <- lapply(r2_vec_list, mean) |> unlist()

  # Optionally render selected gene names red in the legend
  if(is.null(highlight_genes)){
    gene_labeller    <- ggplot2::waiver()
    legend_text_elem <- element_text(size = 14)
  } else {
    if(!requireNamespace("ggtext", quietly = TRUE))
      stop("highlight_genes requires the 'ggtext' package")
    gene_labeller <- function(g) ifelse(g %in% highlight_genes,
                                        paste0("<span style='color:red'>", g, "</span>"), g)
    legend_text_elem <- ggtext::element_markdown(size = 14)
  }

  expr_plots <- list()
  for(j in 1:length(mods)){
    cat("\r", j, "out of", length(these_mods))                  
    # Sort genes by kME
    mod_kme <- datkme[, c(1:3, these_mods[[j]] + 4)]
    mod_kme <- mod_kme[order(mod_kme[, 4], decreasing=T), ]
    mod_kme <- mod_kme[which(mod_kme$topmodposbc == these_mods[[j]]), ]
    
    # Select random samples to graph
    samp_ind1 <- sample(1:nrow(expr_z), 100)
    expr_plot <- as.data.frame(expr_z[samp_ind1, colnames(expr_z) %in% mods[[j]]])
    while(all(apply(expr_plot, 2, var) == 0)){
      samp_ind1 <- sample(1:nrow(expr_z), 100)
      expr_plot <- as.data.frame(expr_z[samp_ind1, colnames(expr_z) %in% mods[[j]]])
    }
    
    # Select gene expr vectors 
    expr_plot <- expr_plot[,colnames(expr_plot) %in% mod_kme[1:min(num_of_genes, length(mods[[j]])),2]]
    expr_plot$sample <- c(1:nrow(expr_plot))
    expr_plot <- pivot_longer(as.data.frame(expr_plot),cols=!sample, names_to = "gene", values_to = "value")
    expr_plot$gene <- factor(expr_plot$gene, levels=mod_kme[1:length(unique(expr_plot$gene)),2])
    
    # graph
    if(length(mods[[j]]) >= 10){
      legt <- "Top 10 genes:"
    } else {
      legt <- "Genes:"
    }
    exprsubt <- paste0("# of module genes: ", length(mods[[j]]), "\n", sprintf("%.0f", vevec[j] * 100), "% variance explained by module eigengene")
    expr_plots[[j]] <- ggplot(expr_plot, aes(x=sample, y=value, color=gene )) +
      theme_light() +
      geom_line(linewidth = 0.2, alpha=0.7) +
      labs(title = paste0("Module ", j),
           subtitle = exprsubt,
           x = "Bulk sample (n = 100 random samples)",
           y = "Expression z-score",
           colour = legt
           ) + 
      scale_color_manual(values = cols2, labels = gene_labeller) +
      theme(text = element_text(size = 16, family = "sans"),
            legend.box.margin = margin(0, 0, 0, -10),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
            axis.text.x = element_text(size = 14), 
            axis.text.y = element_text(size = 14), 
            axis.title.x = element_text(size = 14), 
            axis.title.y = element_text(size = 14), 
            legend.position = "right",
            legend.title = element_text(size = 14),
            legend.text = legend_text_elem,
            legend.key.size = unit(0.2, "cm")
            ) +
      guides(color = guide_legend(override.aes = list(linewidth = 2)))
  }
  return(expr_plots)
} # EOF

make_gsea_plots <- function(gsea,
                            gsea_broad,
                            these_mods,
                            force_width = F # force plotting area to be 1 inch wide
                            ){

  cat("Creating plotting objects (GSEA)...\n")

  # process broad and combine GSEA data
  broad_overview <- fread(system.file("extdata", "broad_geneset_overview.csv", package = "CoPA"), data.table = F) |>
    dplyr::filter(!catType %in% c("c1", "c3", "c4", "c6", "c7")) |>
    dplyr::filter(!grepl("^HP_", setNames))
  gsea_broad <- gsea_broad[gsea_broad$SetID %in% broad_overview$setIDs, ]
  gsea_comb <- rbind(gsea, gsea_broad)

  # calculate fdr cutoff
  gseavec <- unlist(list(gsea_comb[,4:ncol(gsea_comb)]))
  gsq <- qvalue::qvalue(gseavec)
  gsea_cutFDR <- max(gsq$pvalues[gsq$qvalues<.05])
  
  # -log10 transform gsea p-values
  for(k in 4:ncol(gsea_comb)){
    gsea_comb[,k] <- -log10(gsea_comb[,k])
    gsea_comb[gsea_comb[,k]==Inf,k] <- 308
  }
  gsea_cutBC <- -log10(0.05/(nrow(gsea_comb) * (ncol(gsea_comb)-3)))
  gsea_cutFDR <- -log10(gsea_cutFDR)

  # Add breaks to geneset names (for graphing)
  # insert_char_every_n <- function(string, char_to_insert = "\n", n = 50) {
  #   gsub(paste0("(.{", n, "})"), paste0("\\1", char_to_insert), string)
  # }
  insert_char_after_underscore <- function(s, char_to_insert = "\n", n = 50) {
    library(stringr)
    # Find all underscore positions
    underscore_positions <- str_locate_all(s, "_")[[1]][, "start"]

    # Find the first underscore position after n characters
    insertion_point <- NA
    for (pos in underscore_positions) {
      if (pos > n) {
        insertion_point <- pos
        break
      }
    }

    if (is.na(insertion_point)) {
      # No underscore found after n characters, return original string
      return(s)
    } else {
      # Split the string and insert the character
      part1 <- str_sub(s, 1, insertion_point)
      part2 <- str_sub(s, insertion_point + 1, str_length(s))
      return(paste0(part1, char_to_insert, part2))
    }
  }
  # new_names <- unlist(lapply(gsea_comb$SetName, insert_char_after_underscore))
  # new_names <- unlist(lapply(new_names, \(x) insert_char_after_underscore(x, n = 100)))
  # new_names <- unlist(lapply(new_names, \(x) insert_char_after_underscore(x, n = 150)))
  new_names <- unlist(lapply(gsea_comb$SetName, \(x){
    if(nchar(x) >= 100){
      return(insert_char_after_underscore(x, n = nchar(x) / 2))
    } else {
      return(x)
    }
  }))
  gsea_comb$SetName_breaks <- new_names
  gsea_comb <- gsea_comb |> relocate(SetName_breaks, .after = SetName)

  gsea_plots <- list()
 # gsea_plots_breaks <- list()
  for(j in seq_along(these_mods)){
    cat("\r", j, "out of", length(these_mods))                  
    ## enrichment featuring most significant genesets from our lab collection + sn genesets + broad
    # order gsea data by p-value
    gsea_ind <- which(colnames(gsea_comb) == paste0("X", these_mods[j]))
    if(length(gsea_ind)>0){
      gsea_plot <- gsea_comb[order(gsea_comb[,gsea_ind], decreasing=T),]
      gsea_plot <- gsea_plot[1:10, c(1:4, gsea_ind)]
      colnames(gsea_plot)[5] <- "pval"
      gsea_plot$SetName <- factor(gsea_plot$SetName,levels=unique(gsea_plot$SetName))
      gsea_plot$SetName_breaks <- factor(gsea_plot$SetName_breaks,levels=unique(gsea_plot$SetName_breaks))
      ysizevar <- (max(nchar(as.character(gsea_plot$SetName))) - 31)/58 * -4 + 10
      ysizevar <- max(ysizevar, 6)
      ysizevar <- min(ysizevar, 10)

      plotdf <- gsea_plot |> 
        dplyr::slice(1:5) |>
        dplyr::mutate(SetName = factor(SetName, levels = rev(unique(SetName))),
                      SetName_breaks = factor(SetName_breaks, levels = rev(unique(SetName_breaks)))) 
      gsea_plots[[j]] <- ggplot(plotdf, aes(x = SetName_breaks, y = pval)) +
           theme_light() + 
           theme(text = element_text(family = "sans"),
                                     axis.text.x = element_text(size = 14), 
                                     axis.text.y = element_text(size = 14), 
                                     axis.title.x = element_text(size = 14), 
                                     axis.title.y = element_text(size = 14), 
                                     plot.subtitle = element_blank(),
                                     legend.key.size = unit(0.2, "cm"), 
                                     legend.title = element_blank()) +
          geom_bar(stat="identity") + 
          coord_flip() +
          theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5, size = 10, color = "black"),
                axis.text.y = element_text(hjust = 1, vjust = 0.5, size = ysizevar, color = "black"),
                axis.title.y = element_text(size = 11, color = "black"),
                axis.title.x = element_text(size = 11, color = "black"),
                plot.title = element_blank()) + 
          labs(x = "", y = bquote(-log[10]~"(p-val)"), title = "Geneset enrichment") +
        #  scale_x_discrete(limits = levels(plotdf$SetName)) +
          geom_hline(yintercept = gsea_cutFDR, color = "red") 
      if(force_width){
        gsea_plots[[j]] <- gsea_plots[[j]] + 
          force_panelsizes(rows = unit(1, "in"), cols = unit(1, "in"))  
      }

    } else {
      gsea_plots[[j]] <- ggplot() +
        theme_void() +
        geom_text(aes(x=1,y=1,label="[NA]"), size=24)
    }
  }
  return(gsea_plots)
}