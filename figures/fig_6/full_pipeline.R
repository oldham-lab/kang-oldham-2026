fig6_full_pipeline <- function(all_plots,
                               mods,
                               these_mods,
                               sea_dir,
                               mit_dir,
                               snapshot_out_dir,
                               keys,
                               title_vec,
                               title_vec_dataset,
                               dx_list,
                               save_suffix_vec,
                               target_mods = NULL,
                               flat_output = FALSE){

    # Load indices (native, REI) and SE
    index_dirs <- c("log_native", "log_REI")
    index_save_names <- c("native_log", "REI")
    index_xaxis_names <- c("Mean expression (log UMI counts + 1)", "Relative expression index")
    #for(d in 1:2){
    d <- 1 # native_log

    #### 0. Create various vectors and lists for all comparisons

    # Helper functions for loading
    make_name <- function(path, base, stem_prefix) {
    rel   <- sub(paste0(base, "/"), "", path)
    parts <- strsplit(rel, "/")[[1]]
    region <- parts[1]                                    # PFC or MTC
    stem   <- tools::file_path_sans_ext(basename(path))   # e.g. genomewide_means_allAD
    stem   <- sub(paste0("^", stem_prefix, "_"), "", stem) # strip leading prefix
    paste(region, stem, sep = "_")
    }

    make_se_name <- function(path, base) {
    rel    <- sub(paste0(base, "/"), "", path)
    parts  <- strsplit(rel, "/")[[1]]
    region    <- parts[1]                                   # PFC or MTC
    subfolder <- parts[3]                                   # log_native / log_normByMean / log_REI
    stem      <- tools::file_path_sans_ext(basename(path)) # se_{group}_{modtype}
    stem      <- sub("^se_", "", stem)                      # {group}_{modtype}
    paste(region, stem, subfolder, sep = "_")
    }

    make_dcopa_name <- function(path, base) {
    rel    <- sub(paste0(base, "/"), "", path)
    parts  <- strsplit(rel, "/")[[1]]
    region <- parts[1]                                     # PFC or MTC
    stem   <- tools::file_path_sans_ext(basename(path))    # e.g. allAD_vs_Con_bulk_megaset_output_table
    stem   <- sub("_output_table$", "", stem)              # allAD_vs_Con_bulk_megaset
    paste(region, stem, sep = "_")
    }

    ### Load all SEA data
    sea_files <- list.files(sea_dir, 
                            pattern = "\\.csv$",
                            recursive = TRUE, full.names = TRUE)
    # SEA genomewide means
    sea_means_files <- sea_files[grepl("/means/", sea_files)]
    sea_means_list <- setNames(
    lapply(sea_means_files, \(x) fread(x, data.table = F)),
    sapply(sea_means_files, make_name, base = sea_dir, stem_prefix = "genomewide_means")
    )

    # SEA log native projections.
    # Exclude "_seed" files (seed-gene-only projections, added for a separate analysis):
    # Fig 6 uses the full-module topmodposbc projection, and a "_seed" sibling would make
    # the per-key grepl() lookups below match two files (-> recursive [[ ]] -> crash).
    sea_mod_means_files <- sea_files[grepl("/mod_means/log_native/", sea_files) & !grepl("_seed", sea_files)]
    sea_mod_means_list <- setNames(
    lapply(sea_mod_means_files, \(x) fread(x, data.table = F)),
    sapply(sea_mod_means_files, make_name, base = sea_dir, stem_prefix = "mod_means")
    )

    # SEA log native projections (SE)
    sea_se_files <- sea_files[grepl("/se/log_native/", sea_files) & !grepl("_seed", sea_files)]
    sea_se_list <- setNames(
    lapply(sea_se_files, \(x) fread(x, data.table = F)),
    sapply(sea_se_files, make_se_name, base = sea_dir)
    )

    # SEA dcopa output
    sea_dcopa_files <- sea_files[grepl("/euclidean_distances/", sea_files)]
    sea_dcopa_list <- setNames(
    lapply(sea_dcopa_files, \(x) fread(x, data.table = F)),
    sapply(sea_dcopa_files, make_dcopa_name, base = sea_dir)
    )

    ### Load all MIT data
    mit_files <- list.files(mit_dir, 
                            pattern = "\\.csv$",
                            recursive = TRUE, full.names = TRUE)

    # MIT genomewide means
    mit_means_files <- mit_files[grepl("/means/", mit_files)]
    mit_means_list <- setNames(
    lapply(mit_means_files, \(x) fread(x, data.table = F)),
    sapply(mit_means_files, make_name, base = mit_dir, stem_prefix = "genomewide_means")
    )

    # MIT log native projections (exclude "_seed" siblings, see SEA note above)
    mit_mod_means_files <- mit_files[grepl("/mod_means/log_native/", mit_files) & !grepl("_seed", mit_files)]
    mit_mod_means_list <- setNames(
    lapply(mit_mod_means_files, \(x) fread(x, data.table = F)),
    sapply(mit_mod_means_files, make_name, base = mit_dir, stem_prefix = "mod_means")
    )

    # MIT log native projections (SE)
    mit_se_files <- mit_files[grepl("/se/log_native/", mit_files) & !grepl("_seed", mit_files)]
    mit_se_list <- setNames(
    lapply(mit_se_files, \(x) fread(x, data.table = F)),
    sapply(mit_se_files, make_se_name, base = mit_dir)
    )

    # MIT dcopa output
    mit_dcopa_files <- mit_files[grepl("/euclidean_distances/", mit_files)]
    mit_dcopa_list <- setNames(
    lapply(mit_dcopa_files, \(x) fread(x, data.table = F)),
    sapply(mit_dcopa_files, make_dcopa_name, base = mit_dir)
    )

    ##### Run loop for all comparisons
    for(i in 1:length(keys)){
    #for(i in 1:8){ # Need to create separate objects for ROSMAP (new plot objects, new input set of mods)

    # Load genomewide means by ct
    mean_sum_gab <- list(
        sea_means_list[[which(grepl(keys[[i]][1], names(sea_means_list)) & grepl(keys[[i]][2], names(sea_means_list)))]],
        sea_means_list[[which(grepl(keys[[i]][1], names(sea_means_list)) & grepl(keys[[i]][3], names(sea_means_list)))]]
    ) |> 
        lapply(\(x) x |> tibble::column_to_rownames(names(x)[1]) |> as.data.frame()) |>
        #lapply(\(x) x |> tibble::column_to_rownames("gene_ids") |> as.data.frame()) |>
        setNames(dx_list[[i]])

    mean_sum_mit <- list(
        mit_means_list[[which(grepl(keys[[i]][1], names(mit_means_list)) & grepl(keys[[i]][2], names(mit_means_list)))]],
        mit_means_list[[which(grepl(keys[[i]][1], names(mit_means_list)) & grepl(keys[[i]][3], names(mit_means_list)))]]
    ) |> 
        lapply(\(x) x |> tibble::column_to_rownames(names(x)[1]) |> as.data.frame()) |>
        #lapply(\(x) x |> tibble::column_to_rownames("gene_sym") |> as.data.frame()) |>
        setNames(dx_list[[i]])

    # Load dcopa data
    dcopa_gab <- sea_dcopa_list[[which(grepl(keys[[i]][1], names(sea_dcopa_list)) & grepl(keys[[i]][2], names(sea_dcopa_list)) & grepl(keys[[i]][3], names(sea_dcopa_list)) & grepl(keys[[i]][4], names(sea_dcopa_list)))]]
    dcopa_mit <- mit_dcopa_list[[which(grepl(keys[[i]][1], names(mit_dcopa_list)) & grepl(keys[[i]][2], names(mit_dcopa_list)) & grepl(keys[[i]][3], names(mit_dcopa_list)) & grepl(keys[[i]][4], names(mit_dcopa_list)))]]
    dcopa_shared <- dplyr::inner_join(dcopa_gab, dcopa_mit, by = dplyr::join_by(mod, Direction, Consistency), suffix = c("1", "2"), relationship = "many-to-many") |>
        dplyr::filter(Celltype1 == Celltype2,
                        sig_FDR1,
                        sig_FDR2,
                        Consistency %in% c(0, 1)) |>
        dplyr::arrange(mod)

    # Restrict the per-module loop to user-specified modules (if provided).
    # Default (NULL) keeps original behaviour: all dCoPA-shared modules.
    loop_mods <- unique(dcopa_shared$mod)
    if(!is.null(target_mods)){
        missing_mods <- setdiff(target_mods, loop_mods)
        if(length(missing_mods) > 0)
            warning("Requested module(s) not dCoPA-significant in both datasets for key ",
                    i, " (skipped): ", paste(missing_mods, collapse = ", "))
        loop_mods <- loop_mods[loop_mods %in% target_mods]
    }

    ## Load indices
    # Gabitto AD indices
    gabad <- sea_mod_means_list[[which(grepl(keys[[i]][1], names(sea_mod_means_list)) & grepl(keys[[i]][2], names(sea_mod_means_list)) & grepl(keys[[i]][4], names(sea_mod_means_list)))]]
    gabad_se <- sea_se_list[[which(grepl(keys[[i]][1], names(sea_se_list)) & grepl(keys[[i]][2], names(sea_se_list)) & grepl(keys[[i]][4], names(sea_se_list)))]]
    
    # Gabitto control indices
    gabcon <- sea_mod_means_list[[which(grepl(keys[[i]][1], names(sea_mod_means_list)) & grepl(keys[[i]][3], names(sea_mod_means_list)) & grepl(keys[[i]][4], names(sea_mod_means_list)))]]
    gabcon_se <- sea_se_list[[which(grepl(keys[[i]][1], names(sea_se_list)) & grepl(keys[[i]][3], names(sea_se_list)) & grepl(keys[[i]][4], names(sea_se_list)))]]
    
    # Liu AD indices
    mitad <- mit_mod_means_list[[which(grepl(keys[[i]][1], names(mit_mod_means_list)) & grepl(keys[[i]][2], names(mit_mod_means_list)) & grepl(keys[[i]][4], names(mit_mod_means_list)))]]
    mitad_se <- mit_se_list[[which(grepl(keys[[i]][1], names(mit_se_list)) & grepl(keys[[i]][2], names(mit_se_list)) & grepl(keys[[i]][4], names(mit_se_list)))]]
    
    # Liu control indices
    mitcon <- mit_mod_means_list[[which(grepl(keys[[i]][1], names(mit_mod_means_list)) & grepl(keys[[i]][3], names(mit_mod_means_list)) & grepl(keys[[i]][4], names(mit_mod_means_list)))]]
    mitcon_se <- mit_se_list[[which(grepl(keys[[i]][1], names(mit_se_list)) & grepl(keys[[i]][3], names(mit_se_list)) & grepl(keys[[i]][4], names(mit_se_list)))]]
    
    # Find intersection of all cts over all datasets
    allcts <- intersect(intersect(colnames(gabcon), colnames(gabad)), intersect(colnames(mitcon), colnames(mitad)))

    # Order cell types to match the Fig. 4 barplots: excitatory neurons by layer,
    # then inhibitory interneurons, then non-neuronal (overriding the data's default
    # alphabetical column order). Both AD (Astrocyte/Endothelial/...) and SCZ/brainSCOPE
    # (Astro/Endo/...) namings are listed; any cell type not in this reference (e.g.
    # brainSCOPE Immune/PC/SMC) is appended at the end in its original order so none drop.
    fig4_ct_order <- c(
        "L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
        "Chandelier", "Pvalb", "Sst", "Sst Chodl", "Lamp5 Lhx6", "Lamp5", "Pax6", "Sncg", "Vip",
        "Astrocyte", "Astro", "Endothelial", "Endo", "Microglia-PVM", "Micro",
        "Oligodendrocyte", "Oligo", "OPC", "VLMC"
    )
    allcts <- c(intersect(fig4_ct_order, allcts), setdiff(allcts, fig4_ct_order))

    # Align columns
    gabad     <- gabad     |> dplyr::select(all_of(allcts))
    gabad_se  <- gabad_se  |> dplyr::select(all_of(allcts))
    gabcon    <- gabcon    |> dplyr::select(all_of(allcts))
    gabcon_se <- gabcon_se |> dplyr::select(all_of(allcts))
    mitad     <- mitad     |> dplyr::select(all_of(allcts))
    mitad_se  <- mitad_se  |> dplyr::select(all_of(allcts))
    mitcon    <- mitcon    |> dplyr::select(all_of(allcts))
    mitcon_se <- mitcon_se |> dplyr::select(all_of(allcts))

    # Capitalize celltypes
    CT_RENAME <- c(
        "Lamp5"     = "LAMP5",
        "Lamp5 Lhx6"= "LAMP5 LHX6",
        "Pax6"      = "PAX6",
        "Pvalb"     = "PVALB",
        "Sncg"      = "SNCG",
        "Sst"       = "SST",
        "Vip"       = "VIP",
        "Sst Chodl" = "SST CHODL"
    )

    allcts_cap <- allcts |> dplyr::recode(!!!CT_RENAME)

    # Capitalize index titles
    colnames(gabad)     <- allcts_cap
    colnames(gabad_se)  <- allcts_cap
    colnames(gabcon)    <- allcts_cap
    colnames(gabcon_se) <- allcts_cap
    colnames(mitad)     <- allcts_cap
    colnames(mitad_se)  <- allcts_cap
    colnames(mitcon)    <- allcts_cap
    colnames(mitcon_se) <- allcts_cap

    # Capitalize dcopa_shared celltype columns
    dcopa_shared <- dcopa_shared |> 
        mutate(Celltype1 = dplyr::recode(Celltype1, !!!CT_RENAME),
            Celltype2 = dplyr::recode(Celltype2, !!!CT_RENAME))
    dcopa_gab <- dcopa_gab |> 
        mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))
    dcopa_mit <- dcopa_mit |> 
        mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

    # Capitalize mean_sum_* column names
    mean_sum_gab <- lapply(mean_sum_gab, \(x){
        colnames(x) <- dplyr::recode(colnames(x), !!!CT_RENAME)
        return(x)
    })

    mean_sum_mit <- lapply(mean_sum_mit, \(x){
        colnames(x) <- dplyr::recode(colnames(x), !!!CT_RENAME)
        return(x)
    })

    # Plot projection indices
    projplotlist1 <- lapply(loop_mods, \(j){ # Gabitto

        # Fetch projection index data for the module
        plotdf1 <- data.frame("ct" = rep(allcts_cap, 2),
                            "dataset" = c(rep(paste0(title_vec_dataset[1], ", ", title_vec[i]), 2 * length(allcts))),
                            "ind" = c(unlist(gabcon[j, ]), 
                                        unlist(gabad[j, ])),
                            "ind_se" = c(unlist(gabcon_se[j, ]), 
                                        unlist(gabad_se[j, ])),
                            "dx" = c(rep(dx_list[[i]][2], length(allcts)),
                                    rep(dx_list[[i]][1], length(allcts)))) |>
        dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                        ct = factor(ct, levels = allcts_cap))

        # Fetch dCoPA significance data for given module
        highlight_mat1 <- dcopa_shared |> 
        filter(mod == j, sig_FDR1 == T) |>
        dplyr::select(Celltype1, P_value1) |>
        filter(!duplicated(Celltype1))
        
        # Plot grouped barplot of projection indices with asterisks for dCoPA significance
        plotdf1 <- plotdf1 |>
        left_join(highlight_mat1, by = join_by(ct == Celltype1)) |>
        mutate(
            label = case_when(
            P_value1 > 0.05 ~ "",
            (0.01 < P_value1 & P_value1 <= 0.05) ~ "*",
            (0.001 < P_value1 & P_value1 <= 0.01) ~ "**",
            (0.0001 < P_value1 & P_value1 <= 0.001) ~ "***",
            P_value1 <= 0.0001 ~ "****",
            .default = NA
            )) |>
        dplyr::select(!P_value1)

        plotdf2 <- data.frame("ct" = rep(allcts_cap, 2),
                            "dataset" = c(rep(paste0(title_vec_dataset[2], ", ", title_vec[i]), 2 * length(allcts))),
                            "ind" = c(unlist(mitcon[j, ]),
                                        unlist(mitad[j, ])),
                            "ind_se" = c(unlist(mitcon_se[j, ]),
                                        unlist(mitad_se[j, ])),
                            "dx" = c(rep(dx_list[[i]][2], length(allcts)),
                                    rep(dx_list[[i]][1], length(allcts)))) |>
        dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                        ct = factor(ct, levels = allcts_cap))
        
        # Fetch dCoPA significance data for given module
        highlight_mat2 <- dcopa_shared |> 
        filter(mod == j, sig_FDR2 == T) |>
        dplyr::select(Celltype2, P_value2) |>
        filter(!duplicated(Celltype2))

        plotdf2 <- plotdf2 |>
            left_join(highlight_mat2, by = join_by(ct == Celltype2)) |>
            mutate(
            label = case_when(
                P_value2 > 0.05 ~ "",
                (0.01 < P_value2 & P_value2 <= 0.05) ~ "*",
                (0.001 < P_value2 & P_value2 <= 0.01) ~ "**",
                (0.0001 < P_value2 & P_value2 <= 0.001) ~ "***",
                P_value2 <= 0.0001 ~ "****",
                .default = NA
            )) |>
        dplyr::select(!P_value2)

        plotdf <- rbind(plotdf1, plotdf2) |>
        mutate(dx = factor(dx, levels = rev(dx_list[[i]])),
               # the left_joins above coerce ct from factor to character, dropping the
               # allcts_cap level order; re-apply it so the x-axis follows the Fig. 4
               # cell-type order (and ct_levels/axis_colors below align to it too).
               ct = factor(ct, levels = allcts_cap))

        # Remove duplicate labels (assign label to whichever dx has highest value)
        plotdf_label <- plotdf |> 
            filter(!is.na(label)) |>
            group_by(ct, dataset) |>
            slice_max(ind, n = 1)

        dodge_width <- 0.6
        plot_max_y <- plotdf$ind[which.max(plotdf$ind)] + plotdf$ind_se[which.max(plotdf$ind)] * 2
        ct_levels <- sort(unique(plotdf$ct))            
        arrow_cts <- unique(plotdf_label$ct)
        axis_colors <- ifelse(ct_levels %in% arrow_cts, "red", "black") 
        
        p <- ggplot(plotdf, aes(x = ct, y = ind, fill = dx)) +
        theme_classic() +
        geom_col(width = 0.6, position = position_dodge(width = dodge_width), alpha = 0.5) +
        # geom_text(data = plotdf_label, 
        #           aes(x = ct, y = ind, label = label), 
        #           vjust = -0.5#,
        #           #y = plot_max_y + 0.1
        #           ) + 
        geom_errorbar(aes(ymin = ind - 2 * ind_se,
                            ymax = ind + 2 * ind_se),
                            width = 0.2,
                            linewidth = 0.3,
                            position = position_dodge(width = dodge_width)) +
        scale_fill_manual(values = c("#90D5FF", "#FFA500")) +
        theme(text = element_text(family = "sans", color = "black", size = 12),
                legend.position = "bottom", 
                legend.title = element_blank(),
                legend.margin = margin(-0.8, 0, 0, 0, "cm"),
                axis.title.y = element_text(size = 10, margin = margin(0, 5, 0, 0)),
                axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, color = axis_colors),
                strip.text = element_text(color = "black"),
                strip.background = element_rect(fill = "white"),
                plot.margin = margin(5, 0, 0, 0)) +
        labs(y = index_xaxis_names[d], x = "") +
        scale_y_continuous(#breaks = c(0, 0.5, 1)#, 
                            limits = c(0, plot_max_y + 0.27)
                            ) +
        facet_wrap(~dataset, nrow = 2, ncol = 1, scales = "free_y") +
        geom_text(
            data = plotdf_label,
            aes(x = ct, y = ind + 2 * ind_se + 0.15, label =     
        "↓"),                                                  
            position = position_dodge(width = dodge_width),      
            size = 3.5,                                          
            vjust = 0,
            color = "red"                                            
        ) #+
        # geom_segment(                                      
        #   data = plotdf_label,
        #   aes(x = ct, xend = ct,                               
        #       y = ind + 2 * ind_se + 0.06,   # start above error bar                                              
        #       yend = ind + 2 * ind_se + 0.02), # end closer to bar                                                    
        #   position = position_dodge(width = dodge_width),
        #   arrow = arrow(length = unit(0.15, "cm"), type = "closed"),                                             
        #   linewidth = 0.4,                                     
        #   color = "black"                                      
        # )           
        return(p)
    })

    ## Plot boxplots comparing individual module gene expression changes
    boxplotlist_gab <- lapply(loop_mods, \(j){
        #cat(j, " ")
        highlight_mat <- dcopa_shared |> 
        filter(mod == j, sig_FDR2 == T) |>
        dplyr::select(Celltype1, P_value1)
        summary_mat <- mean_sum_gab

        if(length(unique(highlight_mat[,1]) > 2)){
        stripsize <- 8
        } else {
        stripsize <- 10
        }
    
        out_df <- mapply(\(x, y){
        df <- x[rownames(x) %in% mods[[which(these_mods == j)]], colnames(x) %in% highlight_mat[,1], drop = F] 
        out <- data.frame("type" = y, "gene" = rownames(df), df, check.names = F) |>
            pivot_longer(!type:gene, names_to = "ct", values_to = "vals")
        return(out)
        }, summary_mat[1:2], names(summary_mat)[1:2], SIMPLIFY = F) |>
        do.call(what = "rbind") |>
        mutate(type = factor(type, levels = rev(dx_list[[i]]))) |>
        left_join(highlight_mat, by = join_by(ct == Celltype1)) |>
        mutate(
            label = case_when(
            P_value1 > 0.05 ~ "",
            (0.01 < P_value1 & P_value1 <= 0.05) ~ "*",
            (0.001 < P_value1 & P_value1 <= 0.01) ~ "**",
            (0.0001 < P_value1 & P_value1 <= 0.001) ~ "***",
            P_value1 <= 0.0001 ~ "****",
            .default = NA
            )) |>
        dplyr::select(!P_value1)

        outdf_label <- out_df |> 
        filter(!is.na(label)) |>
        group_by(ct) |>
        slice_max(vals, n = 1)

        outdf_label <- out_df |>                               
        group_by(ct) |>                                      
        summarise(                                           
            vals = max(vals, na.rm = TRUE),
            label = unique(label[!is.na(label) & label != ""])[1]                                                
        ) |>
        filter(!is.na(label))  
        
        p <- ggplot(out_df, aes(x = type, y = vals)) +
                theme_classic() +
                geom_point(color = "lightgrey", alpha = 0.7, size = 0.4) +
                geom_line(aes(group = gene), color = "grey40", alpha = 0.7, linewidth = 0.2) +
                geom_boxplot(aes(color = type), linewidth = 0.4, notch = F, outlier.size = 0.4, alpha = 0.7, fill = NA) +
                geom_text(
                data = outdf_label,
                aes(y = vals + 0.1, label = label),
                x = 1.5,
                size = 3,
                vjust = 0
                ) +
                facet_wrap(~ct, nrow = 1) +
                labs(x = "", y = "Mean expression\n(log UMI counts + 1)", title = title_vec_dataset[1]) +
                scale_color_manual(values = c("#90D5FF", "#FFA500")) +
                scale_y_continuous(limits = c(NA, max(out_df$vals, na.rm = TRUE) + 0.4)) +
                theme(legend.position = "none",
                    plot.title = element_text(size = 12, hjust = 0.5),
                    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 1),
                    axis.title.x = element_blank(),
                    axis.title.y = element_text(size = 12),
                    axis.text.y = element_text(size = 12),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.text = element_text(size = stripsize)
                    )
        return(p)
    })

    boxplotlist_liu <- lapply(loop_mods, \(j){
        highlight_mat <- dcopa_shared |> 
        filter(mod == j, sig_FDR2 == T) |>
        dplyr::select(Celltype2, P_value2)
        summary_mat <- mean_sum_mit

        if(length(unique(highlight_mat[,1]) > 2)){
        stripsize <- 8
        } else {
        stripsize <- 10
        }

        out_df <- mapply(\(x, y){
        df <- x[rownames(x) %in% mods[[which(these_mods == j)]], colnames(x) %in% highlight_mat[,1], drop = F] 
        out <- data.frame("type" = y, "gene" = rownames(df), df, check.names = F) |>
            pivot_longer(!type:gene, names_to = "ct", values_to = "vals")
        return(out)
        }, summary_mat[1:2], names(summary_mat)[1:2], SIMPLIFY = F) |>
        do.call(what = "rbind") |>
        mutate(type = factor(type, levels = rev(dx_list[[i]]))) |>
        left_join(highlight_mat, by = join_by(ct == Celltype2)) |>
        mutate(
            label = case_when(
            P_value2 > 0.05 ~ "",
            (0.01 < P_value2 & P_value2 <= 0.05) ~ "*",
            (0.001 < P_value2 & P_value2 <= 0.01) ~ "**",
            (0.0001 < P_value2 & P_value2 <= 0.001) ~ "***",
            P_value2 <= 0.0001 ~ "****",
            .default = NA
            )) |>
        dplyr::select(!P_value2)

        outdf_label <- out_df |> 
        filter(!is.na(label)) |>
        group_by(ct) |>
        slice_max(vals, n = 1)

        outdf_label <- out_df |>                               
        group_by(ct) |>                                      
        summarise(                                           
            vals = max(vals, na.rm = TRUE),
            label = unique(label[!is.na(label) & label != ""])[1]                                                
        ) |>
        filter(!is.na(label))  
            
        p <- ggplot(out_df, aes(x = type, y = vals)) +
        theme_classic() +
        geom_point(color = "lightgrey", alpha = 0.7, size = 0.4) +
        geom_line(aes(group = gene), color = "grey40", alpha = 0.7, linewidth = 0.2) +
        geom_boxplot(aes(color = type), linewidth = 0.4, notch = F, outlier.size = 0.4, alpha = 0.7, fill = NA) +
        geom_text(
            data = outdf_label,
            aes(y = vals + 0.1, label = label),
            x = 1.5,
            size = 3,
            vjust = 0
        ) +
        facet_wrap(~ct, nrow = 1) +
        labs(x = "", y = "Mean expression\n(log UMI counts + 1)", title = title_vec_dataset[2]) +
        scale_color_manual(values = c("#90D5FF", "#FFA500")) +
        scale_y_continuous(limits = c(NA, max(out_df$vals, na.rm = TRUE) + 0.4)) +
        theme(legend.position = "none",
                plot.title = element_text(size = 12, hjust = 0.5),
                axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 1),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.text = element_text(size = stripsize)
                )

        return(p)
    })

    # Construct and save module snapshot

    # flat_output: write <mod>_<suffix>.svg straight into snapshot_out_dir.
    # default: nest under a per-key subfolder as <mod>.svg.
    if(flat_output){
        out_dir <- snapshot_out_dir
        out_file <- function(m) paste0(m, "_", save_suffix_vec[i], ".svg")
    } else {
        out_dir <- file.path(snapshot_out_dir, save_suffix_vec[i])
        out_file <- function(m) paste0(m, ".svg")
    }
    if(!dir.exists(out_dir))
        dir.create(out_dir, recursive = T)

    for(j in seq_along(loop_mods)){
        boxplot_combined <- plot_grid(boxplotlist_gab[[j]], boxplotlist_liu[[j]], nrow = 1)
        pall <- suppressMessages(plot_grid(all_plots[[1]][[which(these_mods == loop_mods[j])]],
                                           all_plots[[2]][[which(these_mods == loop_mods[j])]],
                                           projplotlist1[[j]],
                                           NULL,
                                           boxplot_combined,
                                           nrow = 5, 
                                           rel_heights = c(0.5, 0.3, 0.8, 0.05, 0.6)))
        ggsave(pall, 
            file = file.path(out_dir, out_file(loop_mods[j])),
            device = svglite::svglite, 
            width = 6, 
            height = 14, 
            bg = "white")
        cat(j, " ")
    }
    cat(i, "\n")
    }


}