# Gather lists of ad-associated genes and examine enrichment in modules
# /home/gugene/data_other/ad_genesets/

library(data.table)
library(tidyverse)
library(qs)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "generic_enrichment_function.r"))
options(bitmapType = 'cairo')

# Load modules (bulk)
bulk_seed <- qs::qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
bulk_bc <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
    (\(x) split(x$Gene, x$topmodposbc))()
sn_genelist <- rownames(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
all_genes <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)[,2]

# Load lists of significant modules shared between Morabito and SEAAD2024
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/all_module_index_list.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/pos_module_index_list.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/COPA_conserve/morabito_vs_SEAAD2024/modules_shared_between_all_datasets/neg_module_index_list.qs")))

# # Load ad genesets
# adlist <- list()
# adlist[[1]] <- fread(data.table=F,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "ad_genesets/advp.variant.records.hg38.tsv")) |> 
#     (\(x) split(x$`Pubmed PMID`, x$nearest_gene_symb))() |> lapply(length) |> lapply(\(x) x[!duplicated(x)]) |> unlist() |> sort(decreasing=T) |> (\(x) x[x>=10])() |> names()
# adlist[[2]] <- fread(data.table=F,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "ad_genesets/Bellenguez_gwas-association-downloaded_2025-05-23-accessionId_GCST90027158.tsv")) |> pull(MAPPED_GENE)
# adlist[[3]] <- fread(data.table=F,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "ad_genesets/wightman2021_table1.csv")) |> pull(Gene)
# adlist[[4]] <- fread(data.table=F,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "phenopedia/oct2024/brain.csv")) |> dplyr::filter(disease=="Alzheimer Disease") |>
#     dplyr::filter(n_publications>=3) |> pull(gene)
# adlist[[5]] <- c("APOE", "ACE", "CHRNB2", "CST3", "ESR1", "GAPDHS", "IDE", "MTHFR", "NCSTN", "PRNP", "PSEN1", "TF", "TFAM", "TNF") # Alzgene
# names(adlist) <- c("ADVP", "Bellenguez et al", "Wightman et al", "Phenopedia", "Alzgene")

# # Clean genenames
# adlist[[1]][adlist[[1]]=="APOC4"] <- "APOC4-APOC2"
# adlist[[1]][adlist[[1]]=="LINC02210-CRHR1"] <- "CRHR1"
# adlist[[1]] <- adlist[[1]][adlist[[1]]!=""]
# rem <- adlist[[2]][!adlist[[2]] %in% sn_genelist] 
# rem2 <- rem |>
#     strsplit(" - ") |> unlist() |> strsplit(", ") |> unlist() |> (\(x) x[x %in% sn_genelist])()
# adlist[[2]] <- adlist[[2]][!adlist[[2]] %in% rem] |> (\(x) c(x, rem2))()
# adlist[[3]][adlist[[3]]=="INPPD5"] <- "INPP5D"
# adlist[[3]] <- adlist[[3]] |> strsplit("/") |> unlist() |> (\(x) x[x %in% sn_genelist])()
# adlist[[4]][adlist[[4]]=="PVRL2"] <- "NECTIN2"
#qsave(adlist, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "ad_genesets/five_ad_genesets.qs"))
adlist <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "ad_genesets/five_ad_genesets.qs"))

# Run enrichment
sigi <- euc_dist$all[[1]] |> unlist() |> unique() # all, subclass, fdr - the most enriched mods don't show up in sig mods with fdr cutoff
#sigi <- euc_dist$neg[[1]] |> unlist() |> unique() # all, subclass, nom - but they do show up with nom cutoff
#sigi <- euc_dist$neg[[1]] |> unlist() |> unique() # all, subclass, nom - but they do show up with nom cutoff

# sig mods only
enr1 <- GSHG_custom(bulk_bc[sigi], 
                    adlist, 
                    all_genes)
# sig mods removed
enr2 <- GSHG_custom(bulk_bc[-sigi], 
                    adlist, 
                    all_genes)

plotdf <- list(enr1, enr2) |>
    lapply(\(x) data.frame("set"=x[,1], "pval"=apply(x[-1],1,min))) |>
    (\(.) mapply(\(x,y){x$type=y;return(x)}, ., c("sig", "no_sig"), SIMPLIFY=F))() |>
    do.call(what="rbind")|>
    mutate(pval=-log10(pval),
           set=factor(set, levels=rev(unique(set))),
           type=factor(type,levels=c("no_sig", "sig")))


p <- 
    ggplot(plotdf, aes(y=set, x=pval,fill=type)) +
        theme_bw() + 
        geom_bar(stat="identity", position="dodge") + 
        geom_vline(xintercept=-log10(0.05/10), color="blue",linetype="dashed", linewidth=1) + 
        theme(text=element_text(size=30),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="bottom",
            legend.title=element_blank()) +
        labs(y="", x=bquote({"Enrichment p-value (log"[10]~")"})) +
        guides(fill = guide_legend(nrow = 2)) +
        scale_fill_manual(values=c("no_sig"="grey", "sig"="red"),breaks=c("sig","no_sig"), labels=c("sig"="Significant modules","no_sig"="Other modules"))
ggsave(p,file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/ad_enrich_BCmods_allsigmodsNom_morabitoSEAshared.png"), height=6, width=10)



# Repeat but with only SEAAD2024 significant modules
euc_dist_sea <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_all.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_positive.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_negative.qs")))
sigi_sea <- euc_dist_sea$all[[1]] |> unlist() |> unique() # all, subclass, fdr - the most enriched mods don't show up in sig mods with fdr cutoff
# sig mods only
enr1sea <- GSHG_custom(bulk_bc[sigi_sea], 
                    adlist, 
                    all_genes)
# sig mods removed
enr2sea <- GSHG_custom(bulk_bc[-sigi_sea], 
                    adlist, 
                    all_genes)

plotdfsea <- list(enr1sea, enr2sea) |>
    lapply(\(x) data.frame("set"=x[,1], "pval"=apply(x[-1],1,min))) |>
    (\(.) mapply(\(x,y){x$type=y;return(x)}, ., c("sig", "no_sig"), SIMPLIFY=F))() |>
    do.call(what="rbind")|>
    mutate(pval=-log10(pval),
           set=factor(set, levels=rev(unique(set))),
           type=factor(type,levels=c("no_sig", "sig")))

psea <- 
    ggplot(plotdfsea, aes(y=set, x=pval,fill=type)) +
        theme_bw() + 
        geom_bar(stat="identity", position="dodge") + 
        geom_vline(xintercept=-log10(0.05/10), color="blue",linetype="dashed", linewidth=1) + 
        theme(text=element_text(size=30),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="bottom",
            legend.title=element_blank()) +
        labs(y="", x=bquote({"Enrichment p-value (log"[10]~")"})) +
        guides(fill = guide_legend(nrow = 2)) +
        scale_fill_manual(values=c("no_sig"="grey", "sig"="red"),breaks=c("sig","no_sig"), labels=c("sig"="Significant modules","no_sig"="Other modules"))
ggsave(psea,file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/ad_enrich_BCmods_allsigmodsNom_SEAAD2024only.png"), height=6,width=10)

