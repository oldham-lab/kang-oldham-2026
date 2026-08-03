library(gt)
library(dplyr)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v8/gsea_func_optimized.R"))

create_summary_table <- function(g_overlaps,
                                 ad_db,
                                 column_4_header = "AD-associated genes",
                                 out_string = "panel_G_gsea_summary_table_mtg",
                                 save_dir){

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
    ##########
  
    # Run GSEA
    b_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                         broad = T)
    g_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                        broad = F)
    gsea_out <- rbind(b_out, g_out)

    # Calculate FDR
    gsea_pval_to_fdr <- function(gsea_out) {                              
        meta_cols <- c("SetID", "SetName", "SetSize")                       
        pval_cols <- setdiff(colnames(gsea_out), meta_cols)                 
        fdr_out <- gsea_out                                                 
        fdr_out[, pval_cols] <- lapply(gsea_out[, pval_cols, drop = FALSE], 
                                        p.adjust, method = "BH")            
        fdr_out                                                             
    }                                                                     
    gsea_out_fdr <- gsea_pval_to_fdr(gsea_out)

    # Significant genesets per celltype (p < GSEA_PTHRESH)
    # GSEA_PTHRESH <- 0.05
    # gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
    #   pvals <- gsea_out[[3 + i]]
    #   gsea_out$SetName[!is.na(pvals) & pvals < GSEA_PTHRESH]
    # }), celltypes)

    # ── AD gene database: genes in g_overlaps with documented AD associations ────
    # Each entry: gene = list(n_refs = approx. no. relevant papers, ref = bibliography no.)
    # Ordered alphabetically within each celltype when displayed.
    # Bibliography: gsea_summary_bibliography_mtg.md

    # ad_db <- list(
    # # PI3K–Akt signalling
    # AKT3     = list(n=7,  ref=1),
    # PIK3CA   = list(n=5,  ref=1),
    # # Endosomal APP sorting / Rab5
    # APPL1    = list(n=5,  ref=2),
    # # APP cytoplasmic tail / secretory trafficking
    # APPBP2   = list(n=5,  ref=3),
    # # Mitophagy receptor (PINK1/Parkin-independent)
    # BNIP3    = list(n=4,  ref=4),
    # # Autophagy / ALS-FTD endosomal pathway
    # C9orf72  = list(n=4,  ref=5),
    # # Kinesin-1 adaptor / APP axonal transport
    # CLSTN1   = list(n=5,  ref=6),
    # # Dystrophic neurites / APP vulnerability
    # CLSTN3   = list(n=3,  ref=7),
    # # Endocytosis / clathrin
    # CLTC     = list(n=7,  ref=8),
    # # COPI / APP retrograde trafficking
    # COPA     = list(n=3,  ref=9),
    # # Ubiquitin-proteasome system (SCF-E3 ligase)
    # CUL1     = list(n=3,  ref=10),
    # KLHL7    = list(n=2,  ref=10),
    # # Clathrin uncoating / endocytic recycling
    # DNAJC6   = list(n=3,  ref=11),
    # # Mitochondrial fission and fusion
    # DNM1L    = list(n=8,  ref=12),
    # OPA1     = list(n=5,  ref=12),
    # # Retrograde axonal transport / dynein
    # DYNC1H1  = list(n=5,  ref=13),
    # DYNC1I1  = list(n=5,  ref=13),
    # PAFAH1B1 = list(n=6,  ref=13),
    # # ER stress / ERAD
    # EDEM3    = list(n=3,  ref=14),
    # # Stress granules / tau condensates
    # G3BP2    = list(n=3,  ref=15),
    # # Autophagy / LC3-II
    # MAP1LC3B = list(n=8,  ref=16),
    # # DENN/MADD / neuroprotection
    # MADD     = list(n=4,  ref=17),
    # # JNK / tau hyperphosphorylation
    # MAPK8    = list(n=4,  ref=18),
    # # Neprilysin / Abeta degradation
    # MME      = list(n=5,  ref=19),
    # # mTOR / autophagy / tau clearance
    # MTOR     = list(n=5,  ref=20),
    # # Neurexin / synaptic organisation
    # NRXN1    = list(n=4,  ref=21),
    # # Oxidative stress / neuroprotection
    # OXR1     = list(n=3,  ref=22),
    # # Endolysosomal PI(3,5)P2 synthesis
    # PIKFYVE  = list(n=4,  ref=23),
    # # Tau phosphatase (PP2A)
    # PPP2R2B  = list(n=8,  ref=24),
    # # Tau phosphatase (calcineurin)
    # PPP3CB   = list(n=5,  ref=25),
    # PPP3R1   = list(n=3,  ref=25),
    # # BACE1-APP scaffolding
    # RANBP9   = list(n=5,  ref=26),
    # # Selective vulnerability / NFT-prone neurons
    # RORB     = list(n=4,  ref=27),
    # # ER morphology / BACE1 restriction
    # RTN1     = list(n=3,  ref=28),
    # RTN3     = list(n=5,  ref=28),
    # # PI(4)P lipid homeostasis / APP trafficking
    # SACM1L   = list(n=2,  ref=29),
    # # COPII / APP ER-to-Golgi trafficking
    # SEC23A   = list(n=2,  ref=30),
    # # Endophilin A1 / synaptic Aβ injury
    # SH3GL2   = list(n=3,  ref=31),
    # # Synaptic vesicle exocytosis / SNARE
    # SNAP25   = list(n=5,  ref=32),
    # # Cytoskeletal integrity / calpain cleavage
    # SPTAN1   = list(n=3,  ref=33),
    # # Store-operated Ca2+ entry / spine stability
    # STIM2    = list(n=4,  ref=34),
    # # Synaptic phosphoinositide turnover
    # SYNJ1    = list(n=9,  ref=35),
    # # Synaptic vesicle density / AD PET biomarker
    # SV2A     = list(n=5,  ref=36),
    # # Tau kinase (TAOK / MARK cascade)
    # TAOK1    = list(n=5,  ref=37),
    # # Neuroinflammation / NF-kB / TLR
    # TBK1     = list(n=4,  ref=38),
    # # Lysosomal / TDP-43 pathology
    # TMEM106B = list(n=9,  ref=39),
    # # Selective autophagy / TLR adaptor
    # TOLLIP   = list(n=5,  ref=40),
    # # Tau deubiquitination / aggregation
    # USP10    = list(n=3,  ref=41),
    # # Proteasomal tau clearance
    # USP14    = list(n=5,  ref=42),
    # # ER-mitochondria contact sites (MAM)
    # VAPB     = list(n=4,  ref=43),
    # # Mitochondria / Abeta interaction
    # VDAC1    = list(n=5,  ref=44),
    # VDAC3    = list(n=4,  ref=45),
    # # Retromer / BACE1 endosomal trafficking
    # VPS35    = list(n=5,  ref=46),
    # # Ca2+ sensor / AD CSF biomarker
    # VSNL1    = list(n=6,  ref=47),
    # # 14-3-3 / phospho-tau interaction
    # YWHAB    = list(n=5,  ref=48),
    # YWHAZ    = list(n=5,  ref=49)
    # )

    # ── Build table ─────────────────────────────────────────────────
    celltypes  <- unique(g_overlaps$Celltype)

    # Significant genesets per celltype (BH-adjusted p < GSEA_FDR_THRESH)
    # GSEA_FDR_THRESH <- 0.05
    # gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
    #   pvals <- gsea_out[[3 + i]]
    #   fdr   <- p.adjust(pvals, method = "BH")
    #   gsea_out$SetName[!is.na(fdr) & fdr < GSEA_FDR_THRESH]
    # }), celltypes)
    # qsave(gsea_sig, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6/gsea_sig_mtg.qs"))

    build_row <- function(ct) {
        genes  <- g_overlaps$Gene[g_overlaps$Celltype == ct]
        n_genes <- length(genes)

        # Top 3 GO gene sets for this celltype (columns assumed in same order as celltypes)
        ct_idx  <- which(celltypes == ct)
        pvals   <- gsea_out[[3 + ct_idx]]
        go_mask <- grepl("^GO", gsea_out$SetName)
        if (!is.null(pvals)) {
            go_pvals <- ifelse(go_mask, pvals, NA)
            sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
            n_sig    <- sum(sig_mask)
            if (n_sig == 0) {
            top_str <- "—"
            } else {
            top_idx  <- order(go_pvals)[1:min(3, n_sig)]
            top_sets <- gsea_out$SetName[top_idx]
            top_pval <- pvals[top_idx]
            top_str  <- paste0(
                seq_along(top_sets), ". ", top_sets,
                " (p=", formatC(top_pval, format = "e", digits = 2), ")",
                collapse = "<br>"
            )
            }
        }

        # AD genes: match uppercase gene names against ad_db, sort alphabetically
        genes_upper <- toupper(genes)
        hits        <- genes[genes_upper %in% names(ad_db)]
        hits_up     <- genes_upper[genes_upper %in% names(ad_db)]

        if (length(hits) > 0) {
            ref_nums <- sapply(hits_up, function(g) ad_db[[g]]$ref)
            ord      <- order(hits_up)
            hits_s   <- hits[ord]
            refs_s   <- ref_nums[ord]
            ad_str   <- paste0(hits_s, "<sup>", refs_s, "</sup>", collapse = ", ")
        } else {
            ad_str <- "—"
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
      arrange(Celltype) |>
      dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

    # ── Format with gt ─────────────────────────────────────────────────────────────
    GT_SCALE <- 0.60  # uniform scale factor for PDF table size (1 = original)
    # cols_width evaluates formula RHS via rlang tidy eval against the global env,
    # so local variables (including GT_SCALE) are invisible. rlang::inject() with !!
    # substitutes the values before gt ever sees the formulas.

    gt_table <- summary_df |>
    gt() |>
    cols_label(
        Celltype     = "Cell type",
        N_genes      = "Genes (n)",
        Top_genesets = "Top 3 GO gene sets (by p-value)",
        AD_genes     = column_4_header
    ) |>
    fmt_markdown(columns = AD_genes) |>
    fmt(columns = Top_genesets, fns = function(x) x) |>
    (\(tbl) rlang::inject(cols_width(tbl,
        N_genes      ~ px(!!(45  * GT_SCALE)),
        Top_genesets ~ px(!!(450 * GT_SCALE)),
        AD_genes     ~ px(!!(262 * GT_SCALE))
    )))() |>
    tab_style(
        style     = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) |>
    tab_style(
        style     = cell_text(align = "center"),
        locations = list(cells_body(columns = N_genes),
                        cells_column_labels(columns = N_genes))
    ) |>
    tab_options(
        table.font.names      = c("Arial", "Helvetica", "sans-serif"),
        table.font.size       = 11 * GT_SCALE,
        data_row.padding      = px(6 * GT_SCALE),
        column_labels.padding = px(8 * GT_SCALE)
    )

    out_path <- file.path(save_dir, paste0(out_string, ".html"))
    gtsave(gt_table, out_path)
    message("Saved: ", out_path)

    pdf_path <- file.path(save_dir, paste0(out_string, ".pdf"))
    tryCatch({
    gtsave(gt_table, pdf_path)
    message("Saved: ", pdf_path)
    }, error = function(e) {
    message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
    })

    bib_md   <- file.path(save_dir, paste0(out_string, ".md"))
    bib_docx <- file.path(save_dir, paste0(out_string, ".docx"))
    tryCatch({
    rmarkdown::pandoc_convert(bib_md, to = "docx", output = bib_docx)
    message("Saved: ", bib_docx)
    }, error = function(e) {
    message("Bibliography docx conversion failed (requires pandoc): ", conditionMessage(e))
    })

    # ── FDR-corrected version (MTG) ─────────────────────────────────────────────
    build_row_fdr <- function(ct) {
    genes   <- g_overlaps$Gene[g_overlaps$Celltype == ct]
    n_genes <- length(genes)

    ct_idx  <- which(celltypes == ct)
    pvals   <- gsea_out_fdr[[3 + ct_idx]]
    go_mask <- grepl("^GO", gsea_out_fdr$SetName)
    if (!is.null(pvals)) {
        go_pvals <- ifelse(go_mask, pvals, NA)
        sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
        n_sig    <- sum(sig_mask)
        if (n_sig == 0) {
        top_str <- "—"
        } else {
        top_idx  <- order(go_pvals)[1:min(3, n_sig)]
        top_sets <- gsea_out_fdr$SetName[top_idx]
        top_pval <- pvals[top_idx]
        top_str  <- paste0(
            seq_along(top_sets), ". ", top_sets,
            " (padj=", formatC(top_pval, format = "e", digits = 2), ")",
            collapse = "<br>"
        )
        }
    }

    genes_upper <- toupper(genes)
    hits        <- genes[genes_upper %in% names(ad_db)]
    hits_up     <- genes_upper[genes_upper %in% names(ad_db)]

    if (length(hits) > 0) {
        ref_nums <- sapply(hits_up, function(g) ad_db[[g]]$ref)
        ord      <- order(hits_up)
        hits_s   <- hits[ord]
        refs_s   <- ref_nums[ord]
        ad_str   <- paste0(hits_s, "<sup>", refs_s, "</sup>", collapse = ", ")
    } else {
        ad_str <- "—"
    }

    data.frame(
        Celltype     = ct,
        N_genes      = n_genes,
        Top_genesets = top_str,
        AD_genes     = ad_str,
        stringsAsFactors = FALSE
    )
    }

    summary_df_fdr <- do.call(rbind, lapply(celltypes, build_row_fdr)) |>
    arrange(Celltype) |>
    dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

    gt_table_fdr <- summary_df_fdr |>
    gt() |>
    cols_label(
        Celltype     = "Cell type",
        N_genes      = "Genes (n)",
        Top_genesets = "Top 3 GO gene sets (by FDR-adjusted p-value)",
        AD_genes     = column_4_header
    ) |>
    fmt_markdown(columns = AD_genes) |>
    fmt(columns = Top_genesets, fns = function(x) x) |>
    (\(tbl) rlang::inject(cols_width(tbl,
        N_genes      ~ px(!!(45  * GT_SCALE)),
        Top_genesets ~ px(!!(450 * GT_SCALE)),
        AD_genes     ~ px(!!(262 * GT_SCALE))
    )))() |>
    tab_style(
        style     = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) |>
    tab_style(
        style     = cell_text(align = "center"),
        locations = list(cells_body(columns = N_genes),
                        cells_column_labels(columns = N_genes))
    ) |>
    tab_options(
        table.font.names      = c("Arial", "Helvetica", "sans-serif"),
        table.font.size       = 11 * GT_SCALE,
        data_row.padding      = px(6 * GT_SCALE),
        column_labels.padding = px(8 * GT_SCALE)
    )

    out_path <- file.path(save_dir, paste0(out_string, ".html"))
    gtsave(gt_table_fdr, out_path)
    message("Saved: ", out_path)

    pdf_path <- file.path(save_dir, paste0(out_string, ".pdf"))
    tryCatch({
    gtsave(gt_table_fdr, pdf_path)
    message("Saved: ", pdf_path)
    }, error = function(e) {
    message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
    })
}