# River (alluvial) plots: Liu RNA.Subclass → Gabitto_metacell_labels
# Requires: ggplot2, ggalluvial
#
# Install if needed:
#   install.packages(c("ggplot2", "ggalluvial"))

library(ggplot2)
library(ggalluvial)
library(ggrepel)
library(patchwork)
library(showtext)
showtext_auto()

# Generate n colors with maximally separated hues for adjacent categories.
# Strategy: evenly space hues around the color wheel, then reorder by a
# stride of floor(n/2)+1 so no two neighboring strata share a similar hue.
make_interleaved_colors <- function(n) {
    hues    <- seq(0, 360 * (1 - 1/n), length.out = n)
    stride  <- floor(n / 2) + 1
    idx     <- (seq(0, n - 1) * stride) %% n + 1
    hcl(h = hues[idx], c = 70, l = 60)
}

SAVE_DIR  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s6")
PLOT_W    <- 10 * 0.7   # inches — 30 % narrower than original
PLOT_H    <- 16         # inches
TEXT_SIZE <- 3          # geom_text size (ggplot2 units ≈ mm)
PLOT_W_REF <- 10        # original width stratum_width formula was calibrated for

csv_files <- list.files(SAVE_DIR, pattern = "^liu_gabitto_metacell_labels_.*\\.csv$", full.names = TRUE)

for (path in csv_files) {
    label <- sub("^liu_gabitto_metacell_labels_", "", basename(path))
    label <- sub("\\.csv$", "", label)
    message("\n── Processing: ", label, " ──")

    df <- read.csv(path, row.names = 1)
    df$RNA.Subclass             <- sub("^(Exc|Inh)[[:space:]_-]*", "", df$RNA.Subclass)
    df$RNA.Subclass             <- paste0(df$RNA.Subclass, " (Liu)")
    df$Gabitto_metacell_labels  <- paste0(df$Gabitto_metacell_labels, " (Gab)")

    # Detect metric column (correlation only)
    metric_col  <- "Gabitto_metacell_max_corr"
    metric_name <- "mean_max_corr"
    # else {
    #     metric_col  <- "Gabitto_metacell_min_dist"
    #     metric_name <- "mean_min_dist"
    # }

    # Aggregate counts per (RNA.Subclass, Gabitto_metacell_labels) pair
    agg <- as.data.frame(table(
        Liu_Subclass     = df$RNA.Subclass,
        Gabitto_Metacell = df$Gabitto_metacell_labels
    ))
    agg <- agg[agg$Freq > 0, ]

    subclasses  <- sort(unique(agg$Liu_Subclass),     decreasing = TRUE)
    metacells   <- sort(unique(agg$Gabitto_Metacell), decreasing = TRUE)
    color_map   <- setNames(make_interleaved_colors(length(subclasses)), subclasses)

    agg$Liu_Subclass     <- factor(agg$Liu_Subclass,     levels = subclasses)
    agg$Gabitto_Metacell <- factor(agg$Gabitto_Metacell, levels = metacells)

    # Uniform stratum width scaled to the longest label across both columns.
    # Scaled inversely with PLOT_W so physical bar width stays constant.
    stratum_width <- max(nchar(as.character(c(subclasses, metacells)))) * 0.012 * (PLOT_W_REF / PLOT_W)

    # ── Step 1: probe plot to extract stratum positions ───────────────────────
    p_probe <- ggplot(agg,
                aes(axis1 = Liu_Subclass, axis2 = Gabitto_Metacell, y = Freq)) +
        geom_alluvium(aes(fill = Liu_Subclass), width = stratum_width, alpha = 0.8, reverse = FALSE, discern = TRUE) +
        geom_stratum(width = stratum_width, fill = "grey85", color = "grey40", reverse = FALSE, discern = TRUE) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = TEXT_SIZE, reverse = FALSE, discern = TRUE) +
        scale_x_discrete(limits = c("Liu_Subclass", "Gabitto_Metacell"),
                         expand = c(0.3, 0.3)) +
        scale_fill_manual(values = color_map) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none", panel.grid = element_blank())

    built      <- ggplot_build(p_probe)
    strat_data <- built$data[[2]]   # geom_stratum layer → xmin/xmax/ymin/ymax
    label_data <- built$data[[3]]   # geom_text layer    → x/y/label

    label_data$height <- strat_data$ymax - strat_data$ymin
    label_data$xmin   <- strat_data$xmin
    label_data$xmax   <- strat_data$xmax

    # Threshold: stratum height < font height in data units → outside label
    y_data_range <- max(strat_data$ymax) - min(strat_data$ymin)
    panel_h_mm   <- PLOT_H * 25.4 * 0.82
    thresh       <- TEXT_SIZE / panel_h_mm * y_data_range
    message(sprintf("  Font-height threshold: %.1f data units (y range = %.1f)", thresh, y_data_range))

    inside_df  <- label_data[label_data$height >= thresh, ]
    outside_df <- label_data[label_data$height <  thresh, ]

    outside_df$pt_x    <- ifelse(outside_df$x < 1.5, outside_df$xmin, outside_df$xmax)
    outside_df$nudge_x <- ifelse(outside_df$x < 1.5, -0.08, 0.08)
    outside_df$hjust   <- ifelse(outside_df$x < 1.5, 1, 0)

    # ── Step 2: final plot ────────────────────────────────────────────────────
    p <- ggplot(agg,
                aes(axis1 = Liu_Subclass, axis2 = Gabitto_Metacell, y = Freq)) +
        geom_alluvium(aes(fill = Liu_Subclass), width = stratum_width, alpha = 0.8, reverse = FALSE, discern = TRUE) +
        geom_stratum(width = stratum_width, fill = "grey85", color = "grey40", reverse = FALSE, discern = TRUE) +
        geom_text(data        = inside_df,
                  aes(x = x, y = y, label = label),
                  size        = TEXT_SIZE,
                  inherit.aes = FALSE) +
        geom_text_repel(data               = outside_df,
                        aes(x = pt_x, y = y, label = label, hjust = hjust),
                        nudge_x            = outside_df$nudge_x,
                        direction          = "y",
                        size               = TEXT_SIZE,
                        segment.size       = 0.3,
                        segment.color      = "grey40",
                        min.segment.length = 0,
                        inherit.aes        = FALSE) +
        scale_x_discrete(limits = c("Liu_Subclass", "Gabitto_Metacell"),
                         expand = c(0.3, 0.3)) +
        scale_fill_manual(values = color_map) +
        scale_y_continuous(labels = scales::comma) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none",
              panel.grid      = element_blank(),
              axis.text       = element_blank(),
              axis.ticks      = element_blank(),
              axis.title      = element_blank(),
              plot.title      = element_blank())

    out_pdf <- file.path(SAVE_DIR, paste0("river_plot_", label, ".pdf"))
    out_svg <- file.path(SAVE_DIR, paste0("river_plot_", label, ".svg"))
    ggsave(out_pdf, p, width = PLOT_W, height = PLOT_H)
    ggsave(out_svg, p, width = PLOT_W, height = PLOT_H)
    message("Saved: ", out_pdf)
    message("Saved: ", out_svg)

    # ── Summary table ──────────────────────────────────────────────────────────
    # % of each Liu subclass assigned to each Gabitto metacell
    pct_wide <- reshape(
        as.data.frame(with(df, prop.table(table(RNA.Subclass, Gabitto_metacell_labels), margin = 1) * 100)),
        idvar     = "RNA.Subclass",
        timevar   = "Gabitto_metacell_labels",
        direction = "wide"
    )
    colnames(pct_wide) <- sub("^Freq\\.", "", colnames(pct_wide))
    pct_wide <- pct_wide[order(pct_wide$RNA.Subclass), ]

    # Mean metric per Liu subclass
    mean_metric <- aggregate(df[[metric_col]] ~ df$RNA.Subclass, FUN = mean)
    colnames(mean_metric) <- c("RNA.Subclass", metric_name)
    mean_metric[[metric_name]] <- round(mean_metric[[metric_name]], 4)

    summary_tbl <- merge(pct_wide, mean_metric, by = "RNA.Subclass", sort = TRUE)

    # Footer row: mean metric per Gabitto metacell
    mc_mean_metric <- aggregate(df[[metric_col]] ~ df$Gabitto_metacell_labels, FUN = mean)
    colnames(mc_mean_metric) <- c("Gabitto_metacell_labels", metric_col)
    footer <- data.frame(RNA.Subclass = metric_name)
    footer[[metric_name]] <- NA
    for (mc in mc_mean_metric$Gabitto_metacell_labels) {
        footer[[mc]] <- round(mc_mean_metric[[metric_col]][mc_mean_metric$Gabitto_metacell_labels == mc], 4)
    }
    # Fill any metacell columns missing from footer with NA
    for (col in setdiff(colnames(summary_tbl), colnames(footer))) {
        footer[[col]] <- NA
    }
    summary_tbl <- rbind(summary_tbl, footer[, colnames(summary_tbl)])

    tbl_path <- file.path(SAVE_DIR, paste0("summary_table_", label, ".csv"))
    write.csv(summary_tbl, tbl_path, row.names = FALSE)
    message("Saved: ", tbl_path)

    # ── Barplots: % Liu subclass  /  % assigned to each Gabitto metacell ───────
    liu_pct <- as.data.frame(prop.table(table(RNA.Subclass = df$RNA.Subclass)) * 100)
    liu_pct$RNA.Subclass <- factor(liu_pct$RNA.Subclass,
                                   levels = sort(unique(as.character(liu_pct$RNA.Subclass))))

    mc_pct <- as.data.frame(prop.table(table(Gabitto_Metacell = df$Gabitto_metacell_labels)) * 100)
    mc_pct$Gabitto_Metacell <- factor(mc_pct$Gabitto_Metacell,
                                      levels = sort(unique(as.character(mc_pct$Gabitto_Metacell))))

    bar_theme <- theme_minimal(base_size = 12) +
        theme(axis.text.x  = element_text(angle = 45, hjust = 1),
              panel.grid.major.x = element_blank())

    p_liu <- ggplot(liu_pct, aes(x = RNA.Subclass, y = Freq, fill = RNA.Subclass)) +
        geom_col(show.legend = FALSE) +
        scale_fill_manual(values = setNames(
            make_interleaved_colors(nlevels(liu_pct$RNA.Subclass)),
            levels(liu_pct$RNA.Subclass)
        )) +
        labs(x = NULL, y = "% of Liu cells", title = "Liu subclass composition") +
        bar_theme

    p_mc <- ggplot(mc_pct, aes(x = Gabitto_Metacell, y = Freq, fill = Gabitto_Metacell)) +
        geom_col(show.legend = FALSE) +
        scale_fill_manual(values = setNames(
            make_interleaved_colors(nlevels(mc_pct$Gabitto_Metacell)),
            levels(mc_pct$Gabitto_Metacell)
        )) +
        labs(x = NULL, y = "% of Liu cells", title = "Gabitto metacell assignments") +
        bar_theme

    bar_plot <- p_liu / p_mc +
        plot_annotation(title = paste0("Cell type distributions (", label, ")"))

    bar_pdf <- file.path(SAVE_DIR, paste0("barplot_", label, ".pdf"))
    bar_svg <- file.path(SAVE_DIR, paste0("barplot_", label, ".svg"))
    ggsave(bar_pdf, bar_plot, width = 10, height = 8)
    ggsave(bar_svg, bar_plot, width = 10, height = 8)
    message("Saved: ", bar_pdf)
    message("Saved: ", bar_svg)
}
