# Central path configuration for the kang-oldham-2026 analysis repo.
#
# Scripts read these roots inline via Sys.getenv("<VAR>", "<default>"), so each
# script runs standalone (defaults = the author's original paths) and every root
# is overridable through a local .Renviron (see .Renviron.example). This file
# mirrors those defaults for interactive use / documentation.
#
# Note: /home/gugene/bdata is a symlink to /mnt/bdata, and ~ == /home/gugene,
# so the bdata data is one tree, reached here via DATA_DIR.

# --- repo self (this checkout) ---
REPO_DIR       <- Sys.getenv("REPO_DIR",       unset = normalizePath(getwd(), mustWork = FALSE))

# --- data roots (NOT distributed with the repo; see DATA.md) ---
DATA_DIR       <- Sys.getenv("DATA_DIR", unset = "/mnt/bdata/gugene")            # datasets/, data/greedy_march_pipeline_output/, figures/
MEGASET_DIR    <- Sys.getenv("MEGASET_DIR", unset = "/home/gugene/RNAseq_megaset")  # pseudobulk, ESSA_data, gene maps
DATA_OTHER_DIR <- Sys.getenv("DATA_OTHER_DIR", unset = "/home/gugene/data_other")      # lein DFC, imputed data
EXTFIG_DIR     <- Sys.getenv("EXTFIG_DIR", unset = "/home/gugene/figures")         # figure-1 followup outputs (outside repo)

# --- shared / external tool + data locations (documented, not vendored) ---
SHARED_DIR     <- Sys.getenv("SHARED_DIR", unset = "/home/shared")                 # scsn.expr_data, hg_align_db, code/SampleNetworks
SHARED_DATA_DIR <- Sys.getenv("SHARED_DATA_DIR", unset = "/mnt/bdata/@shared")       # lab-shared scsn.expr_data/human_expr (single-nucleus reference datasets)
ESSA_DIR       <- Sys.getenv("ESSA_DIR", unset = "/home/gugene/code/git/ESSA")   # ESSA DE repo
PSEUDOBULK_DIR <- Sys.getenv("PSEUDOBULK_DIR", unset = "/home/gugene/code/git/Pseudobulk-from-SC-SN-data")
HOME_DATA_DIR  <- Sys.getenv("HOME_DATA_DIR", unset = "/home/gugene/data")          # a few one-off lein_mtg / ABI cell-count inputs

# Note: SHARED_DIR and SHARED_DATA_DIR both contain a `scsn.expr_data/` tree but
# are DIFFERENT directories (neither is a symlink to the other) — do not merge them.

# --- external analysis tools sourced by path (not R packages; see README) ---
FINDMODULES_DIR    <- Sys.getenv("FINDMODULES_DIR", unset = "/home/gugene/code/git/FindModules")        # FindModules_2.0; scripts source FindModules/R/FindModules.R
GSEA_GENERIC_DIR   <- Sys.getenv("GSEA_GENERIC_DIR", unset = "/home/gugene/code/git/GSEA_generic")      # GSEAfxsV3.r / GSEAfxsV3_nonpar_temp.r
SAMPLENETWORK_DIR  <- Sys.getenv("SAMPLENETWORK_DIR", unset = "/home/gugene/code/labcode_old/SampleNetwork")  # SampleNetwork_1.08.r -- see DATA.md, two non-identical copies exist
SCINRB_DIR         <- Sys.getenv("SCINRB_DIR", unset = "/home/gugene/code/git_other/scINRB")           # scINRB imputation (fig_2 related analyses only)
COMPAREMARKERS_DIR <- Sys.getenv("COMPAREMARKERS_DIR", unset = "/home/gugene/code/git/CompareMarkers")  # bulk fidelity table (fig1 followup sensitivity)
BBMAP_DIR          <- Sys.getenv("BBMAP_DIR", unset = "/home/gugene/bin/bbmap/bbmap")                   # bbduk.sh + resources/adapters.fa (read trimming)
PYTHON_BIN         <- Sys.getenv("PYTHON_BIN", unset = "/home/gugene/miniconda3/bin/python")            # interpreter invoked from R for SVG/PPTX rendering

# --- CoPA Shiny app checkout (snapshot generation only, not part of any figure) ---
SHINYAPP_DIR   <- Sys.getenv("SHINYAPP_DIR", unset = "/home/gugene/ShinyApps/copacabana")

# --- throwaway output dir used by exploratory scripts (safe to point anywhere) ---
SCRATCH_DIR    <- Sys.getenv("SCRATCH_DIR", unset = "~/test")

# --- large reference files not shipped in CoPA (download separately) ---
MSIGDB_XML     <- Sys.getenv("MSIGDB_XML", unset = "/home/gugene/code/git/CoPA/data-raw/msigdb_v7.4.xml")  # MSigDB XML for Broad GSEA
HGNC_MAP       <- Sys.getenv("HGNC_MAP",   unset = "/home/gugene/code/git/CoPA/inst/python/dataset_processing_unified/hgnc_prev_to_current.tsv")  # HGNC prev->current symbol map

# --- cell-type AD-vs-control DE result objects (fig_7; ~140MB, NOT tracked) ---
# Regenerate on demand: figures/fig_7/v4/full_DE_pipeline_ADvsCon.R then
# figures/fig_7/v4/DE/concat_DE_results.R write here; fig_7 panels read here.
# Default is a gitignored in-repo dir. See DATA.md.
DE_DIR         <- Sys.getenv("DE_DIR",         unset = file.path(REPO_DIR, "figures/fig_7/DE_old"))

# Method code is a package, not sourced by path:  library(CoPA)
# Shared plot theme is vendored in this repo:      utils/ggplot_theme_settings.R

for (.d in c(DATA_DIR = DATA_DIR, MEGASET_DIR = MEGASET_DIR)) {
  if (!dir.exists(.d)) warning(sprintf("data root does not exist: %s (set it in .Renviron)", .d))
}
rm(.d)
