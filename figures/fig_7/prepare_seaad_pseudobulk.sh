# # MTG
# python prepare_seaad_pseudobulk.py \
#     --input  ${DATA_DIR:-${DATA_DIR:-/mnt/bdata/gugene}}/datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
#     --outdir ${DATA_DIR:-${DATA_DIR:-/mnt/bdata/gugene}}/datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/mtg_data_pseudobulk/ \
#     --outprefix sea_mtg

# DFC
python prepare_seaad_pseudobulk.py \
    --input  ${DATA_DIR:-${DATA_DIR:-/mnt/bdata/gugene}}/datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad \
    --outdir ${DATA_DIR:-${DATA_DIR:-/mnt/bdata/gugene}}/datasets/SN_RNAseq/sea_ad_2024/RNAseq/dfc_data_pseudobulk/ \
    --outprefix sea_dfc
