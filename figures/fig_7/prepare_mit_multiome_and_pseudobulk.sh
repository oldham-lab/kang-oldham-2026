# python prepare_mit_multiome_and_pseudobulk.py \
#     --input         /home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025.h5ad \
#     --outdir        "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/PFC/pfc_data_pseudobulk/" \
#     --region-value  PFC

python prepare_mit_multiome_and_pseudobulk.py \
    --input         ${SHARED_DATA_DIR:-/mnt/bdata/@shared}/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025.h5ad \
    --outdir        "${SHARED_DATA_DIR:-/mnt/bdata/@shared}/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/MTC/mtc_data_pseudobulk/" \
    --region-value  MTC