python calc_projIndexSE_brainSCOPE_v7.py \
      --metadata   ${SHARED_DATA_DIR:-/mnt/bdata/@shared}/scsn.expr_data/human_expr/postnatal/brainSCOPE/PEC2_sample_metadata.txt \
      --matrix_dir ${SHARED_DATA_DIR:-/mnt/bdata/@shared}/scsn.expr_data/human_expr/postnatal/brainSCOPE/snrna_expr_matrices/ \
      --modules    ${DATA_DIR:-${DATA_DIR:-/mnt/bdata/gugene}}/data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/topmodposbc_table.json \
      --output_dir ${DATA_DIR:-${DATA_DIR:-/mnt/bdata/gugene}}/data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/brainSCOPE/sn_proj_indices/log_native/CMC_SZBD_SE \
      --normalize --verbose
