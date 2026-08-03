# instructions for dl: https://anvilproject.org/learn/reference/gtex-v8-free-egress-instructions

# export all PFB
# - CF-GTEx ,bam files, RNA-seq
# filename: export_2023-03-30T08_15_19.avro
pfb to -i ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/GTEx/anvil_GTEx/export_2023-03-30T08_15_19.avro tsv

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/GTEx/anvil_GTEx/

./gen3-client configure --profile=gugenekang --cred=credentials.json --apiendpoint=https://gen3.theanvil.io

./gen3-client download-multiple --profile=gugenekang --manifest=file-manifest_filtered.json --download-path=${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/GTEx/fastq/ --protocol=s3


