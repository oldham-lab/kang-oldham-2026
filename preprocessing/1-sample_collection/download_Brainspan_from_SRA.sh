cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainspan/

prefetch --option-file SRR_Acc_List_filtered.txt --ngc ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/prj_9519_D37443.ngc

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainspan

time for F in SRR35*; do
cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainspan
cd $F
fasterq-dump --ngc ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/prj_9519_D37443.ngc \
$(echo $F).sra
done