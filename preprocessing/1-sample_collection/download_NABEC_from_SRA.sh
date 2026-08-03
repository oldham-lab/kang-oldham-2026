cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/NABEC/

prefetch --option-file SRR_Acc_List.txt --ngc ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/prj_9519_D37443.ngc


cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/NABEC

time for F in SRR5*; do
cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/NABEC
cd $F
fasterq-dump --ngc ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/prj_9519_D37443.ngc \
$(echo $F).sra
done