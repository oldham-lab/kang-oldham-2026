cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/GTEx/bam/

tra="${SHARED_DIR:-/home/shared}/hg_align_db/GRCh38_gencode_transcripts/gencode.v39.transcripts.fa"

# First convert BAM files to FASTQ, then run Kallisto pipeline as we have done for other datasets.
# Pipeline detailed at: https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq

docker pull broadinstitute/gtex_rnaseq:V10

path_to_data=${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/GTEx/bam/

for input_bam in *.bam; do
sample_id=$(echo $input_bam | sed "s/.Aligned.sortedByCoord.out.patched.md.bam//")
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq \
    /bin/bash -c "/src/run_SamToFastq.py /data/$input_bam -p $sample_id -o /data"
done

# There is one sample ("GTEX-1H4P4-0011-R10b-SM-CE6SG") with unpaired reads.

# run kallisto

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/GTEx/fastq/

ref="${DATA_DIR:-/mnt/bdata/gugene}/hg_align_db/kallisto/index/homo_sapiens_index"  
    
time for R1 in *1.fastq.gz; do
R2=$(echo $R1 | sed "s/1\.fastq\.gz/2\.fastq\.gz/")
kallisto quant -i $ref \
-o $(echo $R1 | sed "s/1\.fastq\.gz/trimmed/")_transcripts_quant_kal \
$R1 \
$R2 \
-t 16 \
--bias
done

mkdir kal_quant
mv *quant_kal kal_quant


