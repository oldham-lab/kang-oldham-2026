# move fastq files to same folder

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/PsychENCODE/RNAseq/Yale_ASD/
ref="/home/gugene/bin/bbmap/bbmap/resources/adapters.fa"
outdir="${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/PsychENCODE/RNAseq/Yale_ASD/trimmed"

for R1 in *1_001.fastq.gz; do
R2=$(echo $R1 | sed "s/1_001.fastq.gz/2_001.fastq.gz/")
/home/gugene/bin/bbmap/bbmap/bbduk.sh -Xmx1g \
in1=$R1 in2=$R2 \
out1=${outdir}/$(echo $R1 | sed "s/.fastq/_trimmed.fastq/") \
out2=${outdir}/$(echo $R2 | sed "s/.fastq/_trimmed.fastq/") \
ref=$ref t=10 ktrim=r minlength=60 entropy=0.01 maq=10 tpe tbo
done
  
# FastQC, Multiqc  

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/PsychENCODE/RNAseq/BrainGVEX/trimmed
mkdir fastqc

time for ea in *fastq; do
sem -j 16 fastqc $ea --outdir fastqc
done

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/PsychENCODE/RNAseq/BrainGVEX/trimmed/fastqc
time multiqc .

# Kallisto

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/PsychENCODE/RNAseq/BrainGVEX/trimmed/

ref="${DATA_DIR:-/mnt/bdata/gugene}/hg_align_db/kallisto/index/homo_sapiens_index"  
    
time for R1 in *1_trimmed*; do
R2=$(echo $R1 | sed "s/1_trimmed/2_trimmed/")
kallisto quant -i $ref \
-o $(echo $R1 | sed "s/1_trimmed.fastq/trimmed/")_transcripts_quant_kal \
$R1 \
$R2 \
-t 16 \
--bias
done

mkdir kal_quant
mv *quant_kal kal_quant

