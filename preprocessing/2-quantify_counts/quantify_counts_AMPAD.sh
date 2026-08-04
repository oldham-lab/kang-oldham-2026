# A log of all code used and steps taken in RNAseq preprocessing

##############
# ROSMAP
##############
# Pipeline (from Rebecca)
## QC:
cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/ROSMAP/fastq/
mkdir fastqc

for ea in *.fastq.gz; do 
sem -j 16 fastqc $ea --outdir fastqc
done

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/ROSMAP/fastq/fastqc
multiqc .

## Trim:

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/ROSMAP/fastq/

#ref="${SHARED_DIR:-/home/shared}/programs/bbmap/resources/adapters.fa"
ref="${BBMAP_DIR:-/home/gugene/bin/bbmap/bbmap}/resources/adapters.fa"

outdir="${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed"

for R1 in *r1*fastq.gz; do
R2=$(echo $R1 | sed "s/r1/r2/")
${BBMAP_DIR:-/home/gugene/bin/bbmap/bbmap}/bbduk.sh -Xmx1g \
in1=$R1 in2=$R2 \
out1=${outdir}/$(echo $R1 | sed "s/.fastq/_trimmed.fastq/") \
out2=${outdir}/$(echo $R2 | sed "s/.fastq/_trimmed.fastq/") \
ref=$ref t=16 ktrim=r minlength=60 entropy=0.01 maq=10 tpe tbo
done
# -Xmx1g: restrict to 1Gb of memory, fine for most settings
# t: number of threads
# ktrim=r: trim to right
# entropy: filter reads with entropy below value. range from 0 to 1. 0 = homopolymer(i.e. AAAA)
# maq: discard read with average quality below value
# tbo: trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences)
# tpe: trim both reads to the same length (in the event that an adapter kmer was only detected in one of them).

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed
mkdir fastqc

for ea in *fastq.gz; do
sem -j 16 fastqc $ea --outdir fastqc
done

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed/fastqc
multiqc .

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed/

ref="${DATA_DIR:-/mnt/bdata/gugene}/hg_align_db/kallisto/index/homo_sapiens_index"  
    
time for R1 in *r1*; do
R2=$(echo $R1 | sed "s/r1/r2/")
kallisto quant -i $ref \
-o $(echo $R1 | sed "s/r1_trimmed.fastq.gz/trimmed/")_transcripts_quant_kal \
$R1 \
$R2 \
-t 16 \
--bias
done

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed/
mkdir kal_quant
mv *quant* kal_quant

##############
# MSBB
##############
# cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/MSBB/fastq/
# mkdir fastqc
# 
# for ea in *.fastq.gz; do 
# sem -j 16 fastqc $ea --outdir fastqc
# done
# 
# cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/MSBB/fastq/fastqc
# multiqc .

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/MSBB/fastq/

ref="${BBMAP_DIR:-/home/gugene/bin/bbmap/bbmap}/resources/adapters.fa"

outdir="${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/MSBB/fastq/trimmed"
mkdir $outdir

for R1 in *fastq.gz; do
${BBMAP_DIR:-/home/gugene/bin/bbmap/bbmap}/bbduk.sh -Xmx1g \
in=$R1 \
out1=${outdir}/$(echo $R1 | sed "s/fastq/1_trimmed.fastq/") \
out2=${outdir}/$(echo $R1 | sed "s/fastq/2_trimmed.fastq/") \
ref=$ref t=16 ktrim=r minlength=60 entropy=0.01 maq=10 tpe tbo
done

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/MSBB/fastq/trimmed/
mkdir fastqc
 
for ea in *.fastq.gz; do 
sem -j 16 fastqc $ea --outdir fastqc
done
 
cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/MSBB/fastq/trimmed/fastqc
multiqc .

# Kallisto

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/MSBB/fastq/trimmed/

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

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/AMPAD/MSBB/fastq/trimmed/
mkdir kal_quant
mv *transcripts_quant_kal kal_quant
