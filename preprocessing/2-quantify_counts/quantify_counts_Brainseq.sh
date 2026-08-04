# A log of all code used and steps taken in RNAseq preprocessing

# not all code is R code


# Started by downloading Brainseq data
python ~/code/git/RNAseq_preprocessing/1-sample_collection/download_Brainseq_from_synapse.py ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainseq/
# Total size: ~5.8 TB
 
# Run make_merge_fastq_script.py (in ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainseq) 
# Run merge_fastq.sh to merge separate fastq run in different lanes
  
  
# Pipeline (from Rebecca)
  
  ## QC:

#cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads
cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainseq/RNAseq_merged/


for ea in *.fastq.gz; do 
sem -j 8 fastqc $ea --outdir fastqc
done

#cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads/fastqc
cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainseq/RNAseq_merged/fastqc
multiqc .

## Trim:

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainseq/RNAseq_merged/

#ref="${SHARED_DIR:-/home/shared}/programs/bbmap/resources/adapters.fa"
ref="${BBMAP_DIR:-/home/gugene/bin/bbmap/bbmap}/resources/adapters.fa"

outdir="${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainseq/RNAseq_merged/trimmed"

#for R1 in XZ*R1*fastq.gz; do
for R1 in *R1*fastq.gz; do
R2=$(echo $R1 | sed "s/R1/R2/")
${BBMAP_DIR:-/home/gugene/bin/bbmap/bbmap}/bbduk.sh -Xmx1g \
in1=$R1 in2=$R2 \
out1=${outdir}/$(echo $R1 | sed "s/.fastq/_trimmed.fastq/") \
out2=${outdir}/$(echo $R2 | sed "s/.fastq/_trimmed.fastq/") \
ref=$ref t=10 ktrim=r minlength=60 entropy=0.01 maq=10 tpe tbo
done
# -Xmx1g: restrict to 1Gb of memory, fine for most settings
# t: number of threads
# ktrim=r: trim to right
# entropy: filter reads with entropy below value. range from 0 to 1. 0 = homopolymer(i.e. AAAA)
# maq: discard read with average quality below value
# tbo: trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences)
# tpe: trim both reads to the same length (in the event that an adapter kmer was only detected in one of them).

## QC:

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainseq/RNAseq_merged/trimmed

for ea in *fastq.gz; do
sem -j 16 fastqc $ea --outdir fastqc
done

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainseq/RNAseq_merged/trimmed/fastqc
multiqc .

## STAR & FEATURECOUNTS 

## Make index:
# 
# cd ${SHARED_DIR:-/home/shared}/hg_align_db/GRCm39_gencode_primary
# STAR --runMode genomeGenerate \
# --runThreadN 15 \
# --genomeDir star_index \
# --genomeFastaFiles GRCm39.primary_assembly.genome.fa \
# --sjdbGTFfile gencode.vM28.annotation.gtf \
# --sjdbOverhang 75
# 
# ## Align (v 2.7.8a)
# 
# cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads/trimmed
# 
# ref="${SHARED_DIR:-/home/shared}/hg_align_db/GRCm39_gencode_primary/star_index"
# outdir="/mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/STAR_aligned"
# 
# for R1 in XZ*R1*fastq.gz; do 
# R2=$(echo $R1 | sed "s/R1/R2/")
# STAR \
# --runThreadN 15 \
# --genomeDir $ref  \
# --readFilesCommand zcat \
# --outSAMtype BAM SortedByCoordinate \
# --readFilesIn $R1 $R2 \
# --outFilterType BySJout \
# --outFilterMismatchNmax 999 \
# --outFilterMismatchNoverReadLmax 0.04 \
# --alignIntronMax 1000000 \
# --outFilterMultimapNmax 20 \
# --outFilterScoreMinOverLread 0 \
# --outFilterMatchNminOverLread 0 \
# --outSAMmultNmax 1 \
# --sjdbOverhang 75 \
# --outFileNamePrefix ${outdir}/$(echo $R1 | sed "s/R1_001_trimmed.fastq.gz/trimmed_/")
# done
# 
# ## Align summary:
# 
# cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/STAR_aligned
# 
# annot="${SHARED_DIR:-/home/shared}/hg_align_db/GRCm39_gencode_primary/gencode.vM28.annotation.gtf"
# 
# qualimap="${SHARED_DIR:-/home/shared}/programs/qualimap_v2.2.1/./qualimap"
# 
# for ea in *bam; do
# $qualimap rnaseq \
# -a proportional \
# -outdir qualimap \
# -bam $ea -pe \
# -p strand-specific-forward \
# --java-mem-size=20G \
# -gtf $annot
# done
# 
# cd /home/rebecca/collabs/ghashghaei/bulk_expr/STAR_aligned
# multiqc .
# 
# ## Read quantification:
# 
# cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/STAR_aligned
# 
# annot="${SHARED_DIR:-/home/shared}/hg_align_db/GRCm39_gencode_primary/gencode.vM28.annotation.gtf"
# 
# featurecounts="/opt/subread-2.0.3-Linux-x86_64/bin/featureCounts"
# 
# for ea in *bam; do
# $featurecounts \
# -T 15 \
# -a $annot \
# -o ../featureCounts/$(echo $ea | sed "s/_Aligned.sortedByCoord.out.bam//") \
# -p -s 2 \
# -M -O --fraction \
# --countReadPairs \
# --extraAttributes gene_name \
# $ea
# done
# 
# cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/featureCounts
# multiqc .

## SALMON

## Make decoy for indexing:

## Cat genome to end of reference transcriptome:

#cd ${SHARED_DIR:-/home/shared}/hg_align_db/GRCm39_genome_plus_transcripts
#cat ../GRCm39_gencode_transcripts/gencode.vM28.transcripts.fa ../GRCm39_gencode_primary/GRCm39.primary_assembly.genome.fa > GRCm39.vM28.transcript_primary_assembly.genome.fa
mkdir ${DATA_DIR:-/mnt/bdata/gugene}/hg_align_db/GRCh38_genome_plus_transcripts
cd ${DATA_DIR:-/mnt/bdata/gugene}/hg_align_db/GRCh38_genome_plus_transcripts
cat ${SHARED_DIR:-/home/shared}/hg_align_db/GRCh38_gencode_transcripts/gencode.v39.transcripts.fa ${SHARED_DIR:-/home/shared}/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa > GRCh38.v39.transcript_primary_assembly.genome.fa

## Extract sequence headers from genome:

mkdir salmon_index
cd salmon_index
grep "^>" ${SHARED_DIR:-/home/shared}/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa | cut -d " " -f 1 > GRCh38.primary_assembly_decoys.txt
sed -i.bak -e 's/>//g' GRCh38.primary_assembly_decoys.txt

## Index:

conda activate salmon
salmon index -t ../GRCh38.v39.transcript_primary_assembly.genome.fa -i GRCh38.v39.transcript_primary_assembly.genome_Salmon_index --decoys GRCh38.primary_assembly_decoys.txt -k 31

## Align and quantify:

cd ${DATA_DIR:-/mnt/bdata/gugene}/datasets/RNAseq/Brainseq/RNAseq_merged/trimmed

ref="${DATA_DIR:-/mnt/bdata/gugene}/hg_align_db/GRCh38_genome_plus_transcripts/salmon_index/GRCh38.v39.transcript_primary_assembly.genome_Salmon_index"
annot="${SHARED_DIR:-/home/shared}/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf"

conda activate salmon

for R1 in *R1*; do
R2=$(echo $R1 | sed "s/R1/R2/")
salmon quant -i $ref \
-l A \
-g $annot \
-1 $R1 -2 $R2 \
-p 15 \
--validateMappings \
--gcBias \
--seqBias \
-o ../../Salmon_aligned/decoy/$(echo $R1 | sed "s/R1_trimmed.fastq.gz/trimmed/")_transcripts_quant
done

# -p: number of threads
# --gcBias: correct for gc bias
# --seqBias: Specifically, this model will attempt to correct for random hexamer priming bias, which results in the preferential sequencing of fragments starting with certain nucleotide motifs. 
# -l A: automatically infer library type
# was previously -l ISR