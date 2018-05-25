#!/bin/bash
# Running other multi-mapping resolution methods

# sh multimapping-resolution.sh [read_len]

ref_genome="./genome/reference-genome.fna"

# ----------- Bowtie2 ----------
alignments="./mappings/bowtie/bowtie-mapping-report-all"
outfile="./mappings/bowtie/bowtie"

samtools sort -n ${alignments}.sam -o ${alignments}-id-sorted.sam
samtools view -bS ${alignments}-id-sorted.sam -o ${alignments}-id-sorted.bam

# MMR method

/usr/bin/time -v -o ${outfile}-mmr-time-log.txt mmr -o ${outfile}-mmr.bam -F 3 -p -b -R $1 ${alignments}-id-sorted.bam | tee ${outfile}-mmr-log.txt

sam_file=${outfile}-mmr
samtools sort ${sam_file}.bam -o ${sam_file}-sorted.bam
samtools index ${sam_file}-sorted.bam

echo "\n=== Bowtie + MMR multi-mapping resolution completed! ===\n"

# REMU method

/usr/bin/time -v -o ${outfile}-remu-time-log.txt remu.py -g ${ref_genome} -i ${alignments}-id-sorted.sam -o ${outfile}-remu.sam -r 10 | tee ${outfile}-remu-log.txt

sam_file=${outfile}-remu
samtools view -bS ${sam_file}.sam -o ${sam_file}.bam
samtools sort ${sam_file}.bam -o ${sam_file}-sorted.bam
samtools index ${sam_file}-sorted.bam

echo "\n=== Bowtie + REMU multi-mapping resolution completed! ===\n"

