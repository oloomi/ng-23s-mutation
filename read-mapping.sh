#!/bin/bash

# sh read-mapping.sh

reads_out_dir="./reads/"
reads_file_prefix="reads"

ref_genome="./genome/reference-genome.fna"

# ==== Bowtie2 Read mapping ====

out_dir="./mappings/bowtie/"
file_prefix="bowtie"

# Build index for reference genome
bowtie2-build ${ref_genome} ${out_dir}genome-index

# Single mapping, best-match
bowtie2 -x ${out_dir}genome-index -U ${reads_out_dir}${reads_file_prefix}_1.fq,${reads_out_dir}${reads_file_prefix}_2.fq \
-S ${out_dir}${file_prefix}-mapping-best-match.sam 2> ${out_dir}${file_prefix}-mapping-best-match-log.txt

# Single mapping, report-all
bowtie2 -a -x ${out_dir}genome-index -U ${reads_out_dir}${reads_file_prefix}_1.fq,${reads_out_dir}${reads_file_prefix}_2.fq \
-S ${out_dir}${file_prefix}-mapping-report-all.sam 2> ${out_dir}${file_prefix}-mapping-report-all-log.txt

# Build index for reference genome, masked
ref_genome="./genome/reference-genome-masked.fna"
bowtie2-build ${ref_genome} ${out_dir}masked-genome-index

# Single mapping, best-match, masked genome
bowtie2 -x ${out_dir}masked-genome-index -U ${reads_out_dir}${reads_file_prefix}_1.fq,${reads_out_dir}${reads_file_prefix}_2.fq \
-S ${out_dir}${file_prefix}-mapping-best-match-masked.sam 2> ${out_dir}${file_prefix}-mapping-best-match-masked-log.txt

# Paired mapping, best-match
#bowtie2 -x ${out_dir}genome-index -1 ${reads_out_dir}${reads_file_prefix}_1.fq -2 ${reads_out_dir}${reads_file_prefix}_2.fq \
#-S ${out_dir}${file_prefix}-mapping-best-match.sam 2> ${out_dir}${file_prefix}-mapping-best-match-log.txt

# Paired mapping, report-all
#bowtie2 -a -x ${out_dir}genome-index -1 ${reads_out_dir}${reads_file_prefix}_1.fq -2 ${reads_out_dir}${reads_file_prefix}_2.fq \
#-S ${out_dir}${file_prefix}-mapping-report-all.sam 2> ${out_dir}${file_prefix}-mapping-report-all-log.txt

echo "\n=== Bowtie2 read mapping completed! ===\n"

# ==== SAMTools on Bowtie2 ====

# Creating BAM files for best-match
samtools view -bS ${out_dir}${file_prefix}-mapping-best-match.sam -o ${out_dir}${file_prefix}-mapping-best-match.bam
samtools sort ${out_dir}${file_prefix}-mapping-best-match.bam -o ${out_dir}${file_prefix}-mapping-best-match-sorted.bam
samtools index ${out_dir}${file_prefix}-mapping-best-match-sorted.bam

# Creating BAM files for report-all
samtools view -bS ${out_dir}${file_prefix}-mapping-report-all.sam -o ${out_dir}${file_prefix}-mapping-report-all.bam
samtools sort ${out_dir}${file_prefix}-mapping-report-all.bam -o ${out_dir}${file_prefix}-mapping-report-all-sorted.bam
samtools index ${out_dir}${file_prefix}-mapping-report-all-sorted.bam

# Creating BAM files for best-match, masked
samtools view -bS ${out_dir}${file_prefix}-mapping-best-match-masked.sam -o ${out_dir}${file_prefix}-mapping-best-match-masked.bam
samtools sort ${out_dir}${file_prefix}-mapping-best-match-masked.bam -o ${out_dir}${file_prefix}-mapping-best-match-masked-sorted.bam
samtools index ${out_dir}${file_prefix}-mapping-best-match-masked-sorted.bam

echo "\n=== Bowtie 2 SAMTools completed! ===\n"
