#!/bin/bash

# sh read-mapping.sh [read_len]

reads_out_dir="./reads/"
reads_file_prefix="reads"

ref_genome="./genome/reference-genome.fna"

# ==== Bowtie2 Read mapping ====

out_dir="./mappings/bowtie/"
file_prefix="bowtie"

# Build index for reference genome
bowtie2-build ${ref_genome} ${out_dir}genome-index

# Single mapping, best-match
bowtie2 -x ${out_dir}genome-index -U ${reads_out_dir}${reads_file_prefix}.fq -S ${out_dir}${file_prefix}-mapping-best-match.sam \
2> ${out_dir}${file_prefix}-mapping-best-match-log.txt

# Single mapping, report-all
bowtie2 -a -x ${out_dir}genome-index -U ${reads_out_dir}${reads_file_prefix}.fq -S ${out_dir}${file_prefix}-mapping-report-all.sam \
2> ${out_dir}${file_prefix}-mapping-report-all-log.txt

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

echo "\n=== Bowtie 2 SAMTools completed! ===\n"

# ==== BWA Read mapping ====

out_dir="./mappings/bwa/"
file_prefix="bwa"

# Build index for reference genome
bwa index -p ${out_dir}genome-index ${ref_genome}

# Single mapping, best-match
bwa mem ${out_dir}genome-index ${reads_out_dir}${reads_file_prefix}.fq > ${out_dir}${file_prefix}-mapping-best-match.sam \
2> ${out_dir}${file_prefix}-mapping-best-match-log.txt

# Single mapping, report-all
bwa mem -a ${out_dir}genome-index ${reads_out_dir}${reads_file_prefix}.fq > ${out_dir}${file_prefix}-mapping-report-all.sam \
2> ${out_dir}${file_prefix}-mapping-report-all-log.txt

echo "\n=== BWA read mapping completed! ===\n"

# ==== SAMTools on BWA ====

# Creating BAM files for best-match
samtools view -bS ${out_dir}${file_prefix}-mapping-best-match.sam -o ${out_dir}${file_prefix}-mapping-best-match.bam
samtools sort ${out_dir}${file_prefix}-mapping-best-match.bam -o ${out_dir}${file_prefix}-mapping-best-match-sorted.bam
samtools index ${out_dir}${file_prefix}-mapping-best-match-sorted.bam

# Creating BAM files for report-all
samtools view -bS ${out_dir}${file_prefix}-mapping-report-all.sam -o ${out_dir}${file_prefix}-mapping-report-all.bam
samtools sort ${out_dir}${file_prefix}-mapping-report-all.bam -o ${out_dir}${file_prefix}-mapping-report-all-sorted.bam
samtools index ${out_dir}${file_prefix}-mapping-report-all-sorted.bam

echo "\n=== BWA SAMTools completed! ===\n"


