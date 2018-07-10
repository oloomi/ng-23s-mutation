#!/bin/bash

script_path="$( cd "$(dirname "$0")" ; pwd -P )/"

ref_genome="./genome/reference-genome.fna"

out_path="./variants/"

file_path="./mappings/bowtie/"
alignment_files="bowtie-mapping-best-match-sorted
bowtie-mapping-report-all-sorted
bowtie-mmr-sorted
bowtie-prom-sorted"


for file in ${alignment_files}
do
	echo "\nVariant calling for: $file\n"
	freebayes -f ${ref_genome} -p 1 -F 0.9 -t ${script_path}genes.bed ${file_path}${file}.bam >${out_path}${file}-variants-freebayes.vcf
done

echo "\n=== Variant calling completed! ===\n"
