#!/bin/bash

# sh experiment.sh [ref_genome] [SRA] [read_len] [genes]

script_path="$( cd "$(dirname "$0")" ; pwd -P )/"

make_exp_dir() {
  mkdir -p $1/genome $1/mappings/bowtie $1/mappings/bwa $1/reads $1/results $1/variants
}

get_genome() {
  wget -O - $1 | gunzip -c > ./genome/reference-genome.fna
}

get_reads() {
  fastq-dump --outdir $3/reads/fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-files --clip $1
  gunzip $3/reads/fastq/$1_pass_1.fastq.gz
  gunzip $3/reads/fastq/$1_pass_2.fastq.gz
  fastx_trimmer -l $2 -m $2 -Q33 -i $3/reads/fastq/$1_pass_1.fastq -o $3/reads/reads_1.fq
  fastx_trimmer -l $2 -m $2 -Q33 -i $3/reads/fastq/$1_pass_2.fastq -o $3/reads/reads_2.fq
}


# Experiment

make_exp_dir .
get_genome $1
get_reads $2 $3 .

sh ${script_path}read-mapping.sh
sh ${script_path}multimapping-resolution.sh $3
sh ${script_path}variant-calling.sh

python3 ${script_path}evaluation.py
