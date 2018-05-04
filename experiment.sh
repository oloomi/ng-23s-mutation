#!/bin/bash

make_exp_dir() {
  mkdir -p $1/genome $1/mappings/bowtie $1/mappings/bwa $1/reads $1/results $1/variants
}

get_reads() {
  fastq-dump --outdir $3/reads/fastq --gzip --skip-technical --read-filter pass --dumpbase --split-files --clip $1
  gunzip $3/reads/fastq/$1_pass_1.fastq.gz
  gunzip $3/reads/fastq/$1_pass_2.fastq.gz
  fastx_trimmer -l $2 -m $2 -Q33 -i $3/reads/fastq/$1_pass_1.fastq -o $3/reads/reads_1.fq
  fastx_trimmer -l $2 -m $2 -Q33 -i $3/reads/fastq/$1_pass_2.fastq -o $3/reads/reads_2.fq
}


# Experiment

make_exp_dir .

# reference genome Neisseria gonorrhoeae strain NCCP11945
# https://www.ncbi.nlm.nih.gov/nuccore/CP001050.1
# copy to genome folder

get_reads SRR5827361 130 .

sh read-mapping.sh
sh multimapping-resolution.sh 130
sh variant-calling.sh
