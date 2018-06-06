#!/bin/bash
# sh experiment.sh [ref_genome] [SRA] [read_len] [genes]

# Change to the path to Trimmomatic on your computer
export TRIMMOMATIC=/home/ubuntu/tools/Trimmomatic-0.38/trimmomatic-0.38.jar

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
  java -jar $TRIMMOMATIC PE -phred33 $3/reads/fastq/$1_pass_1.fastq $3/reads/fastq/$1_pass_2.fastq \
  $3/reads/reads_1.fq $3/reads/output_forward_unpaired.fastq $3/reads/reads_2.fq $3/reads/output_reverse_unpaired.fastq \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:75
#  fastx_trimmer -l $2 -m $2 -Q33 -i $3/reads/fastq/$1_pass_1.fastq -o $3/reads/reads_1.fq
#  fastx_trimmer -l $2 -m $2 -Q33 -i $3/reads/fastq/$1_pass_2.fastq -o $3/reads/reads_2.fq
}


# Experiment

make_exp_dir .
get_genome $1
sync
python3 ${script_path}mask_genome.py ./genome/reference-genome.fna $4
get_reads $2 $3 .
sync

sh ${script_path}read-mapping.sh
sh ${script_path}multimapping-resolution.sh $3
sh ${script_path}variant-calling.sh

python3 ${script_path}evaluation.py $4
