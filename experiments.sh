#!/bin/bash

dir=`pwd`/
script_path="$( cd "$(dirname "$0")" ; pwd -P )/"

# Running experiments

# Reference genome Neisseria gonorrhoeae strain NCCP11945
# https://www.ncbi.nlm.nih.gov/nuccore/CP001050.1
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
# grep -E 'NCCP11945' assembly_summary_genbank.txt | cut -f 20

ref_genome="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/020/105/GCA_000020105.1_ASM2010v1/GCA_000020105.1_ASM2010v1_genomic.fna.gz"
genes="[(1263410,1266299),(1620898,1623787),(1724797,1727686),(1956488,1959377)]"

run_experiment() {
  mkdir ${dir}$1
  cd ${dir}$1
  sh ${script_path}pipeline.sh ${ref_genome} $2 $3 ${genes}
}

sample=( "MIC1-0-GCGS039" "MIC2-2-GCGS228" "MIC2-3-GCGS104" "MIC8-4-GCGS120" "Netherlands-69" "Netherlands-29" "Netherlands-28" "Netherlands-277" )
reads=( "ERR191768" "ERR222927" "ERR223610" "ERR223626" "SRR4418280" "SRR4418267" "SRR4418256" "SRR4418263" )
read_len=( 100 100 100 100 250 250 250 250 )

#sample=( "Netherlands-69" "Netherlands-29" "Netherlands-28" "Netherlands-277")
#reads=( "SRR4418280" "SRR4418267" "SRR4418256" "SRR4418263")
#read_len=( 250 250 250 250 )

for ((i=0;i<${#sample[@]};++i)); do
    run_experiment ${sample[i]} ${reads[i]} ${read_len[i]} &
done


# Wind et al. (300bp)
# Netherlands-69 SAMN05901177 SRR4418280 (1/4 mutated, < 30 days)
# Netherlands-29 SAMN05901175 SRR4418267 (3/4 mutated, no exposure)
# Netherlands-28 SAMN05901174 SRR4418256 (3/4 mutated, < 30 days)
# Netherlands-277 SAMN05901189 SRR4418263 (4/4 mutated, no exposure)
# Netherlands-52 SAMN05901176 SRR4418278 (no mutation, no exposure)
# Lee et al. (150bp)
# NZ2015-139 AUSMDU00006903 SRR5827361 (all 4 mutated)
# NZ2015-168 AUSMDU00006931 SRR5827099 (all 4 mutated)

# samtools mpileup -A --ff UNMAP -d 8000 -r CP001050.1:1263703-1263703 ./mappings/bowtie/bowtie-remu-sorted.bam