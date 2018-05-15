#!/bin/bash

dir=`pwd`/
script_path="$( cd "$(dirname "$0")" ; pwd -P )/"

# Running experiments

# Reference genome Neisseria gonorrhoeae strain NCCP11945
# https://www.ncbi.nlm.nih.gov/nuccore/CP001050.1
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
# grep -E 'NCCP11945' assembly_summary_genbank.txt | cut -f 20

ref_genome="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/020/105/GCA_000020105.1_ASM2010v1/GCA_000020105.1_ASM2010v1_genomic.fna.gz"
genes="[(1263410, 1266299), (1620898, 1623787), (1724797, 1727686), (1956488, 1959377)]"

mkdir ${dir}lee-C2611T-1
cd ${dir}lee-C2611T-1
sh ${script_path}pipeline.sh ${ref_genome} SRR5827361 130 ${genes}
