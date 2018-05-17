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

# NZ2015-139 AUSMDU00006903 (all 4 mutated)
mkdir ${dir}lee-C2611T-1
cd ${dir}lee-C2611T-1
sh ${script_path}pipeline.sh ${ref_genome} SRR5827361 130 ${genes}

# NZ2015-168 AUSMDU00006931 (all 4 mutated)
mkdir ${dir}lee-C2611T-2
cd ${dir}lee-C2611T-2
sh ${script_path}pipeline.sh ${ref_genome} SRR5827099 130 ${genes}

# Netherlands-69 SAMN05901177 (1/4 mutated, < 30 days)
mkdir ${dir}wind-C2611T-1
cd ${dir}wind-C2611T-1
sh ${script_path}pipeline.sh ${ref_genome} SRR4418280 200 ${genes}

# Netherlands-29 SAMN05901175 (3/4 mutated, no exposure)
mkdir ${dir}wind-C2611T-2
cd ${dir}wind-C2611T-2
sh ${script_path}pipeline.sh ${ref_genome} SRR4418267 200 ${genes}

# Netherlands-28 SAMN05901174 (3/4 mutated, < 30 days)
mkdir ${dir}wind-C2611T-3
cd ${dir}wind-C2611T-3
sh ${script_path}pipeline.sh ${ref_genome} SRR4418256 200 ${genes}

# Netherlands-277 SAMN05901189 (4/4 mutated, no exposure)
mkdir ${dir}wind-C2611T-4
cd ${dir}wind-C2611T-4
sh ${script_path}pipeline.sh ${ref_genome} SRR4418263 200 ${genes}

# Netherlands-277 SAMN05901176 (no mutation, no exposure)
mkdir ${dir}wind-C2611T-5
cd ${dir}wind-C2611T-5
sh ${script_path}pipeline.sh ${ref_genome} SRR4418278 200 ${genes}