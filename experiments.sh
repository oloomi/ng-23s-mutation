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

# FA1090
#ref_genome="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/845/GCF_000006845.1_ASM684v1/GCF_000006845.1_ASM684v1_genomic.fna.gz"
#genes="[(1116249,1119158),(1258352,1261255),(1649928,1652830),(1873080,1875982)]"

run_experiment() {
  echo "Running experiments for sample {$1} ..."
  mkdir ${dir}$1
  cd ${dir}$1
  sh ${script_path}pipeline.sh ${ref_genome} $2 $3 ${genes}
  sync
  echo "Experiments for sample {$1} completed!"
}

sample=( "MIC1-0-GCGS039" "MIC2-2-GCGS228" "MIC2-3-GCGS104" "MIC8-4-GCGS120" )
reads=( "ERR191768" "ERR222927" "ERR223610" "ERR223626" )
read_len=( 100 100 100 100 )

#sample=( "MIC1-0-GCGS039" "MIC2-2-GCGS228" "MIC2-3-GCGS104" "MIC8-4-GCGS120" "Netherlands-69" "Netherlands-29" "Netherlands-28" "Netherlands-277" )
#reads=( "ERR191768" "ERR222927" "ERR223610" "ERR223626" "SRR4418280" "SRR4418267" "SRR4418256" "SRR4418263" )
#read_len=( 100 100 100 100 250 250 250 250 )

for ((i=0;i<${#sample[@]};++i)); do
    run_experiment ${sample[i]} ${reads[i]} ${read_len[i]} &
done

python3 ${script_path}mutation_count.py "${sample[@]}"
