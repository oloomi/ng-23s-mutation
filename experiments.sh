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
  mkdir -p ${dir}$1
  cd ${dir}$1
  sh ${script_path}pipeline.sh ${ref_genome} $2 $3 ${genes}
  echo "Experiments for sample {$1} completed!"
}

sample=( "GCGS039" "GCGS196" "GCGS102" "GCGS228" "GCGS104" "GCGS142" "GCGS038" "GCGS106" "GCGS118" "GCGS120" "GCGS128" "GCGS136" "GCGS202" "GCGS204" "GCGS224" "GCGS230" "GCGS018" "GCGS096" "GCGS098" "GCGS110" "GCGS156" "GCGS174" "GCGS176" "GCGS198" "GCGS218" "GCGS226" )
reads=( "ERR191768" "ERR222895" "ERR223608" "ERR222927" "ERR223610" "ERR223648" "ERR191767" "ERR223612" "ERR223624" "ERR223626" "ERR223634" "ERR223642" "ERR222901" "ERR222903" "ERR222923" "ERR222929" "ERR191747" "ERR191825" "ERR223604" "ERR223616" "ERR223662" "ERR223680" "ERR223682" "ERR222897" "ERR222917" "ERR222925" )
read_len=( 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 )

#sample=( "GCGS196" "GCGS102" "GCGS228" )
#reads=( "ERR222895" "ERR223608" "ERR222927" )
#read_len=( 100 100 100 )

#sample=( "MIC1-0-GCGS039" "MIC2-2-GCGS228" "MIC2-3-GCGS104" "MIC8-4-GCGS120" )
#reads=( "ERR191768" "ERR222927" "ERR223610" "ERR223626" )
#read_len=( 100 100 100 100 )

#sample=( "MIC1-0-GCGS039" "MIC2-2-GCGS228" "MIC2-3-GCGS104" "MIC8-4-GCGS120" "Netherlands-69" "Netherlands-29" "Netherlands-28" "Netherlands-277" )
#reads=( "ERR191768" "ERR222927" "ERR223610" "ERR223626" "SRR4418280" "SRR4418267" "SRR4418256" "SRR4418263" )
#read_len=( 100 100 100 100 250 250 250 250 )

trap "exit" INT TERM ERR
trap "kill 0" EXIT

for ((i=0;i<${#sample[@]};++i)); do
    run_experiment ${sample[i]} ${reads[i]} ${read_len[i]} &
done

wait
python3 ${script_path}mutation_count.py "${sample[@]}"
echo "All experiments completed!"
