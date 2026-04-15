#!/bin/bash -l
#$ -P challenge2025          # Specify the SCC project name you want to use
#$ -l h_rt=99:59:00
#$ -N download_data          # Give job a name
#$ -o logs/download
#$ -e logs/download
#$ -m bea

start=`date +%s`
echo $HOSTNAME

pulled_data="ncbi_virus_acc"
outpath="/restricted/projectnb/challenge2025/Data/${pulled_data}"

mkdir -p ${outpath}
cd ${outpath}

# Load SRA Toolkit
module load sratoolkit/3.0.10

# Download SRA data from text file
# Code modified from asadprodhan on github
while IFS= read -r accession; do
    prefetch $accession && fasterq-dump $accession --split-files
done < "/restricted/projectnb/challenge2025/markerte/sgRNAtor/Elm/Scripts/ncbi_virus_acc.txt"

# Download SRA Data
#prefetch ${sra_id}
#for sample in `ls .`;
#do
#        echo ${sample}
#        fasterq-dump "${sample}"
#	gzip ${sample}*.fastq
#done

end=`date +%s`
runtime=$((end-start))
echo $runtime

