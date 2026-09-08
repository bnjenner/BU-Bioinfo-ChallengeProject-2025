#!/bin/bash -l
#$ -P challenge2025          # Specify the SCC project name you want to use
#$ -l h_rt=99:59:00
#$ -N download_data          # Give job a name
#$ -o logs/download
#$ -e logs/download
#$ -m bea

start=`date +%s`
echo $HOSTNAME

#sra_id="PRJNA741723"
sra_id="PRJNA699398"
outpath="/restricted/projectnb/challenge2025/Data"

mkdir -p ${outpath}
cd ${outpath}

# Load SRA Toolkit
module load sratoolkit/3.0.10

# Download SRA Data
prefetch ${sra_id}
for sample in `ls .`;
do
        echo ${sample}
        fasterq-dump "${sample}"
	gzip ${sample}*.fastq
done

end=`date +%s`
runtime=$((end-start))
echo $runtime

