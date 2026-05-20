#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -P challenge2025
#$ -N LeTRS
#$ -t 1
#$ -pe omp 8
#$ -o logs/LeTRS_SE
#$ -e logs/LeTRS_SE
#$ -m bea

start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

# Set / Create direcotries
export baseP=/restricted/projectnb/challenge2025/markerte/sgRNAtor/Elm
export cwd=${baseP}/scripts
export pool=0
export primer_path=/restricted/projectnb/challenge2025/Data/primers/artic_primers_v3.bed
declare -a data_arr=("sgenerate_simdata_v1")
# Module Loading
module load miniconda/24.5.0

# Activating conda environment
conda activate LeTRS

# iterating through datasets
for data in "${data_arr[@]}"
do

echo "DATA: ${data}"

# Set / Create Directories for each dataset
export data_path=/restricted/projectnb/challenge2025/Data/${data}
export outP=${baseP}/LeTRS_Output/${data}

[[ -d ${outP} ]] || mkdir -p ${outP}


# Running LeTRS on each data and primer combination
for sample in ${data_path}/final_COV_TAT*.fastq
do
s=${sample##*/}
s=${s%.fastq}
echo "SAMPLE: ${s}"
outP_sample=${outP}/${s}
call="perl /restricted/projectnb/challenge2025/software/LeTRS/LeTRS.pl -mode 'illumina' -fq ${data_path}/${s}.fastq -primer_bed ${primer_path} -pool ${pool} -extractfasta -TRSLindependent -o ${outP_sample}"

echo $call
eval $call

done
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
