#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -P challenge2025
#$ -N LeTRS
#$ -t 1
#$ -pe omp 8
#$ -o logs/LeTRS
#$ -e logs/LeTRS
#$ -m bea

start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

# Set / Create direcotries
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/Elm
export cwd=${baseP}/scripts
export pool=0
export primer_path=/restricted/projectnb/challenge2025/Data/primers
declare -a data_arr=("AA0000144" "AY.103" "B.1.617.2" "beta" "gamma")
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
for sample in ${data_path}/*-r1.fq.gz
do
s=${sample##*/}
s=${s%-r1.fq.gz}
for p in ${primer_path}/*.bed
do
primer=${p%.*}
outP_primer=${outP}/${primer##*/}
[[ -d ${outP_primer} ]] || mkdir -p ${outP_primer}

call="perl /restricted/projectnb/challenge2025/software/LeTRS/LeTRS.pl -mode 'illumina' -fq ${data_path}/${s}-r1.fq.gz:${data_path}/${s}-r2.fq.gz -primer_bed ${p} -pool ${pool} -extractfasta -TRSLindependent -o ${outP_primer}"

echo $call
eval $call

done
done
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
