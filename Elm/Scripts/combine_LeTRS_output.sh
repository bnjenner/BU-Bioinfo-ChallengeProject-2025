#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -P challenge2025
#$ -N LeTRS_combine
#$ -t 1
#$ -pe omp 8
#$ -o logs/LeTRS_combine
#$ -e logs/LeTRS_combine
#$ -m bea

module load python3/3.13.8

python3 /restricted/projectnb/challenge2025/sgRNAtor/Elm/Scripts/combine_LeTRS_output.py