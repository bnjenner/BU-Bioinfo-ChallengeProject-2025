#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -P challenge2025
#$ -N LeTRS_Reads
#$ -t 1
#$ -pe omp 8
#$ -o logs/LeTRS_Reads.out
#$ -e logs/LeTRS_Reads.error
#$ -m bea

module load python3

python3 /restricted/projectnb/challenge2025/sgRNAtor/Elm/Scripts/pull_read_res_LeTRS.py -r /restricted/projectnb/challenge2025/markerte/sgRNAtor/Elm/LeTRS_Output/sgenerate_simdata_v1/ -g /restricted/projectnb/challenge2025/Data/sgenerate_simdata_v1/ -o /restricted/projectnb/challenge2025/markerte/sgRNAtor/Elm/LeTRS_Analysis/sim_read_res.csv --ground-truth-extension "_proportion.txt" --fasta-extension ".fa" 
