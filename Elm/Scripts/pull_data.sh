#!/bin/bash -l
#$ -P challenge2025
#$ -l h_rt=99:59:00
#$ -N pull_data
#$ -o logs/pull_meta
#$ -e logs/pull_meta
#$ -m bea

module load miniconda/24.5.0

conda activate /restricted/projectnb/challenge2025/software/conda_envs/pull_data_entrez

python /restricted/projectnb/challenge2025/markerte/sgRNAtor/Elm/Scripts/pull_data.py