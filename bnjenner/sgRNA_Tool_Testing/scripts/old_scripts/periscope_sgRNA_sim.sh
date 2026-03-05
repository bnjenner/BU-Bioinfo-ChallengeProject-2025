#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N peri_sim
#$ -pe omp 8
#$ -o logs/classify
#$ -e logs/classify
#$ -m bea

start=`date +%s`
echo $HOSTNAME

threads=${OMP_NUM_THREADS}
sample="final_COV_agregate"
echo "SAMPLE: ${sample}"
echo "THREADS: ${threads}"

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/bnjenner
export cwd=${baseP}/scripts
export seqP=${baseP}/00-SimData
export outP=${baseP}/01-Periscope_Sim/${sample}
export refP=/restricted/projectnb/challenge2025/software/periscope/periscope/resources

[[ -d ${outP} ]] || mkdir -p ${outP}

module load miniconda/24.5.0
conda activate /restricted/projectnb/challenge2025/software/conda_envs/periscope

# Periscope Identify sgRNA
call="periscope \
        --fastq ${seqP}/${sample}.fastq \
	--sample ${sample} --output-prefix ${outP}/${sample} \
        --artic-primers V3 --resources ${refP} \
        --technology ont --threads ${threads}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
