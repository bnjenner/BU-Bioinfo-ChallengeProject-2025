#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N periscope
#$ -t 1-3
#$ -o logs/classify
#$ -e logs/classify
#$ -m bea

start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

threads=${NSLOTS}
echo "THREADS: ${threads}"

sample=`sed "${SGE_TASK_ID}q;d" samples_SimArt.txt`
echo "SAMPLE: ${sample}"

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/bnjenner
export cwd=${baseP}/scripts
export seqP=${baseP}/00-SimData
export outP=${baseP}/01-Periscope_SynthIllumina/${sample}
export refP=/restricted/projectnb/challenge2025/software/periscope/periscope/resources

[[ -d ${outP} ]] || mkdir -p ${outP}

module load miniconda/24.5.0
conda activate /restricted/projectnb/challenge2025/software/conda_envs/periscope

R1="${seqP}/${sample}_R1.fq"
R2="${seqP}/${sample}_R2.fq"

ardenix_ART_illumina_100fc_R1.fq
# Periscope Identify sgRNA
call="periscope \
        --fastq ${R1} ${R2} \
	--sample ${sample} --output-prefix ${outP}/${sample} \
        --artic-primers V3 --resources ${refP} \
        --technology illumina --threads ${threads}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
