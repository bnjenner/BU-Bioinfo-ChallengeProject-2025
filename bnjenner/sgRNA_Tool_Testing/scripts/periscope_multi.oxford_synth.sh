#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N periscope_multi
#$ -o logs/classify
#$ -e logs/classify
#$ -m bea

start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

threads=${NSLOTS}
echo "THREADS: ${threads}"

sample="final_COV_agregate"
echo "SAMPLE: ${sample}"

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/bnjenner
export cwd=${baseP}/scripts
export seqP=${baseP}/00-SimData
export outP=${baseP}/01-PeriscopeMulti_SynthOxford/${sample}
export refP=/restricted/projectnb/challenge2025/software/periscope_multifasta/periscope_multi/resources
export tmpP=${cwd}/tmp

[[ -d ${outP} ]] || mkdir -p ${outP}
[[ -d ${tmpP} ]] || mkdir -p ${tmpP}

module load miniconda/24.5.0
conda activate /restricted/projectnb/challenge2025/software/conda_envs/periscope_multifasta/

# Periscope Multi Identify sgRNA
call="periscope_multi \
        --fastq ${seqP}/${sample}.fastq \
	--sample ${sample} --output-prefix ${outP}/${sample} \
        --artic-primers V3 --resources ${refP} --gff ${refP}/covid.gff \
        --technology ont --threads  ${threads}"
echo $call
eval $call


end=`date +%s`
runtime=$((end-start))
echo $runtime
