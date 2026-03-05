#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N periscope_multi
##$ -pe omp 8
#$ -o logs/classify
#$ -e logs/classify
#$ -m bea

start=`date +%s`
echo $HOSTNAME

threads=${OMP_NUM_THREADS}

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/bnjenner
export cwd=${baseP}/scripts
export seqP=${baseP}/00-SimData
export refP=/restricted/projectnb/challenge2025/software/periscope_multifasta/periscope_multi/resources
export tmpP=${cwd}/tmp
[[ -d ${tmpP} ]] || mkdir -p ${tmpP}

module load miniconda/24.5.0
conda activate /restricted/projectnb/challenge2025/software/conda_envs/periscope_multifasta/

for sample in `cat ART_illumina_samples.txt`;
do
	# Create Output Dir
	export outP=${baseP}/01-PeriscopeMulti_SimART/${sample}
	[[ -d ${outP} ]] || mkdir -p ${outP}

	R1="${seqP}/${sample}_R1.fq"
	R2="${seqP}/${sample}_R2.fq"

	# Periscope Multi Identify sgRNA
	call="periscope_multi \
        	--fastq ${R1} ${R2} \
        	--sample ${sample} --output-prefix ${outP}/${sample} \
        	--artic-primers V3 --resources ${refP} --gff ${refP}/covid.gff \
        	--technology illumina --threads  ${threads}"
	echo $call
	eval $call

done

end=`date +%s`
runtime=$((end-start))
echo $runtime

