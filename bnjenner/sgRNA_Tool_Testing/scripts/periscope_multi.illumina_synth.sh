#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N periscope_multi
#$ -o logs/classify
#$ -e logs/classify
#$ -m bea

start=`date +%s`
echo $HOSTNAME

threads=${NSLOTS}
echo "THREADS: ${threads}"

sample_file="samples_SimArt.txt"
echo "SAMPLES: ${sample_file}"

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/bnjenner
export cwd=${baseP}/scripts
export seqP=${baseP}/00-SimData
export refP=/restricted/projectnb/challenge2025/software/periscope_multifasta/periscope_multi/resources

module load miniconda/24.5.0
conda activate /restricted/projectnb/challenge2025/software/conda_envs/periscope_multifasta/

for sample in `cat ${sample_file}`;
do

	# Create Output Dir
	export outP=${baseP}/01-PeriscopeMulti_SynthIllumina/${sample}
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

