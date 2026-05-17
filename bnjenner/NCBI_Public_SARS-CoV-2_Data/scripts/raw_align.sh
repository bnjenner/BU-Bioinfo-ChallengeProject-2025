#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N sgrnaquant
#$ -t 1-75
#$ -pe omp 4
#$ -l mem_per_core=2G
#$ -o logs/raw_align
#$ -e logs/raw_align

start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

threads=4
echo "THREADS: ${threads}"

sample=`sed "${SGE_TASK_ID}q;d" samples.txt`
echo "SAMPLE: ${sample}"

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/bnjenner/sgRNAtor/bnjenner/NCBI_Public_SARS-CoV-2_Data
export cwd=${baseP}/scripts
export seqP=${baseP}/01-sgRNAQuant/${sample}
export outP=${baseP}/01-PreproAligned/${sample}
export refP=${baseP}/References

# References
reference=${refP}/nCoV-2019.reference.fasta
leader=${refP}/leader_seq.fasta
orf=${refP}/sgRNA_template_switch_sites.bed

[[ -d ${outP} ]] || mkdir -p ${outP}

module load miniconda/24.5.0
conda activate /restricted/projectnb/challenge2025/bnjenner/sgRNAtor/sgRNAtor/build

fastq_R1="${seqP}/${sample}_R1.fastq.gz"
fastq_R2="${seqP}/${sample}_R2.fastq.gz"


# Align Preprocessed Reads
call="bwa mem -t ${threads} \
	${reference} \
	${fastq_R1} ${fastq_R2} | \
     samtools sort \
	-o ${outP}/${sample}_raw.sorted.bam"
echo $call
eval $call

# Index Aligned Reads
call="samtools index ${outP}/${sample}_raw.sorted.bam"
echo $call
eval $call


end=`date +%s`
runtime=$((end-start))
echo $runtime
