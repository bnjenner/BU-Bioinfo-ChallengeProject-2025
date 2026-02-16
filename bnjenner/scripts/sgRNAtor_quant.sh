#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N sgrnaquant
#$ -t 1-11
#$ -o logs/sgRNAQuant
#$ -e logs/sgRNAQuant
#$ -m bea

start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

threads=${NSLOTS}
echo "THREADS: ${threads}"

sample=`sed "${SGE_TASK_ID}q;d" samples_PRJNA726840.txt`
echo "SAMPLE: ${sample}"

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/bnjenner/sgRNAtor/bnjenner
export cwd=${baseP}/scripts
export seqP=${baseP}/00-RawData
export outP=${baseP}/01-sgRNAQuant_PRJNA726840/${sample}
export refP=${cwd}/References

# References
reference=${refP}/nCoV-2019.reference.fasta
leader=${refP}/leader_seq.fasta
orf=${refP}/sgRNA_template_switch_sites.bed

[[ -d ${outP} ]] || mkdir -p ${outP}

module load miniconda/24.5.0
conda activate /restricted/projectnb/challenge2025/bnjenner/sgRNAtor/sgRNAtor/build

fastq_R1="${seqP}/${sample}_1.fastq.gz"
fastq_R2="${seqP}/${sample}_2.fastq.gz"

# Run sgRNAQuant
call="sgRNAtor -R ${reference} -L ${leader} -b ${orf} \
        --threads ${threads} -m 10 -e 1 \
        --output-prefix ${outP}/${sample} \
        ${fastq_R1} ${fastq_R2}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
