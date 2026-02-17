#!/bin/bash

# Inputs
sample="SRR31567433"
fastq1="00-RawData/${sample}_1.fq.gz"
fastq2="00-RawData/${sample}_2.fq.gz"

# References
Reference="References/nCoV-2019.reference.fasta"
Leader="References/leader_seq.fasta"
TSS_Bed="References/sgRNA_template_switch_sites.bed"

# Outputs
output="01-sgRNAQuant/${sample}"
mkdir -p ${output}

call="sgRNAtor --reference ${Reference} \
               --leader-fasta ${Leader} \
               --tss-bed ${TSS_Bed} \
               --threads 12 \
               --min-match 10 \
               --max-edit 2 \
               --tss-window 10 \
               --output-prefix ${output}/${sample} \
               ${fastq1} ${fastq2}"
echo $call
eval $call
