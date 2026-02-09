#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N sgRNAID
#$ -t 1-1
#$ -o logs/sgRNAID
#$ -e logs/sgRNAID
#$ -m bea

start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

threads=${NSLOTS}
echo "THREADS: ${threads}"

sample=`sed "${SGE_TASK_ID}q;d" samples_PRJNA726840.txt`
echo "SAMPLE: ${sample}"

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/bnjenner # can also be set to "../"
export cwd=${baseP}/scripts
export seqP=${baseP}/00-RawData
export outP=${baseP}/01-sgRNAQuant/${sample}
export refP=${cwd}/References

[[ -d ${outP} ]] || mkdir -p ${outP}

conda activate /restricted/projectnb/challenge2025/sgRNAtor/sgRNAtor

# Reference Sequences
reference="${refP}/nCoV-2019.reference.fasta"
leader="${refP}/leader_seq.fasta"
leader_len=$(echo -n $(sed "2q;d" ${leader}) | wc -c)

# Input and Output Files
R1="${seqP}/${sample}_1.fastq.gz"
R2="${seqP}/${sample}_2.fastq.gz"
trimmed_R1="${outP}/${sample}_trimmed_R1.fastq.gz"
trimmed_R2="${outP}/${sample}_trimmed_R2.fastq.gz"
untrimmed_R1="${outP}/${sample}_unmatched_R1.fastq.gz"
untrimmed_R2="${outP}/${sample}_unmatched_R2.fastq.gz"
outbam="${outP}/${sample}_sgRNA_aligned.bam"


# Identify sgRNA leader sequence
#	Assumes all leader sequences are the same length
call="bbduk.sh \
        in1=${R1} in2=${R2} \
        outm=${trimmed_R1} outm2=${trimmed_R2} \
	out=${untrimmed_R1} out2=${untrimmed_R2} \
        ref=${leader} k=${leader_len} \
        hdist=0 mincovfraction=1 ordered=t \
        threads=${threads}"
echo $call
eval $call

# Check if Trim was Successful
if [[ ! -f "$trimmed_R1" || ! -f "$trimmed_R2" ]]; then
    echo "Error: $trimmed_R1 and/or $trimmed_R2 not found." >&2
    exit 1
fi

# Align sgRNA sequences
call="bbmap.sh ref=targets.fasta \
	in1=${trimmed_R1} in2=${trimmed_R2} \
	threads=${threads} out=${outbam}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
