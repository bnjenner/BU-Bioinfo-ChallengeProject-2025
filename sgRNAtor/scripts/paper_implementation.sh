#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N sgRNAID
#$ -t 1-17
#$ -t 1
#$ -o logs/paper
#$ -e logs/paper
#$ -m bea

###############################################
#
# Running this script requires a few things:
#     - Specifying project directory
#     - Specifying path to sgRNAtor software directory
#     - Making sure all reference files are in
#	    Reference directory
#
################################################


start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

threads=${NSLOTS:-1}
echo "THREADS: ${threads}"

sample=`sed "${SGE_TASK_ID}q;d" samples.txt`
echo "SAMPLE: ${sample}"


################################################
# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/bnjenner/sgRNAtor/bnjenner/sgRNAtor_Testing # Path to project dir
export cwd=${baseP}/scripts
export seqP=${baseP}/00-RawData
export outP=${baseP}/01-sgRNAQuant/${sample}
export refP=${cwd}/References
export sgRNAtorP=/restricted/projectnb/challenge2025/bnjenner/sgRNAtor/sgRNAtor/ # Path to sgRNAtor software

[[ -d ${outP} ]] || mkdir -p ${outP}

# Activate conda env
conda activate /restricted/projectnb/challenge2025/bnjenner/sgRNAtor/sgRNAtor/build


# Reference Sequences
reference="${refP}/nCoV-2019.reference.fasta"
leader="${refP}/leader_seq.fasta"
trs_bed="${refP}/sgRNA_template_switch_sites.bed"
leader_len=$(echo -n $(sed "2q;d" ${leader}) | wc -c)

# Input and Output Files
R1="${seqP}/${sample}-r1.fq.gz"
R2="${seqP}/${sample}-r2.fq.gz"
trimmed_R1="${outP}/${sample}_trimmed_R1.fastq.gz"
trimmed_R2="${outP}/${sample}_trimmed_R2.fastq.gz"
untrimmed_R1="${outP}/${sample}_unmatched_R1.fastq.gz"
untrimmed_R2="${outP}/${sample}_unmatched_R2.fastq.gz"
outsam="${outP}/${sample}_sgRNA_aligned.sam"
outbam="${outP}/${sample}_sgRNA_aligned.bam"
outcsv="${outP}/${sample}_sgRNA_ORFs.csv"


################################################
# Identify sgRNA leader sequence
#	Assumes all leader sequences are the same length
echo "// Trimming sgRNA Leader Sequences"
call="bbduk.sh \
        in1=${R1} in2=${R2} \
        outm=${trimmed_R1} outm2=${trimmed_R2} \
	out=${untrimmed_R1} out2=${untrimmed_R2} \
        ref=${leader} \
	k=${leader_len} \
        maskmiddle=f \
        ktrim=l \
		ordered=t \
        threads=${threads}"
echo $call
eval $call

# Check if Trim was Successful
if [[ ! -f "$trimmed_R1" || ! -f "$trimmed_R2" ]]; then
    echo "Error: $trimmed_R1 and/or $trimmed_R2 not found." >&2
    exit 1
fi


################################################
# Align sgRNA sequences
echo "// Aligned sgRNA Reads"
call="bbmap.sh ref=targets.fasta \
	maxindel=100 \
	32bit=t \
	mappedonly=t \
	strandedcov=t \
	strictmaxindel=t \
	in1=${trimmed_R1} in2=${trimmed_R2} \
	threads=${threads} out=${outsam}"
echo $call
eval $call

# Convert to BAM file
samtools sort ${outsam} | samtools view -Sb > ${outbam}
rm ${outsam}


# Check if Alignment was Successful
if [[ ! -f "${outbam}" ]]; then
    echo "Error: ${outbam} not found." >&2
    exit 1
fi


################################################
# Assign sgRNAs to ORFs
echo "// Assigning sgRNAs to ORFs"
call="python3 ${sgRNAtorP}/scripts/ProcessBams.py \
			--bed ${trs_bed} --window 10 \
			--output ${outcsv} \
			${outbam}"
echo $call
eval $call


end=`date +%s`
runtime=$((end-start))
echo $runtime

