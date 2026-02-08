#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -P challenge2025
#$ -N sgRNAID
#$ -t 1-15
#$ -o logs/sgRNAID
#$ -e logs/sgRNAID
#$ -m bea

start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

threads=${NSLOTS}
echo "THREADS: ${threads}"

sample=`sed "${SGE_TASK_ID}q;d" samples.txt`
echo "SAMPLE: ${sample}"

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/bnjenner # can also be set to "../"
export cwd=${baseP}/scripts
export seqP=${baseP}/00-RawData
export outP=${baseP}/01-sgRNAQuant/${sample}
export refP=${cwd}/References

[[ -d ${outP} ]] || mkdir -p ${outP}

conda activate /restricted/projectnb/challenge2025/sgRNAtor/sgRNAQuant

R1="${seqP}/${sample}.fastq.gz"

# Periscope Identify sgRNA
call="periscope \
        --fastq ${R1} \
        --sample ${sample} --output-prefix ${outP}/${sample} \
        --artic-primers V3 --resources ${refP} \
        --technology ont --threads ${threads}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime



call="bbduk.sh \
	in1="$R1" in2="$R2" \
	outm1="${OUT_PREFIX}_R1.fastq.gz" \
  	outm2="${OUT_PREFIX}_R2.fastq.gz" \
  	outu1="${OUT_PREFIX}_unmatched_R1.fastq.gz" \
  	outu2="${OUT_PREFIX}_unmatched_R2.fastq.gz" \
  	ref="$REF" \
  	k="$SEQ_LEN" \
  	hdist=0 \
  	mincovfraction=1 \
  	ordered=t \
  	threads=4"
