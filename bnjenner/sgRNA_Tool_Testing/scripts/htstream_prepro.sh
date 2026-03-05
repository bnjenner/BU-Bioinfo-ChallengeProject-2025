#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -P challenge2025
#$ -N htstream
#$ -t 1-11
#$ -pe omp 8
#$ -o logs/prepro
#$ -e logs/prepro
#$ -m bea

start=`date +%s`
echo $HOSTNAME
echo "My SGE_TASK_ID: " $SGE_TASK_ID

sample=`sed "${SGE_TASK_ID}q;d" samples_PRJNA726840.txt`
echo "SAMPLE: ${sample}"

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/bnjenner/
export cwd=${baseP}/scripts
export seqP=${baseP}/00-RawData
export outP=${baseP}/01-HTStream_prepro/${sample}
export qcP=${baseP}/02-MultiQC
export rrna=${baseP}/References/human_rRNA.fasta

[[ -d ${outP} ]] || mkdir -p ${outP}


# HTStream Preprocessing
module load miniconda/24.5.0
conda activate /restricted/projectnb/challenge2025/software/conda_envs/HTStream_v1.4.1
call="hts_Stats -L ${outP}/${sample}.json -N 'initial stats' \
          -1 ${seqP}/${sample}_1.fastq.gz -2 ${seqP}/${sample}_2.fastq.gz | \
      hts_SeqScreener -r -A ${outP}/${sample}.json -N 'screen phix' | \
      hts_SeqScreener -A ${outP}/${sample}.json -N 'count the number of rRNA reads' \
          -r -s ${rrna} | \
      hts_AdapterTrimmer -A ${outP}/${sample}.json -N 'trim adapters' | \
      hts_QWindowTrim -A ${outP}/${sample}.json -N 'quality trim the ends of reads' | \
      hts_NTrimmer -A ${outP}/${sample}.json -N 'remove any remaining N bases' | \
      hts_LengthFilter -A ${outP}/${sample}.json -N 'remove reads < 30 bp' \
          -n -m 30 | \
      hts_Stats -A ${outP}/${sample}.json -N 'final stats' \
          -f ${outP}/${sample}"
echo $call
eval $call


# Additional QC with fastqc
module load fastqc/0.12.1
call="fastqc ${outP}/${sample}_*fastq.gz \
       --outdir ${outP} \
       --dir ${outP}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
