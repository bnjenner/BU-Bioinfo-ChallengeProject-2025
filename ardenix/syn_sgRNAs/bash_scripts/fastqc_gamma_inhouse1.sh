
# Module Loading
module load fastqc/0.12.1
module load miniconda/24.5.0

conda activate /restricted/projectnb/challenge2025/software/conda_envs/HTStream_v1.4.1

export fasts=/restricted/projectnb/challenge2025/Data/gamma/gamma-r1.fq.gz
export qcout=/restricted/projectnb/challenge2025/Data/
export outP=${qcout}/gamma-r1.fq.gz

[[ -d ${outP} ]] || mkdir -p ${outP}

proj="gamma-r1.fq.gz"

# QC with fastqc
call="fastqc ${fasts} \
       --outdir ${outP} \
       --dir ${outP}"
echo $call
eval $call
done

# Create MultiQC Report
call="multiqc -i multiqc_${proj} \
	-o ${outP} -f ${outP}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
