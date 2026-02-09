
# Module Loading
module load fastqc/0.12.1
module load miniconda/24.5.0

conda activate /restricted/projectnb/challenge2025/software/conda_envs/HTStream_v1.4.1

export fasts=/restricted/projectnb/challenge2025/Data/apl_simdata/ART_sGENERATE_nc_sgRNAs/sGEN_nc_sgRNAs1.fq
export qcout=/restricted/projectnb/challenge2025/Data/
export outP=${qcout}/sGEN_nc_sgRNAs1

[[ -d ${outP} ]] || mkdir -p ${outP}

proj="sGEN_nc_sgRNAs1"

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
