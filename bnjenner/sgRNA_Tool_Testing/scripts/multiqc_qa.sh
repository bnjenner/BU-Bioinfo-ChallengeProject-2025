#!/bin/bash -l

start=`date +%s`

# Set / Create Directories
export baseP=/restricted/projectnb/challenge2025/sgRNAtor/bnjenner
export cwd=${baseP}/scripts
export outP=${baseP}/01-HTStream_prepro
export qcP=${baseP}/02-MultiQC


# HTStream Preprocessing
module load miniconda/24.5.0
conda activate /restricted/projectnb/challenge2025/software/conda_envs/HTStream_v1.4.1

call="multiqc -i multiqc_PRJNA726840_and_PRJNA741723 \
       -o ${qcP} -f ${outP}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
