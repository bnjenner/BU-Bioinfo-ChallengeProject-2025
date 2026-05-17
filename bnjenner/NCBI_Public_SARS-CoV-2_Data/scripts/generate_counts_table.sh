#!/bin/bash
sample_file=$1
prefix="NCBI_Public_SARS-CoV-2_Data"
baseP="/restricted/projectnb/challenge2025/bnjenner/sgRNAtor/bnjenner/NCBI_Public_SARS-CoV-2_Data"
input="${baseP}/01-sgRNAQuant"
output="${baseP}/02-sgRNACounts"

mkdir -p ${output}
mkdir -p ${output}/tmp

for sample in `cat ${sample_file}`; do \
    echo ${sample}
    cat ${input}/${sample}/${sample}_ORF_counts.txt | \
	tail -n +2 | cut -f 5 > ${output}/tmp/${sample}.count
done

echo ""
echo "Double Check Order"

ls ${output}/tmp/*.count > ${output}/tmp/tmp.out
paste ${sample_file} ${output}/tmp/tmp.out

temp_samp=`head -1 ${sample_file}`
tail -n +2 ${input}/${temp_samp}/${temp_samp}_ORF_counts.txt \
	 | cut -f1 > ${output}/tmp/geneids.txt

paste ${output}/tmp/geneids.txt ${output}/tmp/*.count > ${output}/tmp/tmp.out
cat <(cat ${sample_file} | paste -s) ${output}/tmp/tmp.out > ${output}/${prefix}_counts.txt
rm -rf ${output}/tmp

