#!/bin/bash
script="../../../software/periscope_multifasta/periscope_multi/scripts/create_multi_reference.py "
ref="../../../software/periscope_multifasta/periscope_multi/resources/nCoV-2019.reference.fasta"
gff="../../../software/periscope_multifasta/periscope_multi/resources/covid.gff"
call="python ${script} \
	--fasta ${ref} --gff ${gff}"
echo $call
eval $call
