#!/bin/bash

art=/restricted/projectnb/challenge2025/software/ART/art_src_MountRainier_Linux/art_illumina



# for ss, HiSeqX TruSeq: does 150bp paired end reads

#lens = (8418 4589 3737 3501 2778 2586 2300 2086 1714 422)
counter = 0

#for fast in /restricted/projectnb/challenge2025/sgRNAtor/ardenix/syn_sgRNAs/fastas/*; do
$art -ss HS25 -i /restricted/projectnb/challenge2025/sgRNAtor/ardenix/syn_sgRNAs/sGENERATE_files/syn_sgRNAs.faa -o /restricted/projectnb/challenge2025/Data/apl_simdata/ardenix_ART_illumina_HS25_100fc/ -l 150 -f 100 -p -m 500 -s 10 -sam 
#counter=$((counter + 1))
#../../software/ART/art_src_MountRainier_Linux/aln2bed.pl syn_sgRNAs-paired_end_com1.aln.bed syn_sgRNAs-paired_end_com1.aln syn_sgRNAs-paired_end_com2.aln