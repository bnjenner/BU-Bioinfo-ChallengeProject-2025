## Guide to Elm's Scripts

1. **LeTRS_multi-primer.sh**: Running LeTRS iterating through primer reference files, with each reference file having its own run. References files are bed files, see "Elm/References" folder on the SCC. Data not on github. Can iterate through multiple samples.
2. **LeTRS_single_end.sh**: Running LeTRS on single end (or nanopore) data. Can iterate through multiple samples.
3. **LeTRS_single_primer.sh**: Running LeTRS using a single primer's bed file. Can iterate through multiple samples. 
4. **combine_LeTRS_output.py**: Combining all LeTRS outputs into megatables. Yields table with source file info.
5. **combine_LeTRS_output.sh**: Shell script to run combine_LeTRS_output.py.
6. **download_data.py**: Python script to download data from NCBI. _Not used after initial week, switched to sra_download.sh_
7. **download_data.sh**: Shell script to run download_data.py
8. **download_data_Prs_Inf.py**: Python script to download persistent infection data from NCBI. _Not used after initial week, switched to sra_download.sh_
9. **download_data_Prs_Inf.sh**: Shell script to run download_data_Prs_Inf.py.
10. **pull_data.py**: Python script to pull metadata about all sequence IDs found in NCBI Virus database. Final product used pull_all_meta()
11. **pull_data.sh**: Shell script to run pull_data.py.
12. **pull_read_res_LeTRS.py**: Script to figure out true/false positive rates from LeTRS runs. Requires ground truth file.
13. **pull_read_res_LeTRS.sh**: Shell script to run corresponding python script.
14. **qc_data.sh**: QC script, _Not used after initial week, switched to qc_data_incl_htstream.sh_
15. **qc_data_incl_htstream.sh**: QC script used for the project, includes htstream commands/filtering as well as multiqc
16. **qc_htstream_mouse.sh**: QC script used to QC mouse transcriptomics data
17. **sra_download.sh**: Download data from the SRA given a text file with a list of sequence IDs
