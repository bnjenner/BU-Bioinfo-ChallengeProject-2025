import sys

print("// Beginning sgRNA Quantification")
quant = quantify.sgRNAquantify(bam = sys.argv[1])
quant.find_template_switches(read_length = sys.argv[2])
quant.assign_TSS_to_orfs(tss_bed = sys.argv[3], window = sys.argv[4])

print("// Writing sgRNA Counts")
quant.write_counts(output_file = sys.argv[5])
print(f"// sgRNAtor Pipeline Complete.")

