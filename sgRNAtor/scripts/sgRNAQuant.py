import argparse
from sgRNAtor import quantify
from sgRNAtor import utils

#################################################
# Argparser
def argparser():
	parser = argparse.ArgumentParser(description="Identification and Quantification pipeline for sgRNA. Performs leader sequence matching and trimming, alignment with BWA, and generates sgRNA counts tables.")
	parser.add_argument("bam", help="Aligned BAM file.")
	parser.add_argument("--tss-bed", "-b", type=str, default=None, required=True, help="Path to sgRNA template switching sites bed file.")
	parser.add_argument("--tss-window", "-w", type=int, default=10, help="Window size for template switching sites (+/- specified number). (default: 10)")
	parser.add_argument("--output-prefix", "-o", type=str, default="sgRNAtor_result", help="Prefix for output files.")
	args = parser.parse_args()

	# Check Input Files
	if not os.path.isfile(args.bam):
		raise RuntimeError(f"// ERROR: Fastq ({args.fastq}) does not exist")

	# TSS BedFile
	if args.tss_bed is None:
		args.tss_bed = os.path.join(curr_path, "../data/sgRNA_template_switch_sites.bed")
	elif not os.path.isfile(args.tss_bed):
		raise RuntimeError(f"// ERROR: Fasta ({args.tss_bed}) does not exist")

	# TSS Window
	if args.tss_window < 0:
		raise RuntimeError(f"// ERROR: Please use a valid window size.")

	return args


#################################################
# Main
def main():

	args = argparser()

	aligned_file = args.bam
	orfs_tsv = f"{args.output_prefix}_ORF_counts.txt"
	sgrnas_tsv = f"{args.output_prefix}_sgRNA_counts.txt"

	print("// Beginning sgRNA Quantification")
	quant = quantify.sgRNAquantify(bam = aligned_file)
	quant.find_template_switches()
	quant.assign_TSS_to_orfs(tss_bed = args.tss_bed, window = args.tss_window)

	print("// Writing ORF Counts")
	quant.write_ORF_counts(output_file = orfs_tsv)
	print("// Writing sgRNA Counts")
	quant.write_sgRNA_counts(output_file = sgrnas_tsv)
	print(f"// sgRNAtor Pipeline Complete.")


#################################################
if __name__ == "__main__":
	main()
