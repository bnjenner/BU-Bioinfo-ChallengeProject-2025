import os
import sys
import argparse
from sgRNAtor import prepro
from sgRNAtor import search
from sgRNAtor import align
from sgRNAtor import quantify
from sgRNAtor import utils
from sgRNAtor import stats

#################################################
# Argparser
def argparser():
	parser = argparse.ArgumentParser(description="Identification and Quantification pipeline for sgRNA. Performs leader sequence matching and trimming, alignment with BWA, and generates sgRNA counts tables.")
	parser.add_argument("fastq", help="Path to the input fastq file (R1 or SE)")
	parser.add_argument("fastq2", help="Path to optional Read 2 fastq file", nargs="?")  # optional positional
	parser.add_argument("--reference", "-R", type=str, required=True, help="Path to genome reference fasta file.")
	parser.add_argument("--leader-fasta", "-L", type=str, required=True, help="Path to leader sequence multi fasta file. All sequences should be no greater than 64 bp long.")
	parser.add_argument("--tss-bed", "-b", type=str, required=True, help="Path to sgRNA template switching sites bed file.")
	parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads to use (default: 1)")
	parser.add_argument("--min-match", "-m", type=int, default=12, help="Minimum length of substring to match (default: 12)")
	parser.add_argument("--max-edit", "-e", type=int, default=2, help="Maximum edit distance for a leader sequence match (default: 2)")
	parser.add_argument("--tss-window", "-w", type=int, default=10, help="Window size for template switching sites (+/- specified number). (default: 10)")
	parser.add_argument("--output-prefix", "-o", type=str, default="sgRNAtor_result", help="Prefix for output files.")
	args = parser.parse_args()

	# Check Input Files
	if not os.path.isfile(args.fastq):
		raise RuntimeError(f"// ERROR: Fastq ({args.fastq}) does not exist")
	if args.fastq2 is not None and not os.path.isfile(args.fastq):
		raise RuntimeError(f"// ERROR: Fastq Read 2 ({args.fastq2}) does not exist")

	# Check Reference Files
	if not os.path.isfile(args.reference):
		raise RuntimeError(f"// ERROR: Fasta ({args.reference}) does not exist")
	if not os.path.isfile(args.leader_fasta):
		raise RuntimeError(f"// ERROR: Fasta ({args.leader_fasta}) does not exist")
	if not os.path.isfile(args.tss_bed):
		raise RuntimeError(f"// ERROR: Bed ({args.tss_bed}) does not exist")

	# Check Parameters
	if args.min_match < 0:
		raise RuntimeError(f"// ERROR: Please use a valid minimum substring match length.")
	if args.max_edit < 0:
		raise RuntimeError(f"// ERROR: Please use a valid maximum edit distance.")
	if args.threads < 0:
		raise RuntimeError(f"// ERROR: Please use a valid number of threads.")
	if args.tss_window < 0:
		raise RuntimeError(f"// ERROR: Please use a valid window size.")

	return args
	

#################################################
# Main
def main():

	args = argparser()

	# Specify Input and Output files
	fastq_files = [args.fastq, args.fastq2]
	orfs_tsv = f"{args.output_prefix}_ORF_counts.txt"
	sgrnas_tsv = f"{args.output_prefix}_sgRNA_counts.txt"
	assignment_tsv = f"{args.output_prefix}_read_assignments.txt"
	summary_tsv = f"{args.output_prefix}_summary.txt"
	
	# Specify PE
	is_PairedEnd = False
	if fastq_files[1] is not None:
		is_PairedEnd = True
	else:
		 fastq_files = fastq_files[:-1]


	# Initialize Stat Collector
	print(f"// sgRNAtor")
	summary = stats.sgRNAstats(sample = args.output_prefix)

	# Create sgRNAsearch Object
	print("// Trimming Sequencing Adapters")
	hts = prepro.preproHTStream()
	hts.trimadapaters(input_fastq = fastq_files,
					  output_prefix = args.output_prefix,
					  threads=1)

	# Create sgRNAsearch Object
	print("// Initializing sgRNAsearch Object")
	sgRNAs = search.sgRNAsearch(fastq_files = hts.output_files,
								leader = args.leader_fasta,
								PE = is_PairedEnd)

	# Find leader sequence
	print("// Beginning sgRNA search")
	sgRNAs.find_sgRNAs(output_prefix = args.output_prefix,
					   min_match = args.min_match,
					   max_edit = args.max_edit,
					   threads = args.threads)
	# Add Stats
	summary.library_size = sgRNAs.library_size
	summary.trs_found    = sgRNAs.matches


	# Align sgRNA (leader-trimmed) sequences
	print("// Beginning BWA Alignment (sgRNA reads)")
	bwa = align.alignBWA(args.reference)
	bwa.align(input_fastq = sgRNAs.output_files,
			  output_prefix = args.output_prefix,
			  threads = args.threads)


	# Identify sgRNA Positions and ORFs
	print("// Beginning sgRNA Quantification")
	quant = quantify.sgRNAquantify(bam = bwa.output_file)
	quant.find_template_switches()
	quant.assign_TSS_to_orfs(tss_bed = args.tss_bed, window = args.tss_window)

	# Add Stats
	summary.aligned_fragments = quant.stat_counts["aligned_fragments"]
	summary.canonical         = quant.stat_counts["canonical"]
	summary.noncanonical      = quant.stat_counts["noncanonical"]
	summary.tss_dict          = quant.tss_dict
	summary.sgRNA_counts      = quant.sgRNA_counts
	summary.unassigned        = quant.unassigned


	# Write output Files
	summary.write_ORF_counts(output_file = orfs_tsv)
	summary.write_sgRNA_counts(output_file = sgrnas_tsv)
	summary.write_read_assignments(output_file = assignment_tsv)
	summary.write_summary(output_file = summary_tsv)
	print(f"// sgRNAtor Pipeline Complete.")


#################################################
if __name__ == "__main__":
	main()