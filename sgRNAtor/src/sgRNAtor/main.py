import os
import sys
import argparse
from sgRNAtor import search
from sgRNAtor import align
from sgRNAtor import quantify
from sgRNAtor import utils

#################################################
# Argparser
def argparser():
	parser = argparse.ArgumentParser(description="Identification and Quantification pipeline for sgRNA. Performs leader sequence matching and trimming, alignment with BWA, and generates sgRNA counts tables.")
	parser.add_argument("fastq", help="Path to the input fastq file (R1 or SE)")
	parser.add_argument("fastq2", help="Path to optional Read 2 fastq file", nargs="?")  # optional positional
	parser.add_argument("--reference", "-R", type=str, default=None, help="Path to optional genome reference fasta file.")
	parser.add_argument("--leader-fasta", "-L", type=str, default=None, help="Path to optional leader sequence multi fasta file.")
	parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads to use (default: 1)")
	parser.add_argument("--min-match", "-m", type=int, default=8, help="Minimum length of substring to match (default: 8)")
	parser.add_argument("--max-edit", "-e", type=int, default=0, help="Maximum edit distance for a leader sequence match (default: 0)")
	parser.add_argument("--output-prefix", "-o", type=str, default="sgRNAtor_result", help="Prefix for output files.")
	parser.add_argument("--force-overwrite", "-f", action='store_true', help="Force overwrite of intermediate files.")
	args = parser.parse_args()

	# Check Input Files
	if not os.path.isfile(args.fastq):
		raise RuntimeError(f"// ERROR: Fastq ({args.fastq}) does not exist")
	if args.fastq2 is not None and not os.path.isfile(args.fastq):
		raise RuntimeError(f"// ERROR: Fastq Read 2 ({args.fastq2}) does not exist")

	# Check Reference Files
	curr_path = os.path.dirname(os.path.abspath(__file__))
	if args.leader_fasta is None:
		args.leader_fasta = os.path.join(curr_path, "../data/leader_seq.fasta")
	elif not os.path.isfile(args.leader_fasta):
		raise RuntimeError(f"// ERROR: Fasta ({args.leader_fasta}) does not exist")

	if args.reference is None:
		args.reference = os.path.join(curr_path, "../data/nCoV-2019.reference.fasta")
	elif not os.path.isfile(args.reference):
		raise RuntimeError(f"// ERROR: Fasta ({args.reference}) does not exist")

	# Check Parameters
	if args.min_match < 0:
		raise RuntimeError(f"// ERROR: Please use a valid minimum substring match length.")

	if args.max_edit < 0:
		raise RuntimeError(f"// ERROR: Please use a valid maximum edit distance.")

	if args.max_edit < 0:
		raise RuntimeError(f"// ERROR: Please use a valid maximum edit distance.")

	return args
	

#################################################
# Main
def main():

	args = argparser()

	# Specify Input and Output files
	fastq_files = [args.fastq, args.fastq2]
	trimmed_files = []
	aligned_file = f"{args.output_prefix}_aligned_sgRNA.bam"
	for i in range(len(fastq_files)):
		if fastq_files[i] is None:
			continue
		read = f"R{i+1}"
		trimmed_files.append(f"{args.output_prefix}_trimmed_{read}.fastq.gz")

	# Specify PE
	is_PairedEnd = False
	if len(trimmed_files) == 2:
		is_PairedEnd = True

	# Create sgRNAsearch Object
	print(f"// sgRNAtor")
	print("// Initializing sgRNAsearch Object")
	sgRNAs = search.sgRNAsearch(fastq_files = fastq_files,
								leader = args.leader_fasta,
								PE = is_PairedEnd)


	# Find leader sequence
	if args.force_overwrite or not utils.files_exist(trimmed_files):
		print("// Beginning sgRNA search")
		sgRNAs.find_sgRNAs(output_files = trimmed_files,
						   threads = args.threads,
						   min_match = args.min_match,
						   max_edit = args.max_edit)
	else:
		print(f"// NOTICE: Trimmed FASTQ Files Found {trimmed_files}. Skipping sgRNA Search.")


	# Align Trimmed Sequences
	if args.force_overwrite or not utils.files_exist(aligned_file):
		print("// Beginning BWA Alignment")
		bwa = align.alignBWA(args.reference)
		bwa.align(input_fastq = trimmed_files, 
				  output_bam = aligned_file,
				  threads = args.threads)
	else:
		print(f"// NOTICE: Aligned BAM File Found {aligned_files}. Skipping alignment.")

	print("// Beginning sgRNA Quantification")
	quant = quantify.sgRNAquantify(bam = aligned_file)
	quant.find_template_switches(read_length = sgRNAs.read_length)
	print(f"// sgRNAtor Pipeline Complete.")


#################################################
if __name__ == "__main__":
	main()