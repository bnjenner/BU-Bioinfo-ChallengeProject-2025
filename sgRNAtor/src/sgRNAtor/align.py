import os
import re
import subprocess
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sgRNAtor import utils

#################################################
# sgRNAs Class
class alignBWA:

	def __init__(self, reference):
		self.reference = reference
		self.output_file = None
		self.gene_counts = {}

	#################################
	# Check Reference Index Exists
	def __index_exists(self):
		extensions = [".bwt", ".pac", ".ann", ".amb", ".sa"]
		return all(os.path.isfile(self.reference + ext) for ext in extensions)


	#################################
	# Parse gene features from GTF (gene-level only, skips subpoly products)
	def _parse_gtf_genes(self, gtf_file):
		genes = {}
		with open(gtf_file) as f:
			for line in f:
				if line.startswith('#'):
					continue
				cols = line.strip().split('\t')
				if len(cols) < 9 or cols[2] != 'gene':
					continue
				start = int(cols[3]) - 1  # convert to 0-based
				end = int(cols[4])
				gene_name = None
				for attr in cols[8].split(';'):
					m = re.match(r'\s*gene\s+"([^"]+)"', attr)
					if m:
						gene_name = m.group(1)
						break
				if gene_name:
					genes[gene_name] = (cols[0], start, end)
		return genes


	#################################
	# Assign aligned reads to ORFs using GTF gene intervals
	def assign_reads(self, gtf_file):
		print(f"// Assigning reads to ORFs using {gtf_file}")
		gene_intervals = self._parse_gtf_genes(gtf_file)

		self.gene_counts = {name: 0 for name in gene_intervals}
		self.gene_counts["unassigned"] = 0

		with pysam.AlignmentFile(self.output_file, "rb") as bam:
			for read in bam:
				if read.is_unmapped or read.is_supplementary or read.is_secondary:
					continue
				assigned = False
				for gene_name, (chrom, start, end) in gene_intervals.items():
					if read.reference_start >= start and read.reference_start < end:
						self.gene_counts[gene_name] += 1
						assigned = True
						break
				if not assigned:
					self.gene_counts["unassigned"] += 1

		print(f"// Read assignment complete")
		return self.gene_counts


	#################################
	# Align Sequences
	def align(self, input_fastq, output_prefix, threads=1, name="sgRNA"):

		# Set Output File
		self.output_file = f"{output_prefix}_aligned_{name}.bam"

		# Check Reference Index
		if not self.__index_exists():
			raise RuntimeError(f"// ERROR: Index for {self.reference} does not exist.")

		# BWA MEM command
		bwa_command = ["bwa", "mem", "-C", "-t", str(threads), self.reference, input_fastq[0]]
		if len(input_fastq) == 2:
			bwa_command.append(input_fastq[1])

		# Open BAM file for writing
		with open(self.output_file, "wb") as bam_out:
			try:
				# Start BWA process
				bwa_proc = subprocess.Popen(
										bwa_command,
										stdout=subprocess.PIPE,
										stderr=subprocess.PIPE
										)

				# Start samtools process, reading from bwa stdout
				samtools_proc = subprocess.Popen(
											["samtools", "view", "-b", "-h", "-"],
											stdin=bwa_proc.stdout,
											stdout=bam_out,
											stderr=subprocess.PIPE
											)

				# Close BWA stdout in parent to avoid hanging
				bwa_proc.stdout.close()

				# Wait for samtools to finish, capture stderr
				samtools_stderr = samtools_proc.communicate()[1]
				bwa_stderr = bwa_proc.communicate()[1]

				# Check return codes
				if bwa_proc.returncode != 0:
					raise RuntimeError(f"// ERROR: BWA Failed:\n{bwa_stderr.decode()}")
				if samtools_proc.returncode != 0:
					raise RuntimeError(f"// ERROR: Samtools Failed:\n{samtools_stderr.decode()}")

			except Exception as e:
				raise RuntimeError(f"Alignment failed: {str(e)}")

		print(f"// Output written to {self.output_file}")
