import os
import subprocess
import pysam
from sgRNAtor import utils

#################################################
# sgRNAs Class
class sgRNAquantify:

	def __init__(self, bam):
		self.bam = bam
		self.reads = {}
		self.sgRNA_counts = {}
		self.tss_dict = None


	#################################
	# Read in TSS ORFs from bed
	def read_TSS_bed(self, tss_bed, window=10):
		self.tss_dict = {}
		with open(tss_bed, "r") as bed:
			for line in bed:
				if line.startswith("#"):
					continue
				cols = line.strip().split("\t")
				
				orf, pos = str(cols[3]), int(cols[1]) - 1
				if pos in self.tss_dict:
					raise RuntimeError(f"// ERROR: TSS bed file has duplicate start positions.")
				
				self.tss_dict[pos] = {"ORF": orf,
									  "Window": (pos-window, pos+window+1),
									  "Counts": 0}


	#################################
	# Find template switching sites
	def assign_TSS_to_orfs(self, tss_bed=None, window=None):

		# TSS not read yet but specified bed and window
		if self.tss_dict is None and tss_bed is not None and window is not None:
			self.read_TSS_bed(tss_bed, window)

		# Assign sgRNAs to ORFs
		for pos, count in self.sgRNA_counts.items():
			for orf, info in self.tss_dict.items():
				if utils.overlap(pos, info["Window"]):
					self.tss_dict[orf]["Counts"] += count
					break


	#################################
	# Find template switching sites
	def find_template_switches(self, threads=1):
		
		# Read in Bam file
		for read in pysam.AlignmentFile(self.bam, "rb"):
			if not read.is_unmapped:

				# Determine R1 or R2
				pair = "R1" if not read.is_read2 else "R2"

				if read.query_name not in self.reads:
					self.reads[f"{read.query_name}"] = {}
				self.reads[f"{read.query_name}"][pair] = {"Pos": read.reference_start,
														  "Length": read.query_length,
														  "Leader": read.has_tag("LS")}

		# Reduce fragments to their TSS sites
		for fragment, reads in self.reads.items():
			
			# Determine if R1 or R2 was trimmed to get TSS site
			#   1-based conversion
			template_switch = None
			if "R1" in reads and reads["R1"]["Leader"]:
				template_switch = int(reads["R1"]["Pos"] + 1)
			elif "R2" in reads and reads["R2"]["Leader"]:
				template_switch = int(reads["R2"]["Pos"] + 1)

			if template_switch is not None:
				if template_switch not in self.sgRNA_counts:
					self.sgRNA_counts[template_switch] = 0
				self.sgRNA_counts[template_switch] += 1

	#################################
	# Output ORF TSV
	def write_ORF_counts(self, output_file):
		with open(output_file, "w") as fo:
			fo.write("ORF\tStart\tStop\tCounts\n")
			for pos, orf in self.tss_dict.items():
				fo.write(f"{orf["ORF"]}\t{orf["Window"][0]}\t{orf["Window"][1]}\t{orf["Counts"]}\n")


	#################################
	# Output sgRNAs TSV
	def write_sgRNA_counts(self, output_file):
		self.sgRNA_counts = dict(sorted(self.sgRNA_counts.items()))
		with open(output_file, "w") as fo:
			fo.write("Pos\tCounts\n")
			for pos, count in self.sgRNA_counts.items():
				fo.write(f"{pos}\t{count}\n")

