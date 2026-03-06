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
		self.ambiguous = 0


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
	def assign_TSS_to_orfs(self, tss_bed=None, window=10):

		# TSS not read yet but specified bed and window
		if self.tss_dict is None and tss_bed is not None and window is not None:
			self.read_TSS_bed(tss_bed, window)

		# Assign sgRNAs to ORFs
		for pos, counts in self.sgRNA_counts.items():
			for orf, info in self.tss_dict.items():
				if utils.overlap(pos, info["Window"]):
					self.tss_dict[orf]["Counts"] += counts["Counts"]
					self.sgRNA_counts[pos]["Assigned"] = orf
					break


	#################################
	# Find template switching sites
	def find_template_switches(self, threads=1, has_tag=False):
		'''
		Parses aligned reads and identifies which read was trimmed and also where the 
		junction site occured. This identifies all junction sites and generates counts
		for them. This will be used later for sgRNA ORF assignment.
		'''
		
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

			template_switch = 0
			for r, attr in reads.items():

				# 1-based conversion 
				tss = int(attr["Pos"] + 1)

				# Grab 3' most TSS site
				if attr["Leader"] and tss > template_switch:
					template_switch = tss

			if template_switch != 0:
				if template_switch not in self.sgRNA_counts:
					self.sgRNA_counts[template_switch] = {"Counts": 0, "Assigned": None}
				self.sgRNA_counts[template_switch]["Counts"] += 1


	#################################
	# Output ORF TSV
	def write_ORF_counts(self, output_file):
		with open(output_file, "w") as fo:
			fo.write("ORF\tStart\tStop\tCounts\n")
			for pos, orf in self.tss_dict.items():
				fo.write(f"{orf["ORF"]}\t{orf["Window"][0]}\t{orf["Window"][1]}\t{orf["Counts"]}\n")
		print(f"// Output written to {output_file}")


	#################################
	# Output sgRNAs TSV
	def write_sgRNA_counts(self, output_file):
		self.sgRNA_counts = dict(sorted(self.sgRNA_counts.items()))
		with open(output_file, "w") as fo:
			fo.write("Pos\tCounts\tAssigned\n")
			for pos, info in self.sgRNA_counts.items():
				fo.write(f"{pos}\t{info["Counts"]}\t{info["Assigned"]}\n")
		print(f"// Output written to {output_file}")

