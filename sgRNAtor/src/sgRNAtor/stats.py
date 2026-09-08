import os
import subprocess
import pysam
from sgRNAtor import utils

#################################################
# sgRNAs Class
class sgRNAstats:

	def __init__(self, sample):
		self.sample = sample

		# Counts
		self.library_size = 0
		self.trs_found = 0
		self.aligned_fragments = 0
		self.canonical = 0
		self.noncanonical = 0

		# Reads and Assignments
		self.tss_dict = {}
		self.sgRNA_counts = {}
		self.unassigned = []

	#################################
	# Output Summary TSV
	def write_summary(self, output_file):
		with open(output_file, "w") as fo:
			fo.write(f"Sample\t{self.sample}\n")
			fo.write(f"Library_Size\t{self.library_size}\n")
			fo.write(f"TRS_Found\t{self.trs_found}\n")
			fo.write(f"Aligned\t{self.aligned_fragments}\n")
			fo.write(f"Canonical\t{self.canonical}\n")
			fo.write(f"Noncanonical\t{self.noncanonical}\n")
		print(f"// Output written to {output_file}")


	#################################
	# Output ORF TSV
	def write_ORF_counts(self, output_file):
		with open(output_file, "w") as fo:
			fo.write("ORF\tStart\tStop\tCounts\tCPMs\n")
			for pos, orf in self.tss_dict.items():
				line = (f"{orf["ORF"]}\t" +
						f"{orf["Window"][0]}\t" + 
						f"{orf["Window"][1]}\t" +
						f"{orf["Counts"]}\t" +
						f"{round(orf["Counts"]/(self.library_size/1e6),3)}\n")
				fo.write(line)
		print(f"// Output written to {output_file}")


	#################################
	# Output sgRNAs TSV
	def write_sgRNA_counts(self, output_file):
		self.sgRNA_counts = dict(sorted(self.sgRNA_counts.items()))
		with open(output_file, "w") as fo:
			fo.write("Pos\tAssigned\tCounts\tCPMs\n")
			for pos, info in self.sgRNA_counts.items():
				line = (f"{pos}\t" + 
						f"{info["Assigned"]}\t" +
						f"{info["Counts"]}\t" +
						f"{round(info["Counts"]/(self.library_size/1e6),3)}\n")
				fo.write(line)
		print(f"// Output written to {output_file}")


	#################################
	# Output Read Assignents
	def write_read_assignments(self, output_file):
		with open(output_file, "w") as fo:
			fo.write("ORF\tReads\n")
			for pos, orf in self.tss_dict.items():
				line = (f"{orf["ORF"]}\t" +
						f"{",".join(orf["Reads"])}\n")
				fo.write(line)
			line = (f"{"unassigned"}\t" +
				    f"{",".join(self.unassigned)}\n")
			fo.write(line)
		print(f"// Output written to {output_file}")

