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

		'''
		sgRNA Locations
		genomic = 67
		      S = 21553
		   orf3 = 25382
		      E = 26237
		      M = 26470
		   orf6 = 27041
		   orf7 = 27385
		   orf8 = 27885
		      N = 28257
		  orf9b = 28280
		     N* = 28878).
		'''

	#################################
	# Align Sequences
	def find_template_switches(self, read_length, threads=1):
		
		bamfile = pysam.AlignmentFile(self.bam, "rb")
	
		for read in bamfile:
			if not read.is_unmapped:

				pair = "R1"
				if read.is_read2:
					pair = "R2"

				if read.query_name not in self.reads:
					self.reads[f"{read.query_name}"] = {}
				self.reads[f"{read.query_name}"][pair] = {"Pos": read.reference_start,
														  "Length": read.query_length}


		for fragment, reads in self.reads.items():
			template_switch = None
			if "R1" in reads and reads["R1"]["Length"] != read_length:
				template_switch = int(reads["R1"]["Pos"])
			elif "R2" in reads and reads["R2"]["Length"] != read_length:
				template_switch = int(reads["R2"]["Pos"])

			if template_switch is not None:
				if template_switch not in self.sgRNA_counts:
					self.sgRNA_counts[template_switch] = 0
				self.sgRNA_counts[template_switch] += 1

		output = []
		for k,v in self.sgRNA_counts.items():
			output.append((k,v))

		
		sorted_output = sorted(output, key=lambda output: output[0])
		for s in sorted_output:
			print(s)