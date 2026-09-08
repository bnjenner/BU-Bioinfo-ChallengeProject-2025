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
		self.tss_dict = {}
		self.unassigned = []
		self.stat_counts = {
			"aligned_fragments": 0,
			"canonical": 0,
			"noncanonical": 0
		}


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
									  "Counts": 0,
									  "Reads": []}


	#################################
	# Find template switching sites
	def assign_TSS_to_orfs(self, tss_bed=None, window=10):

		# TSS not read yet but specified bed and window
		if not self.tss_dict and tss_bed is not None and window is not None:
			self.read_TSS_bed(tss_bed, window)

		# Assign sgRNAs to ORFs
		for pos, counts in self.sgRNA_counts.items():
			_assigned = False
			for orf, info in self.tss_dict.items():

				# Assign Counts and Read IDs
				if utils.overlap(pos, info["Window"]):
					self.tss_dict[orf]["Counts"] += counts["Counts"]
					self.tss_dict[orf]["Reads"].extend(counts["Reads"])
					self.sgRNA_counts[pos]["Assigned"] = orf
					self.stat_counts["canonical"] += counts["Counts"]
					_assigned = True
					break

			if not _assigned:
				self.unassigned.extend(counts["Reads"])
				self.stat_counts["noncanonical"] += counts["Counts"]


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
			if not read.is_unmapped and not read.is_supplementary:

				# Determine R1 or R2
				pair = "R1" if not read.is_read2 else "R2"

				if read.query_name not in self.reads:
					self.reads[f"{read.query_name}"] = {}
				self.reads[f"{read.query_name}"][pair] = {"Pos": read.reference_start,
														  "Length": read.query_length,
														  "Leader": read.has_tag("ls")}

		# Reduce fragments to their TSS sites
		for fragment, reads in self.reads.items():	

			# Add to stats
			self.stat_counts["aligned_fragments"] += 1

			template_switch = 0			
			for r, attr in reads.items():

				# 1-based conversion 
				tss = int(attr["Pos"] + 1)

				# Grab 3' most TSS site
				if attr["Leader"] and tss > template_switch:
					template_switch = tss

			if template_switch != 0:
				if template_switch not in self.sgRNA_counts:
					self.sgRNA_counts[template_switch] = {"Counts": 0, "Assigned": None, "Reads": []}
				self.sgRNA_counts[template_switch]["Counts"] += 1
				self.sgRNA_counts[template_switch]["Reads"].append(fragment)

