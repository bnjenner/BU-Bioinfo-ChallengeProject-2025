import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sgRNAQuant import utils

#################################################
# sgRNAs Class
class alignBWA:

	def __init__(self, reference, threads=1):
		self.reference = reference
		self.threads = threads


	#################################
	# Check Reference Index Exists
	def __index_exists(self):
		extensions = [".bwt", ".pac", ".ann", ".amb", ".sa"]
		return all(os.path.isfile(self.reference + ext) for ext in extensions)


	#################################
	# Align Sequences
	def align(self, input_fastq, output_bam):
		
		# Check Reference Index
		if not self.__index_exists():
			raise RuntimeError(f"// ERROR: Index for {self.reference} does not exist.")

		# Account for Paired End Reads
		bwa_command = ["bwa", "mem", "-C", "-t", str(self.threads), self.reference, input_fastq[0]]
		if len(input_fastq) == 2:
			bwa_command.append(input_fastq[1])

		# Alignment and Bam Compression
		bwa = subprocess.Popen(bwa_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		samtools = subprocess.Popen(
		    ["samtools", "view", "-bS", "-"],
		    stdin=bwa.stdout,
		    stdout=open(output_bam, "wb"),
		    stderr=subprocess.PIPE
		)

		bwa.stdout.close()
		samtools_stderr = samtools.communicate()[1]
		bwa_stderr = bwa.communicate()[1]

		if bwa.returncode != 0:
			raise RuntimeError(f"// ERROR: BWA Failed:\n{bwa_stderr.decode()}")
		if samtools.returncode != 0:
			raise RuntimeError(f"// ERROR: Samtools Failed:\n{samtools_stderr.decode()}")