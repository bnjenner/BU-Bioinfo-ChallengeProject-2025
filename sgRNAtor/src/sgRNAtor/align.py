import os
import subprocess

#################################################
# sgRNAs Class
class alignBWA:

	def __init__(self, reference):
		self.reference = reference
		self.output_file = None

	#################################
	# Check Reference Index Exists
	def __index_exists(self):
		extensions = [".bwt", ".pac", ".ann", ".amb", ".sa"]
		return all(os.path.isfile(self.reference + ext) for ext in extensions)


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
