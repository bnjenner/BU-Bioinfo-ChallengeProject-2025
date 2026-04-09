import os
import subprocess
from sgRNAtor import utils

#################################################
# sgRNAs Class
class preproHTStream:

	def __init__(self):
		self.log_file     = None
		self.output_files = None
		self.stdout_file  = None
		self.stderr_file  = None

	#################################
	# Trim Adapter Sequences
	def trimadapaters(self, input_fastq, output_prefix, threads=1):

		# Output Files
		self.log_file = output_prefix + "_stats.json"
		self.stdout_file = output_prefix + ".stdout"
		self.stderr_file = output_prefix + ".stderr"
		self.output_files = []

		# HTStream command
		stats_command = []
		if len(input_fastq) == 1:
			stats_command = ["hts_Stats", 
								"-L", self.log_file, 
								"-N", "'stats'",
								"-U", input_fastq[0]]
			self.output_files.append(output_prefix + "_SE.fastq.gz")
		else:
			stats_command = ["hts_Stats",
								"-L", self.log_file, 
								"-N", "'stats'",
								"-1", input_fastq[0],
								"-2", input_fastq[1]]
			self.output_files.extend([output_prefix + "_R1.fastq.gz",
									  output_prefix + "_R2.fastq.gz"])

		# Add adapter trim command
		trim_command = ["hts_AdapterTrimmer",
						"-A", self.log_file,
						"-N", "'trim adapters'",
						"-F", "-f", output_prefix]


		try:

			# full_command = " ".join(stats_command) + " | " + " ".join(trim_command)
			# print(f"// Running HTStream: {full_command}")

			# Start hts_Stats process
			stats_proc = subprocess.Popen(
										 stats_command,
										 stdout=subprocess.PIPE,
										 stderr=subprocess.PIPE
										 )

			# Start hts_AdapterTrimmer Process
			trim_proc = subprocess.Popen(
										trim_command,
										stdin=stats_proc.stdout,
										stdout=subprocess.PIPE,
										stderr=subprocess.PIPE
										)

			# Close BWA stdout in parent to avoid hanging
			stats_proc.stdout.close()

			# Wait for samtools to finish, capture stderr
			trim_stderr = trim_proc.communicate()[1]
			stats_stderr = stats_proc.communicate()[1]

			# Check return codes
			if stats_proc.returncode != 0:
				raise RuntimeError(f"// ERROR: HTStream - hts_Stats Failed:\n{stats_stderr.decode()}")
			if trim_proc.returncode != 0:
				raise RuntimeError(f"// ERROR: HTStream - hts_AdapterTrimmer Failed:\n{trim_stderr.decode()}")

		except Exception as e:
			raise RuntimeError(f"Preprocessing failed: {str(e)}")

		print(f"// Output written to {self.output_files}")
