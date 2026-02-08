import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from contextlib import ExitStack
from concurrent.futures import ThreadPoolExecutor, as_completed
from sgRNAQuant import utils

#################################################
# sgRNAs Class
class sgRNAsearch:

	def __init__(self, fastq_files, leader, fastq2=None):
		self.fastq_files = fastq_files
		self.leader = utils.read_fasta(leader)
		self.matches = 0


	#################################
	# Iterate through fastq files
	def __iterate_reads(self, fastq_files):
		'''
		Arbitrary handle for SE or PE reads and handles file closing after iterations
		'''
		with ExitStack() as stack:
			handles = [
				SeqIO.parse(stack.enter_context(utils.gzip_handler(f)), "fastq")
				for f in fastq_files
			]
			for records in zip(*handles):
				yield records


	#################################
	# Sequence Search utility
	def __seq_search(self, lead, read, min_match, max_edit):
		'''
		Performs full pattern search and also slides 3' end of leader 
		over 5' end of read and calcualtes edit distance.
		Returns dictionary (Match, Position)
		'''
		r_end = len(read)
		l_end = len(lead)
		for i in range(r_end, min_match-1, -1):
			if i < l_end:
				_read = read[:i]
				_lead = lead[:i]
			else:
				_read = read[i-l_end:i]
				_lead = lead
			dist = utils.edit_distance(_lead, _read)
			if dist <= max_edit:
				return {"Match": True, "Pos": i}
		return {"Match": False, "Pos": None}


	#################################
	# Find sgRNA auxillary function
	def __find_sgRNAs(self, records, min_match=8, max_edit=0):

		for lead_id, seq in self.leader.items():
			new_records = []
			sgRNA_found = False

			# Iterate through R1 and R2 (or just R1 for SE)
			for i in range(len(records)):
				_record = records[i]
				result = self.__seq_search(seq, _record.seq, min_match, max_edit)

				_id = _record.id
				_read = _record.seq
				_qual = _record.letter_annotations["phred_quality"]
				_desc = _record.id
				
				if result["Match"]:
					trim_pos = result["Pos"]
					_desc = f"{_desc} sgRNA"
					_read = _read[trim_pos:]
					_qual = _qual[trim_pos:]
					sgRNA_found = True

				# New Seq Record
				new_record = SeqRecord(
				    Seq(_read),
				    id=_id,
				    description=_desc,
				    letter_annotations={"phred_quality": _qual}
				)
				new_records.append(new_record)

			if sgRNA_found:
				break

		return {"sgRNA_found": sgRNA_found, "new_records": new_records}




	#################################
	# Find sgRNA main function
	def find_sgRNAs(self, output_files, threads=1, min_match=8, max_edit=0, chunk_size=1000):

		print(f"// Trimming FASTQ Files {self.fastq_files}")

		# Check Inputs Match Outputs
		if len(output_files) != len(self.fastq_files):
			msg = (f"// ERROR: Number of input and output files do not match." +
				   f"//     Input:  {self.fastq_files}" +
				   f"//     Output: {output_files}")
			raise RuntimeError(msg)

		
		# Handle output files		
		out_handles = []
		for i in range(len(self.fastq_files)):
			out_handles.append(gzip.open(output_files[i], "wt"))


		# Multithreaded sgRNA Identification and Trimming
		buffer = {}
		next_seq = 0
		with ThreadPoolExecutor(max_workers=threads) as pool:
			chunk = []
			seq_num_base = 0

			for record in self.__iterate_reads(self.fastq_files):

				chunk.append(record)
				if len(chunk) >= chunk_size:
					futures = {pool.submit(self.__find_sgRNAs, r, min_match, max_edit): seq_num_base + i for i, r in enumerate(chunk)}

					for fut in as_completed(futures):
						seq_num = futures[fut]
						buffer[seq_num] = fut.result()

						# write sequentially ready dicts
						while next_seq in buffer:
							result_dict = buffer.pop(next_seq)
							sgRNA_found = result_dict["sgRNA_found"]
							new_records = result_dict["new_records"]

							# Break for output if sgRNA detected
							if sgRNA_found:
								for i in range(len(new_records)):
									SeqIO.write(new_records[i], out_handles[i], "fastq")
								self.matches += 1
								next_seq += 1

			# Process any remaining reads
			if chunk:
				futures = {pool.submit(self.__find_sgRNAs, r, min_match, max_edit): seq_num_base + i for i, r in enumerate(chunk)}

				for fut in as_completed(futures):
						seq_num = futures[fut]
						buffer[seq_num] = fut.result()

						# write sequentially ready dicts
						while next_seq in buffer:
							result_dict = buffer.pop(next_seq)
							sgRNA_found = result_dict["sgRNA_found"]
							new_records = result_dict["new_records"]

							# Break for output if sgRNA detected
							if sgRNA_found:
								for i in range(len(new_records)):
									SeqIO.write(new_records[i], out_handles[i], "fastq")
								self.matches += 1
								next_seq += 1

		for i in range(len(out_handles)):
			out_handles[i].close
			print(f"// Output written to {output_files[i]}")
