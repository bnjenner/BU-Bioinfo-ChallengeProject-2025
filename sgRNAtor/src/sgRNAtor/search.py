import gzip
import traceback
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from contextlib import ExitStack
from concurrent.futures import ProcessPoolExecutor, as_completed
from sgRNAtor import utils


##################################################################
# sgRNAs Class
class sgRNAsearch:

	def __init__(self, fastq_files, leader, PE=False):
		self.fastq_files = fastq_files
		self.leader = utils.read_fasta(leader)
		self.matches = 0
		self.library_size = 0
		self.PE = PE

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
	# Iterate through fastq files
	def __serialize_reads(self, reads):
		'''
		Passes data instead of obscurred Biopython opbject
		'''
		serial_reads = []
		for r in reads:
			serial_reads.append({"id": r.id, "seq": r.seq, "qual": r.letter_annotations["phred_quality"]})
		return serial_reads


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
				_lead = lead[l_end-i:]
			else:
				_read = read[i-l_end:i]
				_lead = lead
			if _lead == _read:
					return {"Match": True, "Pos": i}
			elif max_edit != 0:
				dist = utils.edit_distance(_lead, _read)
				if dist <= max_edit:
					return {"Match": True, "Pos": i}
		return {"Match": False, "Pos": None}


	#################################
	# Find sgRNA auxillary function
	def __find_sgRNAs(self, records, min_match=8, max_edit=0, PE=None):

		# Set PE Flag
		if PE is None:
			PE = self.PE

		# Iterate through available leader sequences
		for lead_id, seq in self.leader.items():

			results = {"Forward": {"sgRNA_found": False, "new_records": []},
					   "Reverse": {"sgRNA_found": False, "new_records": []}}

			# Libraries are unstranded, try both oreintations
			for strand in results.keys():
				
				# Iterate through R1 and R2 (or just R1 for SE)
				for i in range(len(records)):
					_record = records[i]
					_id =  _record["id"]
					_read = _record["seq"]
					_qual = _record["qual"]
					_desc = _record["id"]

					# Revcomp if PE and Forward & R2 or Reverse & R1
					rev = 0
					if PE and (((i+1)%2 == 0 and strand == "Forward") or ((i+1)%2 == 1 and strand == "Reverse")):
						_read = utils.revcomp(_read)
						_qual = _qual[::-1]
						rev = 1

					result = self.__seq_search(seq, _read, min_match, max_edit)

					if result["Match"]:
						trim_pos = result["Pos"]
						_read = _read[trim_pos:]
						_qual = _qual[trim_pos:]
						_desc = f"{_record["id"]} LS:i:{rev}"
						results[strand]["sgRNA_found"] = True

					# Undo Revcomp if PE and Forward & R2 or Reverse & R1
					if PE and (((i+1)%2 == 0 and strand == "Forward") or ((i+1)%2 == 1 and strand == "Reverse")):
						_read = utils.revcomp(_read)
						_qual = _qual[::-1]

					# New Seq Record
					new_record = SeqRecord(
					    Seq(_read),
					    id=_id,
					    description=_desc,
					    letter_annotations={"phred_quality": _qual}
					)
					results[strand]["new_records"].append(new_record)

				if results[strand]["sgRNA_found"]:
					return results[strand]

		return {"sgRNA_found": False, "new_records":	 []}


	#################################
	# Find sgRNA from chunks of reads
	def find_sgRNAs_chunk(self, records, min_match, max_edit):
		results = []
		try:
			for r in records:
				results.append(self.__find_sgRNAs(r, min_match, max_edit))
			return {"ok": True, "results": results}
		except Exception:
			return {"ok": False, "traceback": traceback.format_exc()}


	#################################
	# Find sgRNA main function
	def find_sgRNAs(self, output_files, threads=1, min_match=8, max_edit=0, chunk_size=10000):
		'''
		This one's a bit of a beast but essentially it handles SE and PE reads without needing separte
		functions, processes them in chunks in a multithreaded fashion, and writes the output in a
		way that preserves the original order. It tries to get the benefits of streaming data
		(not storing all data into memory for reading/writing files) while also not being bottlenecked
		by gzip or single threads. 
		'''

		print(f"// Trimming FASTQ Files {self.fastq_files}")

		# Check Inputs Match Outputs
		if len(output_files) != len(self.fastq_files):
			msg = (f"// ERROR: Number of input and output files do not match." +
				   f"//     Input:  {self.fastq_files}" +
				   f"//     Output: {output_files}")
			raise RuntimeError(msg)

		
		# Create GZIP out file handles	
		out_handles = []
		for i in range(len(self.fastq_files)):
			out_handles.append(gzip.open(output_files[i], "wt"))


		next_seq = 0         # Allows for seq iteration
		seq_num_base = 0     # Keeps track of seq num for chunks
		chunk = []           # Input batch buffer
		futures = {}         # Stores results of thread execution
		out_buffer = {}      # Stores results in order

		# Multithreaded sgRNA Identification and Trimming
		with ProcessPoolExecutor(max_workers=threads) as pool:
			
			# Iterate over FASTQ records
			for record in self.__iterate_reads(self.fastq_files):
				chunk.append(self.__serialize_reads(record))
				self.library_size += 1

				# If desired chunk size reached, execute sgRNA search
				if len(chunk) >= chunk_size:
					fut = pool.submit(self.find_sgRNAs_chunk, chunk, min_match, max_edit)
					futures[fut] = (seq_num_base, len(chunk))
					seq_num_base += len(chunk)
					chunk = []

				# Collect finished jobs opportunistically
				done = [f for f in futures if f.done()]
				for f in done:
					base, size = futures.pop(f)
					res = f.result()

					# Check results
					if not res["ok"]:
						raise RuntimeError("Worker crashed:\n" + res["traceback"])

					# Store results in output buffer for ordered output
					for i, result_dict in enumerate(res["results"]):
						out_buffer[base + i] = result_dict

					# Write sequentially (within and across chunks)
					while next_seq in out_buffer:
						result_dict = out_buffer.pop(next_seq) # Removes from memory as written
						if result_dict["sgRNA_found"]:
							for i, rec in enumerate(result_dict["new_records"]):
								SeqIO.write(rec, out_handles[i], "fastq")
							self.matches += 1
						next_seq += 1

			# Submit remaining records
			if chunk:
				fut = pool.submit(self.find_sgRNAs_chunk, chunk, min_match, max_edit)
				futures[fut] = (seq_num_base, len(chunk))
				seq_num_base += len(chunk)

			# Wait for final completion of thread pool
			for f in as_completed(futures):
				base, size = futures[f]
				res = f.result()

				# Check results
				if not res["ok"]:
					raise RuntimeError("Worker crashed:\n" + res["traceback"])

				# Store results in output buffer for ordered output
				for i, result_dict in enumerate(res["results"]):
					out_buffer[base + i] = result_dict

				# Write sequentially (within and across chunks)
				while next_seq in out_buffer:
					result_dict = out_buffer.pop(next_seq) # Removes from memory as written
					if result_dict["sgRNA_found"]:
						for i, rec in enumerate(result_dict["new_records"]):
							SeqIO.write(rec, out_handles[i], "fastq")
						self.matches += 1
					next_seq += 1

		# Close all opened GZIP output files
		for i in range(len(out_handles)):
			out_handles[i].close
			print(f"// Output written to {output_files[i]}")


		print(f"// sgRNAs found: {self.matches}")