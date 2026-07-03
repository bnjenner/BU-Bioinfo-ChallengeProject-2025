import gzip
import traceback
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from contextlib import ExitStack
from concurrent.futures import ProcessPoolExecutor, as_completed
from sgRNAtor import utils


##################################################################
# Myers (1999) bit-vector approximate search.
# Finds the 5'-most end position in `text` where the FULL `pattern` aligns with
# Levenshtein distance <= max_edit. Free gaps before/after the pattern in the
# text (fitting alignment) and full indel support -- unlike the substitution-only
# bitap. Single O(len(text)) pass while pattern length <= machine word (leaders
# are <= 64 bp). Returns (end_pos, distance) or None.
def _myers_search(pattern, text, max_edit):
	m = len(pattern)
	if m == 0:
		return None
	Peq = {}
	for i, c in enumerate(pattern):
		Peq[c] = Peq.get(c, 0) | (1 << i)
	mask = (1 << m) - 1
	high = 1 << (m - 1)
	Pv = mask
	Mv = 0
	score = m
	for j in range(len(text)):
		Eq = Peq.get(text[j], 0)
		Xv = Eq | Mv
		Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq
		Ph = Mv | (~(Xh | Pv) & mask)
		Mh = Pv & Xh
		if Ph & high:
			score += 1
		elif Mh & high:
			score -= 1
		Ph = (Ph << 1) & mask
		Mh = (Mh << 1) & mask
		Pv = (Mh | (~(Xv | Ph) & mask)) & mask
		Mv = Ph & Xv
		if score <= max_edit:
			return (j, score)
	return None


##################################################################
# Partial 5' overlap: a SUFFIX of `leader` (retaining >= min_match bases) aligned
# to a PREFIX of `text` (anchored at text position 0), Levenshtein <= max_edit,
# with a free text 3' end and a free leader 5' prefix (down to min_match retained
# bases). This is the read-starts-inside-the-leader case that a plain infix search
# cannot represent. Returns (end_pos, distance) for the 5'-most end, or None.
def _suffix_prefix_overlap(leader, text, min_match, max_edit):
	m = len(leader)
	if m == 0 or min_match > m:
		return None
	free = m - min_match                       # leader-prefix positions skippable for free
	n = min(len(text), m + max_edit)
	# prev holds column j=0: D[i][0] = max(0, i - free)
	prev = [max(0, i - free) for i in range(m + 1)]
	for j in range(1, n + 1):
		cur = [0] * (m + 1)
		cur[0] = j                             # D[0][j] = j (text anchored at start)
		cj = text[j - 1]
		for i in range(1, m + 1):
			v = prev[i - 1] + (0 if leader[i - 1] == cj else 1)
			d = prev[i] + 1
			if d < v:
				v = d
			ins = cur[i - 1] + 1
			if ins < v:
				v = ins
			cur[i] = v
		prev = cur
		if prev[m] <= max_edit:
			return (j - 1, prev[m])
	return None


##################################################################
# sgRNAs Class
class sgRNAsearch:

	def __init__(self, fastq_files, leader, PE=False):
		
		# Inputs
		self.fastq_files = fastq_files
		self.leader = utils.read_fasta(leader)
		self.PE = PE

		# Bitap Masks
		self.pattern_mask = {}
		self.initialized_mask = {}

		# Library Stats
		self.matches = 0
		self.library_size = 0

		# Output
		self.output_files = None

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
	# Build leader sequence bitmasks
	def __build_bitmask(self, min_match, max_edit):
		'''
		Builds bitmask of each leader sequence in fasta file. 
		'''
		for header, seq in self.leader.items():

			# Enforce leader seq length
			if len(seq) > 64:
				raise RuntimeError(f"ERROR: Leader Sequence {header} is longer than 64 nucleotides.")

			# Initialize all possible bases (not just those in anchor)
			self.pattern_mask[seq] = {c: ~0 for c in 'ACGTN'}
			for i, c in enumerate(seq):
				self.pattern_mask[seq][c] &= ~(1 << i)

			# Build initial R state (single int; expanded to list per search call)
			self.initialized_mask[seq] = ~0
			for j in range(0, len(seq) - min_match + 1):
				self.initialized_mask[seq] &= ~(1 << j)


	#################################
	# Bitap partial 5' overlap search
	def __bitap_partial_overlap(self, read, lead, min_match, max_edit):
		"""
		Bitap (Shift-Or) search for `lead` overlapping the 5' end of `read`.
		Detects overlaps of length [min_match, len(lead)] with up to max_edit
		substitution mismatches.

		Returns on match: {"Match": True, "overlap_length": int, "mismatches": int,
		                    "anchor_start": int, "read_position": int}
		Returns on failure: {"Match": False}
		"""
		lead = str(lead)
		read = str(read)
		m = len(lead)

		pattern_mask = self.pattern_mask[lead]
		R = [self.initialized_mask[lead]] * (max_edit + 1)


		for i, c in enumerate(read):
			old_R = R[:]
			mask = pattern_mask.get(c, ~0)

			R[0] = (old_R[0] | mask) << 1
			for j in range(1, max_edit + 1):
				R[j] = ((old_R[j] | mask) << 1) & (old_R[j - 1] << 1)

			for j in range(max_edit + 1):
				if 0 == (R[j] & (1 << m)):
					overlap_len = i + 1
					if overlap_len >= min_match:
						return {
							"Match": True,
							"overlap_length": overlap_len,
							"mismatches": j,
							"anchor_start": m - overlap_len,
							"read_position": i,
						}

		return {"Match": False}


	#################################
	# Myers bit-vector leader search (indel-aware, with partial 5' overlap)
	def __myers_leader_search(self, read, lead, min_match, max_edit):
		"""
		Indel-aware leader search. Detects both:
		  - the FULL leader anywhere in the read (free 5' genomic context), and
		  - a partial 5' overlap: a leader 3' suffix (>= min_match) at the read
		    start, for reads/fragments that begin inside the leader.

		Prefers the 5'-most end position; ties go to the full-leader hit. Returns
		the same contract as __bitap_partial_overlap: overlap_length is the number
		of 5' bases to trim (through the end of the leader), leaving the body.
		"""
		read = str(read)
		lead = str(lead)
		m = len(lead)

		full = _myers_search(lead, read, max_edit)                       # full leader, free 5' context
		part = _suffix_prefix_overlap(lead, read, min_match, max_edit)   # partial 3' suffix at read start

		cands = []
		if full is not None:
			cands.append((full[0], 0, full[1]))     # (end, prefer=0 (full), dist)
		if part is not None:
			cands.append((part[0], 1, part[1]))     # (end, prefer=1 (partial), dist)
		if not cands:
			return {"Match": False}

		cands.sort(key=lambda x: (x[0], x[1]))
		end, _prefer, dist = cands[0]
		overlap_length = end + 1
		return {
			"Match": True,
			"overlap_length": overlap_length,
			"mismatches": dist,
			"anchor_start": m - overlap_length,
			"read_position": end,
		}


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

					# Revcomp for Forward & R2 or Reverse & R1 (SE: revcomp R1 on the Reverse pass)
					rev = 0
					if ((i+1)%2 == 0 and strand == "Forward") or ((i+1)%2 == 1 and strand == "Reverse"):
						_read = utils.revcomp(_read)
						_qual = _qual[::-1]
						rev = 1

					# result = self.__bitap_partial_overlap(_read, seq, min_match, max_edit)
					result = self.__myers_leader_search(_read, seq, min_match, max_edit)

					if result["Match"]:
						trim_pos = result["overlap_length"]
						_read = _read[trim_pos:]
						_qual = _qual[trim_pos:]
						_desc = f"{_record["id"]} ls:i:{rev}"
						results[strand]["sgRNA_found"] = True

					# Undo Revcomp for Forward & R2 or Reverse & R1 (restore original orientation)
					if ((i+1)%2 == 0 and strand == "Forward") or ((i+1)%2 == 1 and strand == "Reverse"):
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

		original_records = [
			SeqRecord(
				Seq(str(r["seq"])),
				id=r["id"],
				description=r["id"],
				letter_annotations={"phred_quality": list(r["qual"])}
			)
			for r in records
		]
		return {"sgRNA_found": False, "new_records": original_records}


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
	def find_sgRNAs(self, output_prefix, threads=1, min_match=8, max_edit=0, chunk_size=10000):
		'''
		This one's a bit of a beast but essentially it handles SE and PE reads without needing separte
		functions, processes them in chunks in a multithreaded fashion, and writes the output in a
		way that preserves the original order. It tries to get the benefits of streaming data
		(not storing all data into memory for reading/writing files) while also not being bottlenecked
		by gzip or single threads. 
		'''

		print(f"// Identifying sgRNA Reads in {self.fastq_files}")
		
		# Create GZIP out file handles
		out_handles = []
		for i in range(len(self.fastq_files)):

			if self.output_files is None:
				self.output_files = []

			self.output_files.append(f"{output_prefix}_sgRNA_R{i+1}.fastq.gz")
			out_handles.append(gzip.open(self.output_files[i], "wt"))


		# Build Bitmask for Bitap
		self.__build_bitmask(min_match, max_edit)


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
			out_handles[i].close()
		print(f"// Output written to {self.output_files}")


		print(f"// sgRNAs found: {self.matches}")