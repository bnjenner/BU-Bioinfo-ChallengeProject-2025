import os
import gzip
import editdistance

#################################################
# Read Fasta file
def read_fasta(fasta):
	seq = {}
	with open(fasta, "r") as file:
		curr_id = ""
		for line in file:
			if line.startswith('>'):
				curr_id = line[1:].strip()
				seq[curr_id] = ""
			else:
				seq[curr_id] += line.strip()
	return seq

#################################################
# Check if file is gzipped
def gzip_handler(file):
    try:
        with gzip.open(file, 'rb') as f:
            f.read(1) 
        return gzip.open(file, "rt")
    except (gzip.BadGzipFile):
        return open(file, "r")
    except (EOFError, zlib.error, OSError):
        raise RuntimeError(f"// ERROR: Error checking gzip status on {file}")

 #################################################
# Check if file is gzipped
def files_exist(files):
	return all([os.path.isfile(file) for file in files])

#################################################
# Reverse Compliment sequence
def revcomp(seq: str):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return "".join(complement.get(base, base) for base in reversed(seq))

#################################################
# Checks bounds overlap
def overlap(pos: int, window: tuple):
	if pos >= window[1] or pos < window[0]:
		return False
	return True
