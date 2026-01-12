from pyfaidx import Fasta
from Bio.Seq import Seq

def load_reference_genome(fasta_file):
	"""Loads the reference genome from a FASTA file."""
	return Fasta(fasta_file)

def reverse_complement(seq):
	return str(Seq(seq).reverse_complement())

def extract_flanking_sequence(row, reference, flank_size):
	"""
	Extracts a flanking sequence of specified size around a mutation. The total sequence length will be (2 * flank_size + 1).
	"""
	chrom = row['CHROM']
	start = row['POS_up'] - 1
	end = row['POS_down']

	try:
	# Extract sequence from reference genome
		return reference[chrom][start:end].seq
	except KeyError:
		return None
