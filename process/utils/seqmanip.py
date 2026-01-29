from pyfaidx import Fasta
from Bio.Seq import Seq
import numpy as np
import ast

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

def safe_to_array(val):
	if isinstance(val, str):
		val = ast.literal_eval(val)
	return np.atleast_1d(np.array(val, dtype=float))

def vcf2bed(df, VAR_Type):
	if (VAR_Type =='SNV'):
		df['STOP'] = df['POS']
		df = df[['CHROM','POS','STOP','REF', 'ALT', 'VAR_Type', 'MutationType', 'MutationID','WT_flankseq',  'Sample Names']]
	elif  (VAR_Type =='INS'):
		df['STOP'] = df['POS']
		df = df[['CHROM','POS','STOP','REF', 'ALT', 'VAR_Type', 'MutationType', 'MutationID', 'Flankseq', 'Sample Names']]
	else:
		df['STOP'] = df['POS']+df['REF_len']-1
		df = df[['CHROM','POS','STOP','REF', 'ALT', 'VAR_Type', 'MutationType', 'MutationID', 'Flankseq', 'Sample Names']]
	return df
