import os
import pandas as pd
from sys import argv
import numpy as np
from utils import loadsave, seqmanip

def calculate_flanking_positions(df, flank_len):
	"""
	Calculate downstream and upstream positions for each mutation and add them as new columns to the DataFrame.
	"""
	df['POS_down'] = df['POS'] + flank_len
	df['POS_up'] = df['POS'] - flank_len
	return df

def add_wt_sequences(df, reference, flank_size):
	"""
	Adds both `Ref_Flanking_Sequence` and `Mut_Flanking_Sequence` columns to the DataFrame.
	"""
	df["WT_flankseq"] = df.apply(lambda row: seqmanip.extract_flanking_sequence(row, reference, flank_size)if row["REF"] in ["C", "T"] else seqmanip.reverse_complement(seqmanip.extract_flanking_sequence(row, reference, flank_size)) ,axis=1)
	df['WT_flankseq'] = df['WT_flankseq'].str.upper()

def add_mut_sequences(df, reference, flank_size):
	"""
	Adds a `Mut_flankseq` column to the DataFrame.
	Uses `WT_flankseq` and adjusts for the mutation, maintaining reverse complement orientation if necessary.
	"""
	def get_mut_seq(row):
		wt_seq = row['WT_flankseq']
		if row['REF'] in ['C', 'T']:
			# On the same strand, just replace the middle base
			mid = flank_size  # assuming the mutation is in the middle
			return wt_seq[:mid] + row['ALT'] + wt_seq[mid+1:]
		else:
			# If REF is A or G, sequence was reversed before, so we also reverse complement ALT
			mid = flank_size
			mut_seq = wt_seq[:mid] + row['ALT'] + wt_seq[mid+1:]
			return seqmanip.reverse_complement(mut_seq)

	df['MUT_flankseq'] = df.apply(get_mut_seq, axis=1)
	df['MUT_flankseq'] = df['MUT_flankseq'].str.upper()

def processSNV(df, ref_fasta,  flank_size):
	"""
	Main function to process mutations: add flanking positions, sequences, and trinucleotide contexts.
	"""
	reference = seqmanip.load_reference_genome(ref_fasta)
	df["POS"] = pd.to_numeric(df["POS"])
	df = calculate_flanking_positions(df, flank_size)
	add_wt_sequences(df, reference, flank_size)
	add_mut_sequences(df, reference, flank_size)
	print(df.head())
	return


if __name__ == "__main__":
	input_snv = str(argv[1])
	input_sig = str(argv[2])
	input_kmerfile = str(argv[3])
	ref_fasta = str(argv[4])
	flank_size = int(argv[5])

	snv = loadsave.load_vcf(input_snv)
	#sig = loadsave.load_tsvHeaders(input_sig)
	processSNV(snv, ref_fasta,  flank_size)

	#kmer = loadsave.load_tsvHeaders(input_kmerfile)
	#print(kmer.head())
