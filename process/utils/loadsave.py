'''
Utilities useful for frequent loading and saving of data
'''

import pandas as pd
from pyfaidx import Fasta

def load_vcf(file_path):
	"""Load a CSV file into a DataFrame."""
	return pd.read_csv(file_path, sep='\t', comment='#',names=['CHROM','POS','ID','REF','ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL','TUMOR'], header=None)

def load_tsvHeaders(file_path):
	return pd.read_csv(file_path, sep='\t')

def save_tsv(df, file_path):
	"""Save DataFrame to CSV."""
	df.to_csv(file_path, index=False,  sep='\t')
	print(f"Data saved to {file_path}")
