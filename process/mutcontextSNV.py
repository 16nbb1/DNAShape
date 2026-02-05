import os
import pandas as pd
from sys import argv
import numpy as np
from utils import loadsave, seqmanip

SBSsignature_groups = {
"MMR_deficiency_signatures": ["SBS6", "SBS14","SBS15","SBS20", "SBS21","SBS26","SBS44"],
"POL_deficiency_signatures":["SBS10a", "SBS10b","SBS10c", "SBS10d","SBS28"],
"HR_deficiency_signatures":["SBS3"],
"BER_deficiency_signatures":["SBS30", "SBS36"],
"Treatment_signatures":["SBS11", "SBS25", "SBS31", "SBS32","SBS35", "SBS86", "SBS87", "SBS90", "SBS99"],
"Immunosuppressants_signatures":["SBS32"],
"APOBEC_signatures":["SBS2", "SBS13"],
"Tobacco_signatures":["SBS4", "SBS29", "SBS92"],
"UV_signatures":["SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"],
"AA_signatures":["SBS22a", "SBS22b"],
"Colibactin_signatures":["SBS88"],
"Artifact_signatures":["SBS27", "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57", "SBS58", "SBS59", "SBS60", "SBS95"],
"Lymphoid_signatures":["SBS9", "SBS84", "SBS85"]
}

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


def sigfileManipSNV(df):
	df = df.rename(columns={'Chr':'CHROM', 'Pos':'POS'})
	df['CHROM'] = 'chr' + df['CHROM'].astype(str)
	cols = list(df.filter(like='SBS').columns)
	max_col_name = df[cols].idxmax(axis=1)
	df['maxSBS'] = max_col_name
	signature_to_group = {sig: group for group, sigs in SBSsignature_groups.items() for sig in sigs}
	df["maxSBS_collapse"] = df["maxSBS"].map(signature_to_group).fillna("Unknown")
	return df

def processSNV(df, ref_fasta,  flank_size, sig, kmer, input_snv):
	"""
	Main function to process mutations: add flanking positions, sequences, and trinucleotide contexts.
	"""
	reference = seqmanip.load_reference_genome(ref_fasta)
	df["POS"] = pd.to_numeric(df["POS"])
	df = calculate_flanking_positions(df, flank_size)
	add_wt_sequences(df, reference, flank_size)
	add_mut_sequences(df, reference, flank_size)
	merged = pd.merge(df,sig, on=['CHROM', 'POS'], how='inner')
	merged = merged.merge(kmer.rename(lambda c: f"{c}_WT" if c != "kmer" else c, axis=1),how="left",left_on="WT_flankseq",right_on="kmer").drop(columns=["kmer"])
	merged = merged.merge(kmer.rename(lambda c: f"{c}_MUT" if c != "kmer" else c, axis=1),how="left",left_on="MUT_flankseq",right_on="kmer").drop(columns=["kmer"])

	features = ['Buckle_7', 'HelT_7', 'MGW_7', 'Opening_7', 'ProT_7', 'Rise_7', 'Roll_7', 'Shear_7', 'Shift_7', 'Slide_7', 'Stagger_7', 'Stretch_7', 'Tilt_7', 'Bend_net', 'Bend_max']
	for i in features:
		merged[i+'_eucdist'] = merged.apply(lambda row: np.linalg.norm(seqmanip.safe_to_array(row[i+'_WT']) - seqmanip.safe_to_array(row[i+'_MUT'])),axis=1)
		merged[i+'_paireddiff'] = merged.apply(lambda row: (seqmanip.safe_to_array(row[i+'_WT']) - seqmanip.safe_to_array(row[i+'_MUT'])).tolist(),axis=1)
		merged[i+'_sumdiff'] = merged.apply(lambda row: np.sum(row[i+'_paireddiff']),axis=1)

	merged['VAR_Type'] = 'SNV'
	merged['MutationID'] = merged['CHROM']+'_'+merged['POS'].astype(str) +'_'+merged['MutationType']
	merged['Sample Names'] = input_snv
	return merged


if __name__ == "__main__":
	input_snv = str(argv[1])
	input_sig = str(argv[2])
	input_kmerfile = str(argv[3])
	ref_fasta = str(argv[4])
	flank_size = int(argv[5])
	chr_file = str(argv[6])
	dir_out = str(argv[7])
	dir_bed = str(argv[8])

	snv = loadsave.load_vcf(input_snv)

	sig = loadsave.load_tsvHeaders(input_sig)
	sig = sigfileManipSNV(sig)

	kmer = loadsave.load_tsvHeaders(input_kmerfile)
	if len(snv)>0:
		merged = processSNV(snv, ref_fasta,  flank_size, sig, kmer, input_snv)
		bed = seqmanip.vcf2bed(merged, 'SNV')

		loadsave.save_tsv(merged, dir_out)
		loadsave.save_bed(bed, 'SNV'+dir_bed)
