import os
import pandas as pd
from sys import argv
import numpy as np
from utils import loadsave, seqmanip
from mutcontextSNV import calculate_flanking_positions

IDsignature_groups = {
"MMR_deficiency_signatures": ["ID7"],
"HR_deficiency_signatures":["ID6"],
"Tobacco_signatures":["ID3"],
"UV_signatures":["ID13"],
"AA_signatures":["ID23"],
"Colibactin_signatures":["ID18"],
}

def add_wt_sequencesINDEL(df, reference, flank_size, indeltype):
	"""
	Adds both `Ref_Flanking_Sequence` and `Mut_Flanking_Sequence` columns to the DataFrame.
	"""
	if indeltype=='INS':
		df["Flankseq"] = df.apply(lambda row: seqmanip.extract_flanking_sequence(row, reference, flank_size)if row["REF"] in ["C", "T"] else seqmanip.reverse_complement(seqmanip.extract_flanking_sequence(row, reference, flank_size)) ,axis=1)
	else:
		df["Flankseq"] = df.apply(lambda row: seqmanip.extract_flanking_sequence(row, reference, flank_size)if row["ALT"] in ["C", "T"] else seqmanip.reverse_complement(seqmanip.extract_flanking_sequence(row, reference, flank_size)) ,axis=1)
	df['Flankseq'] = df['Flankseq'].str.upper()

def sigfileManipINDEL(df):
	df = df.rename(columns={'Chr':'CHROM', 'Pos':'POS'})
	df['CHROM'] = 'chr' + df['CHROM'].astype(str)
	cols = list(df.filter(like='ID').columns)
	max_col_name = df[cols].idxmax(axis=1)
	df['maxID'] = max_col_name
	signature_to_group = {sig: group for group, sigs in IDsignature_groups.items() for sig in sigs}
	df["maxID_collapse"] = df["maxID"].map(signature_to_group).fillna("Unknown")
	return df

def splittingINDEL(df):
	df['REF_len'] = df['REF'].str.len()
	df['ALT_len'] = df['ALT'].str.len()

	INS = df[(df.REF_len==1)&(df.ALT_len>1)]
	DEL = df[(df.REF_len>1)&(df.ALT_len==1)]

	return [INS, DEL]

def processINDEL(df, ref_fasta, chr_list, flank_size, sig,kmer, indeltype, input_indel):
	"""
	Main function to process mutations: add flanking positions, sequences, and trinucleotide contexts.
	"""
	reference = seqmanip.load_reference_genome(ref_fasta)

	df = df[df['CHROM'].isin(chr_list)]
	df["POS"] = pd.to_numeric(df["POS"])
	df = calculate_flanking_positions(df, flank_size)

	add_wt_sequencesINDEL(df, reference, flank_size, indeltype)

	merged = pd.merge(df,sig, on=['CHROM', 'POS'], how='inner')
	merged = merged.merge(kmer.rename(lambda c: f"{c}_FLANK" if c != "kmer" else c, axis=1),how="left",left_on="Flankseq",right_on="kmer").drop(columns=["kmer"])
	merged['VAR_Type'] = indeltype
	merged['MutationID'] = merged['CHROM']+'_'+merged['POS'].astype(str) +'_'+merged['MutationType']
	merged['Sample Names'] = input_indel
	return merged

if __name__ == "__main__":
	input_indel = str(argv[1])
	input_sig = str(argv[2])
	input_kmerfile = str(argv[3])
	ref_fasta = str(argv[4])
	flank_size = int(argv[5])
	chr_file = str(argv[6])
	dir_out = str(argv[7])
	dir_bed = str(argv[8])

	indel = loadsave.load_vcf(input_indel)
	INS, DEL = splittingINDEL(indel)

	sig = loadsave.load_tsvHeaders(input_sig)
	sig = sigfileManipINDEL(sig)

	kmer = loadsave.load_tsvHeaders(input_kmerfile)

	chr_list = loadsave.load_txtlist(chr_file)
	mergedINS = processINDEL(INS, ref_fasta, chr_list,  flank_size, sig,kmer, 'INS', input_indel)
	mergedDEL = processINDEL(DEL, ref_fasta, chr_list,  flank_size, sig,kmer, 'DEL', input_indel)

	bedINS = seqmanip.vcf2bed(mergedINS, 'INS')
	bedDEL = seqmanip.vcf2bed(mergedDEL, 'DEL')

	merged = pd.concat([mergedINS, mergedDEL], axis=0, ignore_index=True)
	counting = merged.groupby(['Sample Names', 'Flankseq', 'VAR_Type'], as_index=False).size()

	loadsave.save_tsv(merged, dir_out)
	loadsave.save_tsv(counting, 'COUNT'+dir_out)
	loadsave.save_bed(bedINS, 'INS'+dir_bed)
	loadsave.save_bed(bedDEL, 'DEL'+dir_bed)
