import os
import sys
from sys import argv
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator import install as genInstall

#genInstall.install('GRCh38')

vcf_directory = str(argv[1])
output_directory = str(argv[2])
ref = str(argv[3])

if ref == 'GRCh38.99':
	Analyze.cosmic_fit(samples=vcf_directory,output=output_directory,input_type="vcf",context_type="ID",genome_build="GRCh38",exome=False,export_probabilities_per_mutation=True,collapse_to_SBS96 = False,cosmic_version=3.4)
