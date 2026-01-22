import os
import sys
from sys import argv
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator import install as genInstall

#genInstall.install('GRCh38')
#genInstall.install('GRCh37')
#genInstall.install('c_elegans')

vcf_directory = str(argv[1])
output_directory = str(argv[2])
ref = str(argv[3])
vartype = str(argv[4])

if (ref == 'GRCh38.99'):
	if (vartype=='SNV'):
		Analyze.cosmic_fit(samples=vcf_directory,output=output_directory,input_type="vcf",context_type="96",genome_build="GRCh38",exome=False,export_probabilities_per_mutation=True,cosmic_version=3.4)
	elif (vartype=='INDEL'):
		Analyze.cosmic_fit(samples=vcf_directory,output=output_directory,input_type="vcf",context_type="ID",genome_build="GRCh38",exome=False,export_probabilities_per_mutation=True,collapse_to_SBS96 = False,cosmic_version=3.4)

#elif ref == 'hg19':
#	Analyze.cosmic_fit(samples=vcf_directory,output=output_directory,input_type="vcf",context_type="96",genome_build="GRCh37",exome=False,export_probabilities_per_mutation=True,cosmic_version=3.4)
#elif ref == 'ce11':
#	Analyze.cosmic_fit(samples=vcf_directory,output=output_directory,input_type="vcf",context_type="96",genome_build="c_elegans",exome=False,export_probabilities_per_mutation=True,cosmic_version=3.4)
