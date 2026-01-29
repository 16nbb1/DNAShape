import os
import pandas as pd
from sys import argv
import numpy as np
from sys import argv
from utils import loadsave

def mergingContext(input_loc):
	li=[]
	for i in input_loc:
		print(i)
		df = loadsave.load_tsvHeaders(i)
		li.append(df)
	big = pd.concat(li)
	return big

if __name__ == "__main__":
	output_files = str(argv[1])
	input_files = argv[2:]
	big = mergingContext(input_files)
	loadsave.save_tsv(big, output_files)
