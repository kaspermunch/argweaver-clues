
import sys
import re
import pandas as pd
import h5py

_, clues_file_name = sys.argv

chrom, pos = re.search(r'([^_/]+)_(\d+).h5$', clues_file_name).groups()

h5 = h5py.File(clues_file_name, 'r')
log_likelihood_ratio = h5['logLikelihoodRatios'][h5.attrs['iHat'], h5.attrs['jHat']]
selection_coef = h5.attrs['sHat']

print(chrom, pos, log_likelihood_ratio, selection_coef, sep='\t')