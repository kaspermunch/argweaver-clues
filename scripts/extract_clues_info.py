
import sys
import re
import os
import pandas as pd
import h5py

_, output_file_name, steps_dir, *clues_file_names = sys.argv # pylint: disable=unbalanced-tuple-unpacking

# open output file:
output_file = open(output_file_name, 'w')

# loop over base names of 
for clues_file_name in clues_file_names:

    chrom, start, end, pop, pos, chain = re.search(r'([^_]+)_(\d+)_(\d+)_([^_]+)_(\d+)_(\d+).h5', clues_file_name).groups()
    h5_file_name = os.path.join(steps_dir, clues_file_name)

    if os.path.getsize(h5_file_name) == 0:
        print(h5_file_name, 'IS EMPTY... (probably because I could not get CLUES to run)')
        continue

    h5 = h5py.File(h5_file_name, 'r')
    log_likelihood_ratio = h5['logLikelihoodRatios'][h5.attrs['iHat'], h5.attrs['jHat']]
    selection_coef = h5.attrs['sHat']
    snp_freq = h5.attrs['snp_freq']

    print(chrom, start, end, pop, pos, chain, log_likelihood_ratio, selection_coef, snp_freq, sep='\t', file=output_file)
