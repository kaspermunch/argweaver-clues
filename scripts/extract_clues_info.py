
import sys
import re
import os
import pandas as pd
import h5py

_, output_file_name, clues_dir, *clues_file_names = sys.argv # pylint: disable=unbalanced-tuple-unpacking

output_file = open(output_file_name, 'w')

for clues_file_name in clues_file_names:
    chrom, pos = re.search(r'([^_/]+)_(\d+).h5$', clues_file_name).groups()
    h5 = h5py.File(os.path.join(clues_dir, clues_file_name), 'r')
    log_likelihood_ratio = h5['logLikelihoodRatios'][h5.attrs['iHat'], h5.attrs['jHat']]
    selection_coef = h5.attrs['sHat']

    print(chrom, pos, log_likelihood_ratio, selection_coef, sep='\t', file=output_file)