
from pathlib import Path
import os, sys
from subprocess import Popen, PIPE
from Bio import SeqIO
import pandas as pd

_, hdf_file_name = sys.argv

if os.path.exists(hdf_file_name):
    os.remove(hdf_file_name)

data_dir = '/project/simons/faststorage/data/1000Genomes/'

pop_info_file_name = os.path.join(data_dir, 'metainfo/pop_names.tsv')

# read population identifiers
populations = []
with open(pop_info_file_name) as f:
    for line in f:
        pop, _ = line.split(maxsplit=1)
        populations.append(pop)

chromosomes = list(map(str, range(1, 23))) + ['X']
for chrom in chromosomes:

    print(chrom)

    # name of ancestral sequence file
    ancestral_file_name = os.path.join(data_dir, 'ancestral_alignmens/human_ancestor_GRCh37_e59/human_ancestor_{}.fa'.format(chrom))

    for record in SeqIO.parse(ancestral_file_name, "fasta"):
        ancestral_seq = record.seq

    for pop in populations:
        print('   ', pop)

        with open(f'steps/freq_data/{pop}_{chrom}.txt.frq') as f:

            derived_freq_records = []
            for line in f:

                if line.startswith('CHROM'):
                    continue
                _, pos, _, _, base1_info, base2_info = line.split()
                pos = int(pos)
                base1, freq1 = base1_info.split(':')
                freq1 = float(freq1)
                base2, freq2 = base2_info.split(':')
                freq2 = float(freq2)

                ancestral = ancestral_seq[pos-1]
                if ancestral not in 'ATGC':
                    # only high confidence ancestral calls
                    continue

                if not (base1 == ancestral or base2 == ancestral):
                    # only single derivced sites
                    continue

                if base1 != ancestral_seq[pos-1]:
                    derived_freq_records.append((pos, base1, freq1))
                else:
                    derived_freq_records.append((pos, base2, freq2))

        df = pd.DataFrame.from_records(derived_freq_records, columns=['pos', 'base', 'freq'])
        # print(df)
        df.to_hdf(hdf_file_name, key='{}/chr{}'.format(pop, chrom), mode='a', format='table', data_columns=['pos', 'freq'], complevel=9,complib='blosc')

    #     break
    # break