
import sys
from Bio import SeqIO

_, input_file_name, ancestral_sequence_file_name, output_file_name = sys.argv # pylint: disable=unbalanced-tuple-unpacking

input_file = open(input_file_name)
ancestral_sequence_file = open(ancestral_sequence_file_name)
output_file = open(output_file_name, 'w')

bases = 'ATGC'

ancestral_sequence = dict()
for record in SeqIO.parse(ancestral_sequence_file, "fasta"):
    ancestral_sequence[record.id] = str(record.seq)

header = next(input_file)
print(header.strip() + '\tANCESTRAL', file=output_file)

region_info = next(input_file)
_, chrom, start, end = region_info.split()
print(region_info.strip(), file=output_file)

for line in input_file:
    pos, data = line.split()
    pos = int(pos)

    ancestral_base = ancestral_sequence[chrom][pos-1]

    if ancestral_base in bases:
        print(line.strip() + ancestral_base, file=output_file)


