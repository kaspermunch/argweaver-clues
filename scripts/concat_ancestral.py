
import sys, re
from Bio import SeqIO

_, *file_names = sys.argv

for file_name in file_names:
    with open(file_name) as f:
        for record in SeqIO.parse(f, "fasta"):
            chrom = re.search(r'ANCESTOR_for_chromosome:GRCh37:([^:]+):.*', record.id).group(1)
            record.id = chrom
            record.description = ''
            SeqIO.write([record], sys.stdout, "fasta")

#python ../scripts/concat_ancestral.py ~/simons/faststorage/data/1000Genomes/ancestral_alignmens/human_ancestor_GRCh37_e59/human_ancestor_*.fa > ~/simons/faststorage/data/1000Genomes/ancestral_alignmens/human_ancestor_GRCh37_e59/all_chromosomes.fa