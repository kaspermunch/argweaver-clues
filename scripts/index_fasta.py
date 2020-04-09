

import sys
from pyfaidx import Fasta

_, fasta_file_name = sys.argv

chromosomes = Fasta(fasta_file_name)
