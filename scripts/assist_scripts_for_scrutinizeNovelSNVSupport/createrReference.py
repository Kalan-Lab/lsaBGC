import os
import sys
from Bio import SeqIO
assembly_listing_file = sys.argv[1] # first column sample name ; second column assembly path ; tab-separated

for line in open(assembly_listing_file):
    line = line.strip()
    ls = line.split('\t')
    with open(ls[1]) as of:
        for rec in SeqIO.parse(of, 'fasta'):
            print('>' + ls[0] + '|' + rec.description + '\n' + str(rec.seq))
