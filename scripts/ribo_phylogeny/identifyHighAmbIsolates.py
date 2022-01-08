import os
import sys
from Bio import SeqIO

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    gaps = len([x for x in str(rec.seq) if x == '-'])
    total = len(str(rec.seq))
    
    if float(gaps)/float(total) > 0.25:
        print(rec.id)#+ '\t' + str(float(gaps)/float(total)))
