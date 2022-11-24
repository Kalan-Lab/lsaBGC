import os
import sys
from Bio import SeqIO

valid_chars = set(['A', 'C', 'G', 'T'])
fasta_file = sys.argv[1]
strains_file = sys.argv[2]
ambig_cutoff = float(sys.argv[3])

strains = set([x.strip() for x in open(strains_file).readlines()])

fasta_data = []
for rec in SeqIO.parse(fasta_file, 'fasta'):
        if rec.id in strains:
                fasta_rec = [rec.id] + list(str(rec.seq).upper())
                fasta_data.append(fasta_rec)

fasta_data_tr = []
for i, ls in enumerate(zip(*fasta_data)):
        if i == 0:
                fasta_data_tr.append(ls)
        else:
                chars = set(ls).intersection(valid_chars)
                n_count = len([x for x in ls if not x in valid_chars])
                if float(n_count)/len(ls) < ambig_cutoff and len(chars) > 1:
                        new_ls = ''
                        for b in ls:
                            if not  b in valid_chars:
                                b = '-'
                            new_ls += b
                        fasta_data_tr.append(new_ls)

for rec in zip(*fasta_data_tr):
        print('>' + rec[0] + '\n' + ''.join(rec[1:]))
