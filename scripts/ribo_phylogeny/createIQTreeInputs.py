import os
import sys
from collections import defaultdict
from Bio import SeqIO

scc_algn_dir = os.path.abspath(sys.argv[1]) + '/'
strainlist_file = sys.argv[2]

strains = set([x.strip() for x in open(strainlist_file).readlines()])

gene_data = defaultdict(lambda: '')
start = 1
part_id = 1
for f in os.listdir(scc_algn_dir):
    if f.endswith('.fna'):
        gene_length = None
        visited = set([])
        with open(scc_algn_dir + f) as of:	
            for rec in SeqIO.parse(of, 'fasta'):
                genome = rec.id
                if not genome in strains: continue
                gene_data[genome] += str(rec.seq)
                gene_length = len(str(rec.seq)); visited.add(genome)
        for s in strains:
            if not s in visited:
                gene_data[s] += '-'*gene_length
        stop = start + gene_length - 1
        print('CODON, part%d = %d-%d' % (part_id, start, stop))
        start += gene_length 
        part_id += 1

outf = open('concat_msa.fa', 'w')
for genome in gene_data:
    outf.write('>' + genome + '\n' + gene_data[genome] + '\n')
outf.close()
