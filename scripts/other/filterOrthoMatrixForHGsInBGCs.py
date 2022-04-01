import os
import sys
from Bio import SeqIO

gcf_listings_dir = os.path.abspath(sys.argv[1]) + '/'
orthogroups_matrix_file = sys.argv[2]

all_lt_in_bgcs = set([])
for f in os.listdir(gcf_listings_dir):
    with open(gcf_listings_dir + f) as ogldf:
        for line in ogldf:
            line = line.strip()
            ls = line.split('\t')
            with open(ls[1]) as ogbk:
                for rec in SeqIO.parse(ogbk, 'genbank'):
                    for feature in rec.features:
                        if not feature.type == 'CDS': continue
                        locus_tag = feature.qualifiers.get('locus_tag')[0]
                        all_lt_in_bgcs.add(locus_tag)

with open(orthogroups_matrix_file) as omf:
    for i, line in enumerate(omf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: 
            print(line)
        else:
            hg = ls[0]
            bgc_flag = False
            for lts in ls[1:]:
                for lt in lts.split(','):
                    lt = lt.strip()
                    if lt in all_lt_in_bgcs:
                        bgc_flag = True
                        break
                if bgc_flag:
                    break
            if bgc_flag:
                print(line)
