import os
import sys

itol_keep_strains = set([x.strip() for x in open(sys.argv[1]).readlines()])

for line in open(sys.argv[2]):
    line = line.strip()
    ls = [x for x in line.split('_') if x != '']
    alt_name = ' '.join(ls)
    if alt_name in itol_keep_strains:
        print(line)
