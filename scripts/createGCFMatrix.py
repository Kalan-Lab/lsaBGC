import os
import sys
from collections import defaultdict

gcf_dir = os.path.abspath(sys.argv[1]) + '/' 
sample_set = set([x.strip() for x in open(sys.argv[2]).readlines()])
gcf_det_file = sys.argv[3]

outf1 = 'GCF_to_Class.txt'
outh1 = open(outf1, 'w')

outf2 = 'GCF_Count.Matrix.txt'
outh2 = open(outf2, 'w')

hybrid_gcfs = set([])
gcf_annots = defaultdict(set)
with open(gcf_det_file) as ogdf:
    for i, line in enumerate(ogdf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        annots = set([])
        for a in ls[-1].split('; '):
            annots.add(a.split(':')[0])
        if 'NRPS' in annots and 'NRPS-like' in annots:
            annots.remove('NRPS-like')
        if len(annots) >= 2:
            gcf_annots[ls[0]] = 'hybrid'

            
            an_dict = {}
            tot_count = 0
            for an in ls[-1].split('; '):
                ans = [x.strip() for x in an.split(':')]
                an_dict[ans[0]] = float(ans[1])
                tot_count += float(ans[1])

            for an in an_dict:
                if an_dict[an]/float(tot_count) >= 0.90:
                    gcf_annots[ls[0]] = an 
        else:
            gcf_annots[ls[0]] = list(annots)[0]

outh2.write('X\t' + '\t'.join(sorted(sample_set)) + '\n')
for f in os.listdir(gcf_dir):
    g = f.split('.txt')[0]
    outh1.write(g + '\t ' + gcf_annots[g] + '\n')
    gcf_samples = set([])
    with open(gcf_dir + f) as of:
        for i, line in enumerate(of):
            if i == 0: continue
            line = line.strip()
            ls = line.split('\t')
            if ls[0] in sample_set:
                gcf_samples.add(ls[0]) 
    printlist = [g]
    for s in sorted(sample_set):
        if s in gcf_samples:
            printlist.append('1')
        else:
            printlist.append('0')
    outh2.write('\t'.join(printlist) + '\n')
outh1.close()
outh2.close()
