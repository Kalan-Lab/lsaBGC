import os
import sys

alg_dir = os.path.abspath(sys.argv[1]) + '/' # results directory of comparator BAM
res_dir = os.path.abspath(sys.argv[2]) + '/' # results directory for listing files

for g in os.listdir(alg_dir):
    outf = open(res_dir + g + '.txt', 'w')
    gcf_dir = alg_dir + g + '/'
    for f in os.listdir(gcf_dir):
        if not f.endswith('.sorted.bam'): continue
        s = f.split('.sorted.bam')[0]
        outf.write(s + '\t' + gcf_dir + f + '\n')
    outf.close()
