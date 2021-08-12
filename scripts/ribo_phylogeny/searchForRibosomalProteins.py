import os
import sys
from collections import defaultdict

annotation_dir = os.path.abspath(sys.argv[1]) + '/'
hmm_dir = os.path.abspath(sys.argv[2]) + '/'
result_dir = os.path.abspath(sys.argv[3]) + '/'
for f in os.listdir(annotation_dir):
    if not f.endswith('.faa'): continue
    s = f.split('.faa')[0]
    pep_file = annotation_dir + f 

    if not os.path.isfile(pep_file): continue
    for hmm in os.listdir(hmm_dir):
                if not hmm.endswith("bact.hmm") or hmm.startswith('.'): continue
                hmmdb = hmm_dir + hmm
                outf= result_dir + s + '_with_' + hmm.split('_bact.hmm')[0] + '.out'
                os.system('hmmsearch --tblout %s %s %s' % (outf, hmmdb, pep_file))
