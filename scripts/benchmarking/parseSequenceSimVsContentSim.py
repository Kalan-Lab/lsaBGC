import os
import sys
from Bio import SeqIO
from collections import defaultdict
import statistics
from scipy import stats

popgene_dir = os.path.abspath(sys.argv[1]) + '/'
bgsee_dir = os.path.abspath(sys.argv[2]) + '/'

print('\t'.join(['GCF', 'Sample1', 'Sample2', 'Avg_Matching_Identity', 'Shared_Homologs_Proportion', 'Syntenic_Similarity']) )
for gcf in os.listdir(popgene_dir):
    itol_file = bgsee_dir + gcf + '/BGCs_Visualization.iTol.txt'
    if not os.path.isfile(itol_file): continue



    cod_align_dir = popgene_dir + gcf + '/Codon_Alignments/' 
    sample_cogs = defaultdict(set)
    sample_cog_seqs = defaultdict(dict)
    for f in os.listdir(cod_align_dir):
        cog = f.split('.msa')[0]
        with open(cod_align_dir + f) as of:
            for rec in SeqIO.parse(of, 'fasta'):
                sample = rec.id.split('|')[0]
                sample_cogs[sample].add(cog)
                sample_cog_seqs[sample][cog] = str(rec.seq)

    bad_samples = set([])
    cog_order = defaultdict(lambda: defaultdict(list))
    with open(itol_file) as oif:
        for line in oif:
            if line.startswith("Micrococcus"):
                line = line.strip()
                ls = line.split('\t')
                samp = ls[0]
                if not samp in sample_cogs.keys():
                    if samp[:-2] in sample_cogs.keys():
                        bad_samples.add(samp[:-2])
                    elif samp[:-3] in sample_cogs.keys():
                        bad_samples.add(samp[:-3])
                for i, x in enumerate(ls[2:]):
                    cog_order[ls[0]][x.split('|')[-1]].append(i)

    valid = set(['A', 'C', 'G', 'T'])
    for i, s1 in enumerate(set(sample_cogs.keys()).difference(bad_samples)):
        s1cogs = sample_cogs[s1]
        for j, s2 in enumerate(set(sample_cogs.keys()).difference(bad_samples)):
            if i < j: continue
            s2cogs = sample_cogs[s2]
            shared = s1cogs.intersection(s2cogs)
            shared_prop = float(len(shared))/min(len(s1cogs), len(s2cogs))
            match_ids = []
            for cog in shared:
                s1seq = sample_cog_seqs[s1][cog].upper()
                s2seq = sample_cog_seqs[s2][cog].upper()
                nongap = 0
                match = 0
                for pos, b1 in enumerate(s1seq):
                    b2 = s2seq[pos]
                    if b1 in valid or b2 in valid: nongap += 1
                    if b1 in valid and b2 in valid and b1 == b2:
                        match += 1
                match_id = float(match)/nongap
                match_ids.append(match_id)
            if len(shared) >= 3:
                avg_match_id = statistics.mean(match_ids)

                cog1_order_lists = []
                cog2_order_lists = []
                for cog in shared:
                    for ind1 in cog_order[s1][cog]:
                        for ind2 in cog_order[s2][cog]:
                            cog1_order_lists.append(ind1)
                            cog2_order_lists.append(ind2)

                r, p = stats.pearsonr(cog1_order_lists, cog2_order_lists)
                abs_r = abs(r)
                print('\t'.join([str(x) for x in [gcf, s1, s2, avg_match_id, shared_prop, abs_r]]))
