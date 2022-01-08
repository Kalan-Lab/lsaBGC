import os
import sys
from collections import defaultdict

gcf_popgene_file = 'GCF_Ortholog_Group_Information.txt'
hg_matrix_file = '../../Processing/Process_Core_Genomes/Orthogroups.csv'
cluster_listings_dir = '../../Clustering/Cluster_i14j00r60/GCF_Listings/'

hg_samples = defaultdict(set)
samples = []
with open(hg_matrix_file) as ohmf:
    for i, line in enumerate(ohmf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: 
            samples = ls[1:]
        else:
            hg = ls[0]
            for j, v in enumerate(ls[1:]):
                if v.strip() != '':
                    s = samples[j]
                    hg_samples[hg].add(s)

all_samples = set(samples)
gcf_samples = defaultdict(set)
for f in os.listdir(cluster_listings_dir):
    gcf = f.split('.txt')[0]
    with open(cluster_listings_dir + f) as of:
        for line in of:
            line = line.strip()
            ls = line.split('\t')
            gcf_samples[gcf].add(ls[0])
            
with open(gcf_popgene_file) as opf:
    for i, line in enumerate(opf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gcf = ls[0]
        hg = ls[2]
        
        gcf_samps = gcf_samples[gcf]
        gcf_samps_comp = all_samples.difference(gcf_samps)

        hg_samps = hg_samples[hg]
        hg_samps_comp = all_samples.difference(hg_samps)
        if len(hg_samps.intersection(gcf_samps_comp)) == 0:
            print(gcf + '\t' + hg)
