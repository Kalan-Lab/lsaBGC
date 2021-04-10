import os
import sys

og_orthogroups = '/home/salamzade/lsaBGC_Development/Streptomyces_Dataset/Processing/Orthogroups.csv'
gcf_listing_dir = '/home/salamzade/lsaBGC_Development/Streptomyces_Dataset/BGC_Clustering_Results_i28/GCF_Listings/'
popgene_dir = '/home/salamzade/lsaBGC_Development/Streptomyces_Dataset/GCF_PopGen_Analyses/'

lsabin_dir = "/home/salamzade/Git_Repos/lsaBGC/bin/"
popgene_prog = lsabin_dir + 'lsaBGC-PopGene.py'
see_prog = lsabin_dir + 'lsaBGC-See.py'
expan_prog = lsabin_dir + 'lsaBGC-HMMExpansion.py'

for i, g in enumerate(os.listdir(gcf_listing_dir)):
	if i == 0: continue
	gcf = g.split('.txt')[0]
	gcf_file = gcf_listing_dir + g

	# 3. run popgene
	pg_outdir = popgene_dir + gcf + '/'
	pg_expan_outdir = popgene_dir + gcf + '_Expanded/'
	os.system('python %s -g %s -m %s -o %s -c 30' % (popgene_prog, gcf_file, og_orthogroups, pg_outdir))