import os
import sys

all_bgcs = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/Mluteus_Process_All_V4/All_AntiSMASH_BGCs.txt'
species_tree = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/Processing/SpeciesTree_rooted.txt'
og_orthogroups = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/Processing/Orthogroups.csv'
gcf_listing_dir = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/GCF_Clustering_Results/GCF_Listings/'
expansion_dir = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/GCF_Clustering_Expanded/'
see_dir = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/Visualization_Results/'
popgene_dir = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/GCF_PopGen_Analyses/'

lsabin_dir = "/home/salamzade/Git_Repos/lsaBGC/bin/"
popgene_prog = lsabin_dir + 'lsaBGC-PopGene.py'
see_prog = lsabin_dir + 'lsaBGC-See.py'
expan_prog = lsabin_dir + 'lsaBGC-HMMExpansion.py'

for i, g in enumerate(os.listdir(gcf_listing_dir)):
	if i == 0: continue
	gcf = g.split('.txt')[0]
	gcf_file = gcf_listing_dir + g

	# 1. run hmm expansion
	hmmexpan_outdir = expansion_dir + gcf + '/'
	os.system('python %s -g %s -m %s -a %s -o %s -c 20' % (expan_prog, gcf_file, og_orthogroups, all_bgcs, hmmexpan_outdir))
	hmmexpan_orthogroups = hmmexpan_outdir + 'Orthogroups.expanded.csv'
	hmmexpan_gcf_file = hmmexpan_outdir + 'GCF_Expanded.txt'

	# 2. run see
	see_outdir = see_dir + gcf + '/'
	see_expan_outdir = see_dir + gcf + '_Expanded/'
	os.system('python %s -g %s -m %s -o %s -l %s -s %s -c 20 -p' % (see_prog, gcf_file, og_orthogroups, see_outdir, gcf, species_tree))
	os.system('python %s -g %s -m %s -o %s -l %s -c 20 -p' % (see_prog, hmmexpan_gcf_file, hmmexpan_orthogroups, see_expan_outdir, gcf + '_Expanded'))

	# 3. run popgene
	pg_outdir = popgene_dir + gcf + '/'
	pg_expan_outdir = popgene_dir + gcf + '_Expanded/'
	os.system('python %s -g %s -m %s -o %s -c 20' % (popgene_prog, gcf_file, og_orthogroups, pg_outdir))
	os.system('python %s -g %s -m %s -o %s -c 20' % (popgene_prog, hmmexpan_gcf_file, hmmexpan_orthogroups, pg_expan_outdir))