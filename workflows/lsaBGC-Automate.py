#!/usr/bin/env python

### Program: lsaBGC-Automate.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2021, Kalan-Lab
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import sys
import argparse

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Automate.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Program to parallelize most of lsaBGC programs for each GCF 
	
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-o', '--output_directory', help="Parent output/workspace directory.", required=True)
	parser.add_argument('-g', '--gcf_listing_dir', help='Directory with GCF listing files.', required=True)
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder homolog group by sample matrix.", required=True)
	parser.add_argument('-s', '--species_phylogeny', help="OrthoFinder generated species phylogeny in newick format.")
	parser.add_argument('-a', '--assembly_listing', help="Tab delimited text file. First column is the sample name and the second is the path to its assembly in FASTA format. Please remove troublesome characters in the sample name.", required=True)
	parser.add_argument('-')
	parser.add_argument('-n', '--novelty_input', help="")

	parser.add_argument('-fo')
	parser.add_argument('-c', '--cores', type=int, help="Number of cores to use for MCL step.", required=False, default=1)

	args = parser.parse_args()
	return args

all_bgcs = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/Mluteus_Process_All_V4/All_AntiSMASH_BGCs.txt'
species_tree = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/Processing/SpeciesTree_rooted.txt'
og_orthogroups = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/Processing/Orthogroups.csv'
gcf_listing_dir = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/BGC_Clustering_I14_Jaccard//GCF_Listings/'
expansion_dir = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/GCF_Clustering_Expanded/'
see_dir = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/GCF_Visualization/'
popgene_dir = '/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/GCF_PopGen_Analyses/'

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