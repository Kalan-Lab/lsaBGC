#!/usr/bin/env python

### Program: lsaBGC-Ready.py
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
from Bio import SeqIO
from collections import defaultdict
from time import sleep
from lsaBGC import util
import subprocess
import traceback
import multiprocessing
import math
from ete3 import Tree

os.environ['OMP_NUM_THREADS'] = '4'

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'


def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Easy.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Downloads genomes for a taxa from ncbi using Kai Blin's ncbi-genome-download, performs BGC predictions 
	using GECCO (very-light weight dependency), and then runs the full lsaBGC suite to generate a "quick" (not 
	guaranteed, depends on the number of cores you have) and "easy" (I think its pretty easy - but obviously there 
	are considerations here and population genetics/evolutionary statistics will need to be interpretted with caution
	as there could be alot of replicate genomes for the taxa!!!).
	
	Completed genomes for the taxa will be downloaded (assuming > 2 exist) after which additional instances of 
	GCFs (gene-cluster families - groups of homologous BGCs) will be downloaded 
	
	A SECOND WARNING IN ALL CAPS THAT EVOLUTIONARY STATS SHOULD BE VIEWED WITH CAUTION IN THE ABSENCE OF DEREPLICATION!
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-t', '--taxa_name', help='Name of the taxa of interest.', required=True)
	parser.add_argument('-o', '--output_directory', help='Parent output/workspace directory.', required=True)
	parser.add_argument('-c', '--cores', type=int,
						help="Total number of cores/threads to use for running OrthoFinder2/prodigal.", required=False,
						default=1)
	args = parser.parse_args()
	return args

def lsaBGC_Easy():
# Step 0: Uncompress test_case.tar.gz and cd into it.
rm -rf test_case/
tar -zxvf test_case.tar.gz
cd test_case/

# Step 1: create genome listing inputs for lsaBGC-Ready.py
listAllGenomesInDirectory.py -i Primary_Genomes/ > Primary_Genomes_Listing.txt
listAllGenomesInDirectory.py -i Additional_Genomes > Additional_Genomes_Listing.txt

# Step 2: create BGC prediction Genbank listing input for lsaBGC-Ready.py (using AntiSMASH)
listAllBGCGenbanksInDirectory.py -i Primary_Genome_AntiSMASH_Results/ -p antiSMASH \
	-f  > Primary_Genome_BGC_Genbanks_Listing.txt

# Step 3: run lsaBGC-Ready.py - with clustering of primary genome BGCs, expansion to
#         additional genomes, and phylogeny construction set to automatically run.
FILE=../db/database_location_paths.txt
if [ -f "$FILE" ]; then
    lsaBGC-Ready.py -i Primary_Genomes_Listing.txt -d Additional_Genomes_Listing.txt \
   	-l Primary_Genome_BGC_Genbanks_Listing.txt -p antiSMASH -m BGC_Only \
	  -c 40 -t -a -lc -le -o lsaBGC_Ready_Results/
else
    lsaBGC-Ready.py -i Primary_Genomes_Listing.txt -d Additional_Genomes_Listing.txt \
	  -l Primary_Genome_BGC_Genbanks_Listing.txt -p antiSMASH -m BGC_Only \
	  -c 40 -t -lc -le -o lsaBGC_Ready_Results/
fi

# Step 4: run lsaBGC-AutoAnalyze.py - automatically run analytical programs for
#         visualization, evolutionary stats computation and single row per
#         homolog group/gene tables. Can also perform metagenomic analysis
#         if requested.
lsaBGC-AutoAnalyze.py -i lsaBGC_Ready_Results/Final_Results/Expanded_Sample_Annotation_Files.txt \
	-g lsaBGC_Ready_Results/Final_Results/Expanded_GCF_Listings/ \
	-m lsaBGC_Ready_Results/Final_Results/Expanded_Orthogroups.tsv \
	-s lsaBGC_Ready_Results/Final_Results/GToTree_output.tre \
	-w lsaBGC_Ready_Results/Final_Results/GToTree_Expected_Similarities.txt \
	-k lsaBGC_Ready_Results/Final_Results/Samples_in_GToTree_Tree.txt \
	-u Genome_to_Species_Mapping.txt -c 40 -o lsaBGC_AutoAnalyze_Results/ \
