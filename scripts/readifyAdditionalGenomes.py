#!/usr/bin/env python

### Program: readifyAdditionalGenomes.py
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



def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: readifyAdditionalGenomes.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Prepares additional genomes 
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-d', '--additional_genome_listing',
						help='Tab-delimited, two column file for samples with additional/draft\ngenomes (same format as for the "--genome_listing" argument). The genomes/BGCs of these\nsamples won\'t be used in ortholog-grouping of proteins and clustering of BGCs, but will simply have gene\ncalling run for them. This will enable more sensitive/expanded detection of GCF instances later\nusing lsaBGC-Expansion/AutoExpansion.\nCheck note above about available scripts to automatically create this.',
						required=True)
	parser.add_argument('-o', '--output_directory', help='Parent output/workspace directory.', required=True)
	parser.add_argument('-c', '--cores', type=int,
						help="Total number of cores/threads to use for running OrthoFinder2/prodigal.", required=False,
						default=1)
	args = parser.parse_args()
	return args

def readifyAdditionalGenomes():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	outdir = os.path.abspath(myargs.output_directory) + '/'
	additional_genome_listing_file = myargs.additional_genome_listing
	cores = myargs.cores

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)
	try:
		additional_genome_listing_file = os.path.abspath(additional_genome_listing_file)
		assert (os.path.isfile(additional_genome_listing_file))
	except:
		raise RuntimeError('Issue with reading genome listing file for samples with additional genomic assemblies.')

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)
	logObject.info("Saving parameters for future records.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [additional_genome_listing_file, outdir, cores]
	parameter_names = ["Additional Genome Listing File", "Output Directory", "Number of Cores"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Process Additional Genomes
	additional_sample_annotation_listing_file = outdir + 'Additional_Sample_Annotation_Files.txt'
	
	additional_sample_genomes, additional_format_prediction = util.parseSampleGenomes(additional_genome_listing_file, logObject)
	if additional_format_prediction == 'mixed':
		logObject.error(
			'Format of additional genomes provided is not consistently FASTA or Genbank, please check input.')
		raise RuntimeError(
			'Format of additional genomes provided is not consistently FASTA or Genbank, please check input.')

	additional_proteomes_directory = outdir + 'Predicted_Proteomes_Additional/'
	additional_genbanks_directory = outdir + 'Genomic_Genbanks_Additional/'
	util.setupReadyDirectory([additional_proteomes_directory, additional_genbanks_directory])
	if additional_format_prediction == 'fasta':
		additional_prodigal_outdir = outdir + 'Prodigal_Gene_Calling_Additional/'
		util.setupReadyDirectory([additional_prodigal_outdir])

		# Note, locus tags of length 4 are used within lsaBGC to mark samples with additional genomes where we ultimately
		# find them via lsaBGC-Expansion.
		util.processGenomes(additional_sample_genomes, additional_prodigal_outdir, additional_proteomes_directory,
							additional_genbanks_directory, logObject, cores=cores, locus_tag_length=4)
	else:
		# genomes are provided as Genbanks with CDS features
		gene_name_mapping_outdir = outdir + 'Mapping_of_New_Gene_Names_to_Original/'
		util.setupReadyDirectory([gene_name_mapping_outdir])
		util.processGenomesAsGenbanks(additional_sample_genomes, additional_proteomes_directory,
									  additional_genbanks_directory, gene_name_mapping_outdir, logObject,
									  cores=cores, locus_tag_length=4)

	additional_sample_annotation_listing_handle = open(additional_sample_annotation_listing_file, 'w')
	for f in os.listdir(additional_proteomes_directory):
		sample = f.split('.faa')[0]
		additional_sample_annotation_listing_handle.write(sample + '\t' + additional_genbanks_directory + sample + '.gbk' + '\t' + additional_proteomes_directory + f + '\n')
	additional_sample_annotation_listing_handle.close()

	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	readifyAdditionalGenomes()
