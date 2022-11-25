#!/usr/bin/env python

### Program: lsaBGC-Cluster.py
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
from time import sleep
import argparse
from lsaBGC import util
from lsaBGC.classes.Pan import Pan

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Cluster.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will cluster BGC Genbanks using MCL based on similarity exhibited in ortholog group presence/
	absence data. Clustering uses MCL.""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-b', '--bgc_listings', help='BGC listing file. Tab delimited 2-column file: (1) sample name (2) path to predicted BGC in Genbank format.', required=True)
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
	parser.add_argument('-o', '--output_directory', help="Output directory.", required=True)
	parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is antiSMASH.', default='antiSMASH', required=False)
	parser.add_argument('-c', '--cpus', type=int, help="Number of cpus to use for MCL step.", required=False, default=1)
	parser.add_argument('-i', '--mcl_inflation', type=float, help="Inflation parameter to be used for MCL.", required=False, default=1.4)
	parser.add_argument('-j', '--jaccard_cutoff', type=float, help="Cutoff for Jaccard similarity of homolog groups shared between two BGCs.", required=False, default=50.0)
	parser.add_argument('-r', '--syntenic_correlation_cutoff', type=float, help="Minimum absolute correlation coefficient between two BGCs.", required=False, default=0.0)
	parser.add_argument('-s', '--split_by_annotation', action='store_true', help="Partition BGCs into groups based on annotation first.", required=False, default=False)
	parser.add_argument('-t', '--run_parameter_tests', action='store_true', help="Run tests for selecting the best inflation parameter and jaccard for MCL analysis and exit.", required=False, default=False)
	args = parser.parse_args()
	return args

def lsaBGC_Cluster():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	bgc_listings_file = os.path.abspath(myargs.bgc_listings)
	orthofinder_matrix_file = os.path.abspath(myargs.orthofinder_matrix)
	outdir = os.path.abspath(myargs.output_directory) + '/'

	### vet input files quickly
	try:
		assert(os.path.isfile(orthofinder_matrix_file))
		assert(os.path.isfile(bgc_listings_file))
	except:
		raise RuntimeError('One or more of the input files provided, does not exist. Exiting now ...')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	"""
	PARSE OPTIONAL INPUTS
	"""

	bgc_prediction_software = myargs.bgc_prediction_software.upper()
	cpus = myargs.cpus
	mcl_inflation = myargs.mcl_inflation
	jaccard_cutoff = myargs.jaccard_cutoff
	syntenic_correlation_cutoff = myargs.syntenic_correlation_cutoff
	run_parameter_tests = myargs.run_parameter_tests
	split_by_annotation = myargs.split_by_annotation

	try:
		assert (bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
	except:
		raise RuntimeError('BGC prediction software option is not a valid option.')

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)
	version_string = util.parseVersionFromSetupPy()
	logObject.info('Running lsaBGC version %s' % version_string)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for future records.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [bgc_listings_file, orthofinder_matrix_file, outdir, bgc_prediction_software, cpus, mcl_inflation,
						jaccard_cutoff, syntenic_correlation_cutoff, run_parameter_tests, split_by_annotation]
	parameter_names = ["BGC Listing File", "OrthoFinder Orthogroups.csv File", "Output Directory",
					   "BGC Prediction Method", "cpus", "MCL Inflation Parameter", "Jaccard Similarity Cutoff",
					   "Syntenic Correlation Coefficient Cutoff", "Run Inflation Parameter Tests?",
					   "Split BGCs into Annotation Categories First Prior to Clustering?"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Create Pan object
	Pan_Object = Pan(bgc_listings_file, logObject=logObject)

	# Step 1: Parse BGCs from Listing File
	logObject.info("Starting to process BGC Genbanks from listing file.")
	Pan_Object.readInBGCGenbanks(comprehensive_parsing=False, prediction_method=bgc_prediction_software)
	logObject.info("Successfully parsed BGC Genbanks.")

	# Step 2: Parse OrthoFinder Homolog vs Sample Matrix
	logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
	gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy = util.parseOrthoFinderMatrix(orthofinder_matrix_file, Pan_Object.pan_genes)
	Pan_Object.inputHomologyInformation(gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy)
	logObject.info("Successfully parsed homolog matrix.")

	# Step 3: Calculate overlap in homolog profiles between pairs of BGCs and prepare for MCL
	logObject.info('Calculating overlap in single copy ortholog Groups between BGC GBKs.')
	mcl_outdir = outdir + 'MCL_intermediate_files/'
	if not run_parameter_tests: mcl_outdir = outdir
	elif not os.path.isdir(mcl_outdir): os.system('mkdir %s' % mcl_outdir)
	Pan_Object.calculateBGCPairwiseRelations(mcl_outdir, split_by_annotation=split_by_annotation)
	logObject.info("Successfully calculated pairwise distances between BGCs based on homolog profiles.")

	# Step 4: Run MCL clustering, iterating through multiple inflation parameters if necessary.
	logObject.info('Starting to run MCL for finding Gene Cluster Families (GCFs)!')
	# Create and write to file which will detail the GCFs found from MCL clustering
	Pan_Object.openStatsFile(outdir, run_parameter_tests=run_parameter_tests)
	mcl_inflation_params = [mcl_inflation]
	jaccard_cutoff_params = [jaccard_cutoff]
	if run_parameter_tests:
		mcl_inflation_params = [0.8, 1.4, 2, 2.5, 3, 3.5, 4, 5]
		jaccard_cutoff_params = [0, 20, 30, 50, 75, 90]
	for mip in mcl_inflation_params:
		for jcp in jaccard_cutoff_params:
			Pan_Object.runMCLAndReportGCFs(mip, jcp, syntenic_correlation_cutoff, mcl_outdir, run_parameter_tests=run_parameter_tests, cpus=cpus)

	if run_parameter_tests:
		Pan_Object.plotResultsFromUsingDifferentParameters(outdir)

	logObject.info("Successfully ran MCL clustering analysis to determine GCFs!")

	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	lsaBGC_Cluster()