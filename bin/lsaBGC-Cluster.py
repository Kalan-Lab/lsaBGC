#!/usr/bin/env python

### Program: lsaBGC-Cluster.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

import os
import sys
from time import sleep
import argparse
from lsaBGC import lsaBGC

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Cluster.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will cluster BGCs found by AntiSMASH using MCL based on similarity exhibited in ortholog group presence/
	absence data. Clustering uses MCL.""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-b', '--bgc_listings', help='BGC listing file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
	parser.add_argument('-o', '--output_directory', help="Output directory.", required=True)
	parser.add_argument('-c', '--cores', type=int, help="Number of cores to use for MCL step.", required=False, default=1)
	parser.add_argument('-i', '--mcl_inflation', type=float, help="Inflation parameter to be used for MCL.", required=False, default=1.4)
	parser.add_argument('-j', '--jaccard_cutoff', type=float, help="Cutoff for Jaccard similarity of homolog groups shared between two BGCs.", default=50.0)
	parser.add_argument('-r', '--run_parameter_tests', action='store_true', help="Run tests for selecting the best inflation parameter and jaccard for MCL analysis and exit.", default=False, required=False)
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

	cores = myargs.cores
	mcl_inflation = myargs.mcl_inflation
	jaccard_cutoff = myargs.jacard_cutoff
	run_parameter_tests = myargs.run_inflation_tests


	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = lsaBGC.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for future provedance.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [bgc_listings_file, orthofinder_matrix_file, outdir, cores, mcl_inflation, jaccard_cutoff, run_parameter_tests]
	parameter_names = ["BGC Listing File", "OrthoFinder Orthogroups.csv File", "Output Directory", "Cores",
					   "MCL Inflation Parameter", "Jaccard Similarity Cutoff", "Run Inflation Parameter Tests?"]
	lsaBGC.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Parse BGCs from Listing File
	logObject.info("Starting to process BGC Genbanks from listing file.")
	bgc_sample, bgc_product, bgc_genes, all_genes = lsaBGC.readInBGCGenbanksComprehensive(bgc_listings_file, logObject)
	logObject.info("Successfully parsed BGC Genbanks.")

	# Step 2: Parse OrthoFinder Homolog vs Sample Matrix
	logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
	gene_to_cog, cog_genes, prop_multi_copy = lsaBGC.parseOrthoFinderMatrix(orthofinder_matrix_file, all_genes, calc_prop_multicopy=True)
	logObject.info("Successfully parsed homolog matrix.")

	# Step 3: Calculate overlap in homolog profiles between pairs of BGCs and prepare for MCL
	logObject.info('Calculating overlap in single copy ortholog Groups between BGC GBKs.')
	mcl_outdir = outdir + 'MCL_tmp_files/'
	if not run_parameter_tests: mcl_outdir = outdir
	elif not os.path.isdir(mcl_outdir): os.system('mkdir %s' % mcl_outdir)
	bgc_cogs, pairwise_relations, pair_relations_txt_file = lsaBGC.calculateBGCPairwiseRelations(bgc_genes, gene_to_cog, prop_multi_copy, mcl_outdir, logObject)
	logObject.info("Successfully calculated pairwise distances between BGCs based on homolog profiles.")

	# Step 4: Run MCL clustering, iterating through multiple inflation parameters if necessary.
	logObject.info('Starting to run MCL for finding Gene Cluster Families (GCFs)!')
	# Create and write to file which will detail the GCFs found from MCL clustering
	stats_file = outdir + 'GCF_details.txt'
	sf_handle = open(stats_file, 'w')
	if run_parameter_tests:
		sf_handle.write('\t'.join(['inflation parameter', 'GCF id', 'number of BGCs',
								   'samples with multiple BGCs in GCF', 'size of the SCC', 'mean number of OGs',
								   'stdev for number of OGs', 'min difference', 'max difference',
								   'annotations']) + '\n')
	else:
		sf_handle.write('\t'.join(['GCF id', 'number of BGCs', 'samples with multiple BGCs in GCF',
								   'size of the SCC', 'mean number of OGs', 'stdev for number of OGs',
								   'min difference', 'max difference', 'annotations']) + '\n')

	mcl_inflation_params = [mcl_inflation]
	jaccard_cutoff_params = [jaccard_cutoff]
	if run_parameter_tests:
		mcl_inflation_params = [0.8, 1.4, 2, 2.5, 3, 3.5, 4, 5]
		jaccard_cutoff_params = [0, 20, 30, 50, 75, 90]
	for mip in mcl_inflation_params:
		for jcp in jaccard_cutoff_params:
			lsaBGC.runMCLAndReportGCFs(mip, jcp, mcl_outdir, sf_handle, pairwise_relations, pair_relations_txt_file, bgc_cogs, bgc_product, bgc_sample, run_parameter_tests, cores, logObject)
	sf_handle.close()

	logObject.info("Successfully ran MCL clustering analysis to determine GCFs!")

	# Close logging object and exit
	lsaBGC.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	lsaBGC_Cluster()