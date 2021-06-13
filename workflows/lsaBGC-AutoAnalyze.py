#!/usr/bin/env python

### Program: lsaBGC-AutoAnalyze.py
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
from lsaBGC.classes.GCF import GCF
from lsaBGC.classes.Pan import Pan
from lsaBGC import util

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_NJTREECONSTRUCTION = lsaBGC_main_directory + '/lsaBGC/Rscripts/createNJTree.R'
RSCRIPT_FOR_DEFINECLADES = lsaBGC_main_directory + '/lsaBGC/Rscripts/defineClades.R'

annotated_gcfs = set(['GCF_6', 'GCF_8', 'GCF_9', 'GCF_10', 'GCF_11', 'GCF_12'])

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Automate.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Program to parallelize most of lsaBGC programs for each GCF 
	
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-o', '--output_directory', help="Parent output/workspace directory.", required=True)
	parser.add_argument('-l', '--input_listing', help="Path to tab delimited file listing: (1) sample name (2) path to Prokka Genbank and (3) path to Prokka predicted proteome. This file is produced by lsaBGC-Process.py.", required=False, default=None)
	parser.add_argument('-g', '--gcf_listing_dir', help='Directory with GCF listing files.', required=True)
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder homolog group by sample matrix.", required=True)
	parser.add_argument('-p', '--population_analysis',  action='store_true', help="Whether to construct species phylogeny and use it to determine populations.", default=False, required=False)
	parser.add_argument('-s', '--lineage_phylogeny', help="Path to species phylogeny. If not provided a MASH based neighborjoining tree will be constructed and used.", default=None, required=False)
	parser.add_argument('-i', '--discovary_analysis_id', help="Identifier for novelty SNV mining analysis. Not providing this parameter will avoid running lsaBGC-DiscoVary step.", required=False, default=None)
	parser.add_argument('-n', '--discovary_input_listing', help="Path to tab delimited file listing: (1) sample name (2) path to forward reads and (3) path to reverse reads.", required=False, default=None)
	parser.add_argument('-c', '--cores', type=int, help="Total number of cores to use.", required=False, default=1)

	args = parser.parse_args()
	return args

def writeToOpenHandle(gcf_results_file, combined_results_handle, include_header):
	with open(gcf_results_file) as ogrf:
		for i, line in enumerate(ogrf):
			if i == 0 and include_header:
				combined_results_handle.write(line)
			elif i != 0:
				combined_results_handle.write(line)

def lsaBGC_AutoAnalyze():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	outdir = os.path.abspath(myargs.output_directory) + '/'
	gcf_listing_dir = os.path.abspath(myargs.gcf_listing_dir) + '/'
	original_orthofinder_matrix_file = os.path.abspath(myargs.orthofinder_matrix)
	input_listing_file = os.path.abspath(myargs.input_listing)

	try:
		assert(os.path.isdir(gcf_listing_dir) and os.path.isfile(original_orthofinder_matrix_file) and os.path.isfile(input_listing_file))
	except:
		raise RuntimeError('Input directory with GCF listings does not exist or the OrthoFinder does not exist. Exiting now ...')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Continuing in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	"""
	PARSE OPTIONAL INPUTS
	"""

	lineage_phylogeny_file = myargs.lineage_phylogeny
	population_analysis = myargs.population_analysis
	discovary_analysis_id = myargs.discovary_analysis_id
	discovary_input_listing = myargs.discovary_input_listing
	cores = myargs.cores

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for easier determination of results' provenance in the future.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [gcf_listing_dir, input_listing_file, original_orthofinder_matrix_file, outdir,
						lineage_phylogeny_file, population_analysis, discovary_analysis_id, discovary_input_listing,
						cores]
	parameter_names = ["GCF Listings Directory", "Listing File of Prokka Annotation Files for Initial Set of Samples",
					   "OrthoFinder Homolog Matrix", "Output Directory", "Phylogeny File in Newick Format",
					   "Delineate Populations and Perform Population Genetics Analytics", "DiscoVary Analysis ID",
					   "DiscoVary Sequencing Data Location Specification File", "Cores"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	Pan_Object = Pan(input_listing_file, logObject=logObject)
	logObject.info("Converting Genbanks from Expansion listing file into FASTA per sample.")
	gw_fasta_dir = outdir + 'Sample_Assemblies_in_FASTA/'
	if not os.path.isdir(gw_fasta_dir): os.system('mkdir %s' % gw_fasta_dir)
	gw_fasta_listing_file = outdir + 'Genome_FASTA_Listings.txt'
	Pan_Object.convertGenbanksIntoFastas(gw_fasta_dir, gw_fasta_listing_file)
	logObject.info("Successfully performed conversions.")

	logObject.info("Running MASH Analysis Between Genomes.")
	gw_pairwise_differences = util.calculateMashPairwiseDifferences(gw_fasta_listing_file, outdir, 'genome_wide',
																		10000, cores, logObject)
	logObject.info("Ran MASH Analysis Between Genomes.")
	mash_matrix_file = outdir + 'MASH_Distance_Matrix.txt'
	mash_matrix_handle = open(mash_matrix_file, 'w')
	mash_matrix_handle.write('Sample/Sample\t' + '\t'.join([s for s in sorted(gw_pairwise_differences)]) + '\n')

	for s1 in sorted(gw_pairwise_differences):
		printlist = [s1]
		for s2 in sorted(gw_pairwise_differences):
			printlist.append(str(gw_pairwise_differences[s1][s2]))
		mash_matrix_handle.write('\t'.join(printlist) + '\n')
	mash_matrix_handle.close()

	if not lineage_phylogeny_file:
		# Run MASH Analysis Between Genomic Assemblies
		logObject.info("Running MASH Analysis Between Genomes.")
		mash_nj_tree = outdir + 'MASH_NeighborJoining_Tree.nwk'

		# create neighbor-joining tree
		cmd = ['Rscript', RSCRIPT_FOR_NJTREECONSTRUCTION, mash_matrix_file, mash_nj_tree]
		try:
			util.run_cmd(cmd, logObject)
			assert(util.is_newick(mash_nj_tree))
			lineage_phylogeny_file = mash_nj_tree
		except Exception as e:
			logObject.error("Had issues with creating neighbor joining tree and defining populations using treestructure.")
			raise RuntimeError("Had issues with creating neighbor joining tree and defining populations using treestructure.")

	population_listing_file = None
	if population_analysis:
		population_listing_file = outdir + 'Populations_Defined.txt'
		populations_on_nj_tree_pdf = outdir + 'Populations_on_NJ_Tree.pdf'
		# create neighbor-joining tree
		cmd = ['Rscript', RSCRIPT_FOR_DEFINECLADES, lineage_phylogeny_file, mash_matrix_file, population_listing_file, populations_on_nj_tree_pdf]
		try:
			util.run_cmd(cmd, logObject)
		except Exception as e:
			logObject.error(
				"Had issues with creating neighbor joining tree and defining populations using treestructure.")
			raise RuntimeError(
				"Had issues with creating neighbor joining tree and defining populations using treestructure.")

	see_outdir = outdir + 'See/'
	pop_outdir = outdir + 'PopGene/'
	div_outdir = outdir + 'Divergence/'
	if not os.path.isdir(see_outdir): os.system('mkdir %s' % see_outdir)
	if not os.path.isdir(pop_outdir): os.system('mkdir %s' % pop_outdir)
	if not os.path.isdir(div_outdir): os.system('mkdir %s' % div_outdir)

	if discovary_analysis_id and discovary_input_listing:
		dis_outdir = outdir + 'DiscoVary_' + '_'.join(discovary_analysis_id.split()) + '/'
		if not os.path.isdir(dis_outdir): os.system('mkdir %s' % dis_outdir)

	for g in os.listdir(gcf_listing_dir):
		gcf_id = g.split('.txt')[0]
		gcf_listing_file = gcf_listing_dir + g
		orthofinder_matrix_file = original_orthofinder_matrix_file
		logObject.info("Beginning analysis for GCF %s" % gcf_id)
		sys.stderr.write("Beginning analysis for GCF %s\n" % gcf_id)

		# 1. Run lsaBGC-See.py
		gcf_see_outdir = see_outdir + gcf_id + '/'
		if not os.path.isdir(gcf_see_outdir):
			os.system('mkdir %s' % gcf_see_outdir)
			cmd = ['lsaBGC-See.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-o', gcf_see_outdir,
				   '-i', gcf_id, '-s', lineage_phylogeny_file, '-p', '-c', str(cores)]
			try:
				util.run_cmd(cmd, logObject)
			except Exception as e:
				logObject.warning("lsaBGC-See.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-See.py was unsuccessful for GCF %s\n" % gcf_id)

		# 2. Run lsaBGC-PopGene.py
		gcf_pop_outdir = pop_outdir + gcf_id + '/'
		if not os.path.isdir(gcf_pop_outdir):
			os.system('mkdir %s' % gcf_pop_outdir)
			cmd = ['lsaBGC-PopGene.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-o', gcf_pop_outdir,
				   '-i', gcf_id, '-c', str(cores)]
			if population_listing_file:
				cmd += ['-p', population_listing_file]
			try:
				util.run_cmd(cmd, logObject, stderr=sys.stderr, stdout=sys.stdout)
			except Exception as e:
				logObject.warning("lsaBGC-PopGene.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-PopGene.py was unsuccessful for GCF %s\n" % gcf_id)

		# 3. Run lsaBGC-Divergence.py
		gcf_div_outdir = div_outdir + gcf_id + '/'
		if not os.path.isdir(gcf_div_outdir):
			os.system('mkdir %s' % gcf_div_outdir)
			cmd = ['lsaBGC-Divergence.py', '-g', gcf_listing_file, '-l', input_listing_file, '-o', gcf_div_outdir,
				   '-i', gcf_id, '-c', str(cores)]
			try:
				util.run_cmd(cmd, logObject)
			except Exception as e:
				logObject.warning("lsaBGC-Divergence.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-Divergence.py was unsuccessful for GCF %s\n" % gcf_id)

		# 4. Run lsaBGC-DiscoVary.py
		if discovary_analysis_id and discovary_input_listing:
			gcf_dis_outdir = dis_outdir + gcf_id + '/'
			if gcf_id in annotated_gcfs: #not os.path.isdir(gcf_dis_outdir):
				os.system('mkdir %s' % gcf_dis_outdir)
				cmd = ['lsaBGC-DiscoVary.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-o',
					   gcf_dis_outdir, '-i', gcf_id, '-c', str(cores), '-p', discovary_input_listing, '-a',
					   gcf_pop_outdir + 'Codon_Alignments_Listings.txt']
				try:
					util.run_cmd(cmd, logObject, stderr=sys.stderr, stdout=sys.stdout)
				except Exception as e:
					logObject.warning("lsaBGC-DiscoVary.py was unsuccessful for GCF %s" % gcf_id)
					sys.stderr.write("Warning: lsaBGC-DiscoVary.py was unsuccessful for GCF %s\n" % gcf_id)

	combined_orthoresults_refined_handle = open(outdir + 'GCF_Ortholog_Group_Information.txt', 'w')
	combined_orthoresults_pop_refined_handle = open(outdir + 'Population_GCF_Ortholog_Group_Information.txt', 'w')
	combined_orthoresults_unrefined_handle = open(outdir + 'GCF_Ortholog_Group_Information_Mad_Refined.txt', 'w')
	combined_orthoresults_pop_unrefined_handle = open(outdir + 'Population_GCF_Ortholog_Group_Information_Mad_Refined.txt', 'w')
	combined_divergence_results_handle = open(outdir + 'GCF_Divergences.txt', 'w')
	for i, g in enumerate(os.listdir(gcf_listing_dir)):
		gcf_id = g.split('.txt')[0]
		gcf_pop_outdir = pop_outdir + gcf_id + '/'
		gcf_div_outdir = div_outdir + gcf_id + '/'
		gcf_pop_results_refined = gcf_pop_outdir + 'Ortholog_Group_Information_MAD_Refined.txt'
		gcf_pop_results_unrefined = gcf_pop_outdir + 'Ortholog_Group_Information.txt'

		jref = 0
		junref = 0
		for f in os.listdir(gcf_pop_outdir):
			include_header = False
			popgene_result_file = gcf_pop_outdir + f
			if 'Ortholog_Group_Information' in f and 'Pop-' in f and not 'Ortholog_Group_Information_MAD_Refined' in f:
				if jref == 0 and i == 0: include_header = True
				writeToOpenHandle(popgene_result_file, combined_orthoresults_pop_refined_handle, include_header)
				jref += 1
			if 'Ortholog_Group_Information_MAD_Refined' in f and 'Pop-' in f:
				if junref == 0 and i == 0: include_header = True
				writeToOpenHandle(popgene_result_file, combined_orthoresults_pop_unrefined_handle, include_header)
				junref += 1

		gcf_div_results = gcf_div_outdir + 'Relative_Divergence_Report.txt'

		include_header = False
		if i == 0: include_header = True
		writeToOpenHandle(gcf_pop_results_refined, combined_orthoresults_refined_handle, include_header)
		writeToOpenHandle(gcf_pop_results_unrefined, combined_orthoresults_unrefined_handle, include_header)
		writeToOpenHandle(gcf_div_results, combined_divergence_results_handle, include_header)

	combined_orthoresults_refined_handle.close()
	combined_orthoresults_unrefined_handle.close()
	combined_divergence_results_handle.close()

	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	lsaBGC_AutoAnalyze()
