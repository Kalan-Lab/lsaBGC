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
from time import sleep
import argparse
from lsaBGC.classes.GCF import GCF
from lsaBGC.classes.Pan import Pan
from lsaBGC import util

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_TREESTRUCTURE = lsaBGC_main_directory + '/lsaBGC/Rscripts/createNJTreeAndDefineClades.R'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Automate.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Program to parallelize most of lsaBGC programs for each GCF 
	
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-o', '--output_directory', help="Parent output/workspace directory.", required=True)
	parser.add_argument('-l', '--initial_input_listing', type=str, help="Tab delimited text file. First column is the sample name and the second is the path to its assembly in FASTA format. Please remove troublesome characters in the sample name.", required=True)
	parser.add_argument('-g', '--gcf_listing_dir', help='Directory with GCF listing files.', required=True)
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder homolog group by sample matrix.", required=True)
	parser.add_argument('-e', '--expansion_input_listing', help="Path to tab delimited file listing: (1) sample name (2) path to Prokka Genbank and (3) path to Prokka predicted proteome. This file is produced by lsaBGC-Process.py.", required=False, default=None)
	parser.add_argument('-p', '--population_analysis',  action='store_true', help="Whether to construct species phylogeny and use it to determine populations.", default=False, required=False)
	parser.add_argument('-i', '--discovary_analysis_id', help="Identifier for novelty SNV mining analysis. Not providing this parameter will avoid running lsaBGC-DiscoVary step.", required=False, default=None)
	parser.add_argument('-n', '--discovary_input_listing', help="Path to tab delimited file listing: (1) sample name (2) path to forward reads and (3) path to reverse reads.", required=False, default=None)
	parser.add_argument('-c', '--cores', type=int, help="Total number of cores to use.", required=False, default=1)

	args = parser.parse_args()
	return args


def lsaBGC_Automate():
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
	initial_input_listing_file = os.path.abspath(myargs.initial_input_listing)

	try:
		assert(os.path.isdir(gcf_listing_dir) and os.path.isfile(original_orthofinder_matrix_file) and os.path.isfile(initial_input_listing_file))
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

	expansion_listing_file = myargs.expansion_input_listing
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
	parameter_values = [gcf_listing_dir, initial_input_listing_file, original_orthofinder_matrix_file, outdir,
						expansion_listing_file, population_analysis, discovary_analysis_id, discovary_input_listing,
						cores]
	parameter_names = ["GCF Listings Directory", "Listing File of Prokka Annotation Files for Initial Set of Samples",
					   "OrthoFinder Homolog Matrix", "Output Directory",
					   "Listing File of Prokka Annotation Files for Expansion/Additional Set of Samples",
					   "Delineate Populations and Perform Population Genetics Analytics", "DiscoVary Analysis ID",
					   "DiscoVary Sequencing Data Location Specification File", "Cores"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	"""
	if expansion_listing_file:
		try:
			assert (os.path.isfile(expansion_listing_file))
		except:
			raise RuntimeError('One or more of the input files provided, does not exist. Exiting now ...')

		# Create Pan object
		Pan_Object = Pan(expansion_listing_file, logObject=logObject)

		# Extract Genbank Sequences into FASTA
		logObject.info("Converting Genbanks from Expansion listing file into FASTA per sample.")
		gw_fasta_dir = outdir + 'Sample_Expansion_FASTAs/'
		if not os.path.isdir(gw_fasta_dir): os.system('mkdir %s' % gw_fasta_dir)
		gw_fasta_listing_file = outdir + 'Expansion_Genome_FASTA_Listings.txt'
		Pan_Object.convertGenbanksIntoFastas(gw_fasta_dir, gw_fasta_listing_file)
		logObject.info("Successfully performed conversions.")

		try:
			new_assembly_listing_file = outdir + assembly_listing_file.split('/')[-1]
			os.system('cp %s %s' % (assembly_listing_file, new_assembly_listing_file))
			assembly_listing_file = new_assembly_listing_file

			assembly_listing_handle = open(assembly_listing_file, 'a+')
			with open(gw_fasta_listing_file) as ogflf:
				for line in ogflf:
					assembly_listing_handle.write(line)
			assembly_listing_handle.close()
			logObject.info("Successfully updated genomes listing")
		except:
			raise RuntimeError('Issues with updating genomes listing')
	"""

	mash_nj_tree = None
	population_listing_file = None
	if population_analysis:
		# Run MASH Analysis Between Genomic Assemblies
		logObject.info("Running MASH Analysis Between Genomes.")
		gw_pairwise_differences = util.calculateMashPairwiseDifferences(assembly_listing_file, outdir, 'genome_wide',
																		10000, cores, logObject)
		logObject.info("Ran MASH Analysis Between Genomes.")

		mash_matrix_file = outdir + 'MASH_Distance_Matrix.txt'
		mash_matrix_handle = open(mash_matrix_file, 'w')
		mash_matrix_handle.write('Sample/Sample\t' + '\t'.join([s.replace(' ', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_") for s in sorted(gw_pairwise_differences)]) + '\n')

		for s1 in sorted(gw_pairwise_differences):
			printlist = [s1.replace(' ', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_")]
			for s2 in sorted(gw_pairwise_differences):
				printlist.append(str(gw_pairwise_differences[s1][s2]))
			mash_matrix_handle.write('\t'.join(printlist) + '\n')
		mash_matrix_handle.close()

		mash_nj_tree = outdir + 'MASH_NeighborJoining_Tree.nwk'
		population_listing_file = outdir + 'Populations_Defined.txt'
		populations_on_nj_tree_pdf = outdir + 'Populations_on_NJ_Tree.pdf'
		cmd = ['Rscript', RSCRIPT_FOR_TREESTRUCTURE, mash_matrix_file, mash_nj_tree, populations_on_nj_tree_pdf, population_listing_file]
		try:
			util.run_cmd(cmd, logObject)
		except Exception as e:
			logObject.error("Had issues with creating neighbor joining tree and defining populations using treestructure.")
			raise RuntimeError("Had issues with creating neighbor joining tree and defining populations using treestructure.")

	see_outdir = outdir + 'See/'
	pop_outdir = outdir + 'PopGene/'
	div_outdir = outdir + 'Divergence/'
	if not os.path.isdir(see_outdir): os.system('mkdir %s' % see_outdir)
	if not os.path.isdir(pop_outdir): os.system('mkdir %s' % pop_outdir)
	if not os.path.isdir(div_outdir): os.system('mkdir %s' % div_outdir)

	if expansion_listing_file:
		exp_outdir = outdir + 'Expansion/'
		if not os.path.isdir(exp_outdir): os.system('mkdir %s' % exp_outdir)

	if discovary_analysis_id and discovary_input_listing:
		dis_outdir = outdir + 'DiscoVary_' + '_'.join(discovary_analysis_id.split()) + '/'
		if not os.path.isdir(dis_outdir): os.system('mkdir %s' % dis_outdir)

	mash_nj_tree = "/home/salamzade/lsaBGC_Development/Micrococcus_luteus_Dataset/Processing_All/SpeciesTree_rooted.nwk"
	for g in os.listdir(gcf_listing_dir):
		gcf_id = g.split('.txt')[0]
		gcf_listing_file = gcf_listing_dir + g
		orthofinder_matrix_file = original_orthofinder_matrix_file
		logObject.info("Beginning analysis for GCF %s" % gcf_id)
		sys.stderr.write("Beginning analysis for GCF %s\n" % gcf_id)

		# 1. Run lsaBGC-Expansion.py
		if expansion_listing_file:
			gcf_exp_outdir = exp_outdir + gcf_id + '/'
			if True: #not os.path.isdir(gcf_exp_outdir):
				os.system('mkdir %s' % gcf_exp_outdir)
				cmd = ['lsaBGC-Expansion.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-l',
					   initial_input_listing_file, '-e', expansion_listing_file, '-o', gcf_exp_outdir, '-i', gcf_id,
					   '-c', str(cores)]
				try:
					util.run_cmd(cmd, logObject, stderr=sys.stderr, stdout=sys.stdout)
					assert (os.path.isfile(gcf_exp_outdir + 'Orthogroups.expanded.csv'))
					assert (os.path.isfile(gcf_exp_outdir + 'GCF_Expanded.txt'))
					orthofinder_matrix_file = gcf_exp_outdir + 'Orthogroups.expanded.csv'
					gcf_listing_file = gcf_exp_outdir + 'GCF_Expanded.txt'
				except Exception as e:
					logObject.warning("lsaBGC-Expansion.py was unsuccessful, skipping over GCF %s" % gcf_id)
					sys.stderr.write("lsaBGC-Expansion.py was unsuccessful, skipping over GCF %s\n" % gcf_id)
					continue
			else:
				try:
					assert(os.path.isfile(gcf_exp_outdir + 'Orthogroups.expanded.csv'))
					assert(os.path.isfile(gcf_exp_outdir + 'GCF_Expanded.txt'))
					orthofinder_matrix_file = gcf_exp_outdir + 'Orthogroups.expanded.csv'
					gcf_listing_file = gcf_exp_outdir + 'GCF_Expanded.txt'
					sys.stderr.write("lsaBGC-Expansion.py results directory already exists thus continuing using existing results.")
				except Exception as e:
					logObject.warning("lsaBGC-Expansion.py was unsuccessful, skipping over GCF %s" % gcf_id)
					sys.stderr.write("Warning: lsaBGC-Expansion.py was unsuccessful, skipping over GCF %s\n" % gcf_id)
					continue

		# 2. Run lsaBGC-See.py
		gcf_see_outdir = see_outdir + gcf_id + '/'
		if True: #not os.path.isdir(gcf_see_outdir):
			os.system('mkdir %s' % gcf_see_outdir)
			cmd = ['lsaBGC-See.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-o', gcf_see_outdir,
				   '-i', gcf_id, '-p', '-c', str(cores)]
			if mash_nj_tree:
				cmd += ['-s', mash_nj_tree]
			try:
				util.run_cmd(cmd, logObject)
			except Exception as e:
				logObject.warning("lsaBGC-See.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-See.py was unsuccessful for GCF %s\n" % gcf_id)

		"""
		# 3. Run lsaBGC-PopGene.py
		gcf_pop_outdir = pop_outdir + gcf_id + '/'
		if not os.path.isdir(gcf_pop_outdir):
			os.system('mkdir %s' % gcf_pop_outdir)
			cmd = ['lsaBGC-PopGene.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-o', gcf_pop_outdir,
				   '-i', gcf_id, '-c', str(cores)]
			if population_listing_file:
				cmd += ['-p', population_listing_file]
			try:
				util.run_cmd(cmd, logObject)
			except Exception as e:
				logObject.warning("lsaBGC-PopGene.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-PopGene.py was unsuccessful for GCF %s\n" % gcf_id)

		# 4. Run lsaBGC-Divergence.py
		gcf_div_outdir = div_outdir + gcf_id + '/'
		if not os.path.isdir(gcf_div_outdir):
			os.system('mkdir %s' % gcf_div_outdir)
			cmd = ['lsaBGC-Divergence.py', '-g', gcf_listing_file, '-a', assembly_listing_file, '-o', gcf_div_outdir,
				   '-i', gcf_id, '-c', str(cores)]
			if expansion_listing_file:
				cmd += ['-e', expansion_listing_file]
			try:
				util.run_cmd(cmd, logObject)
			except Exception as e:
				logObject.warning("lsaBGC-Divergence.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-Divergence.py was unsuccessful for GCF %s\n" % gcf_id)

		# 5. Run lsaBGC-DiscoVary.py
		if discovary_analysis_id and discovary_input_listing:
			gcf_dis_outdir = dis_outdir + gcf_id + '/'
			if not os.path.isdir(gcf_dis_outdir):
				os.system('mkdir %s' % gcf_dis_outdir)
				cmd = ['lsaBGC-DiscoVary.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-o',
					   gcf_dis_outdir, '-i', gcf_id, '-c', str(cores), '-p', discovary_input_listing, '-a',
					   gcf_pop_outdir + 'Codon_Alignments_Listings.txt']
				try:
					util.run_cmd(cmd, logObject)
				except Exception as e:
					logObject.warning("lsaBGC-DiscoVary.py was unsuccessful for GCF %s" % gcf_id)
					sys.stderr.write("Warning: lsaBGC-DiscoVary.py was unsuccessful for GCF %s\n" % gcf_id)
	"""
	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	lsaBGC_Automate()