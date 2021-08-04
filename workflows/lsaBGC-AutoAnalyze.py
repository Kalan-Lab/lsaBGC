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
#	 list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#	 this list of conditions and the following disclaimer in the documentation
#	 and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#	 contributors may be used to endorse or promote products derived from
#	 this software without specific prior written permission.
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
from operator import itemgetter
import argparse
from collections import defaultdict
from ete3 import Tree
from lsaBGC.classes.Pan import Pan
from lsaBGC import util

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_NJTREECONSTRUCTION = lsaBGC_main_directory + '/lsaBGC/Rscripts/createNJTree.R'
RSCRIPT_FOR_DEFINECLADES_FROM_PHYLO = lsaBGC_main_directory + '/lsaBGC/Rscripts/defineCladesFromPhylo.R'
RSCRIPT_FOR_BIGPICTUREHEATMAP = lsaBGC_main_directory + '/lsaBGC/Rscripts/plotBigPictureHeatmap.R'
RSCRIPT_FOR_GCFGENEPLOTS = lsaBGC_main_directory + '/lsaBGC/Rscripts/gcfGenePlots.R'
RSCRIPT_FOR_DIVERGENCEPLOT = lsaBGC_main_directory + '/lsaBGC/Rscripts/divergencePlot.R'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-AutoAnalyze.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Program to parallelize most of lsaBGC programs for each GCF. 
	
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-o', '--output_directory', help="Parent output/workspace directory.", required=True)
	parser.add_argument('-l', '--input_listing', help="Path to tab delimited file listing: (1) sample name (2) path to Prokka Genbank and (3) path to Prokka predicted proteome. This file is produced by lsaBGC-Process.py.", required=False, default=None)
	parser.add_argument('-g', '--gcf_listing_dir', help='Directory with GCF listing files.', required=True)
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder homolog group by sample matrix.", required=True)
	parser.add_argument('-k', '--sample_set', help="Sample set to keep in analysis. Should be file with one sample id per line.", required=False)
	parser.add_argument('-s', '--lineage_phylogeny', help="Path to species phylogeny. If not provided a MASH based neighborjoining tree will be constructed and used.", default=None, required=False)
	parser.add_argument('-p', '--population_analysis',	action='store_true', help="Whether to construct species phylogeny and use it to determine populations.", default=False, required=False)
	parser.add_argument('-ps', '--num_populations', type=int, help='If population analysis specified, what is the number of populations to . Use the script determinePopulationK.py to see how populations will look with k set to different values.', required=False, default=4)
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

	sample_set_file = myargs.sample_set
	lineage_phylogeny_file = myargs.lineage_phylogeny
	population_analysis = myargs.population_analysis
	num_populations = myargs.num_populations
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
						sample_set_file, num_populations, cores]
	parameter_names = ["GCF Listings Directory", "Listing File of Prokka Annotation Files for Initial Set of Samples",
					   "OrthoFinder Homolog Matrix", "Output Directory", "Phylogeny File in Newick Format",
					   "Delineate Populations and Perform Population Genetics Analytics", "DiscoVary Analysis ID",
					   "DiscoVary Sequencing Data Location Specification File", "Sample Retention Set",
					   "Number of Populations", "Cores"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 0: (Optional) Parse sample set retention specifications file, if provided by the user.
	sample_retention_set = util.getSampleRetentionSet(sample_set_file)

	Pan_Object = Pan(input_listing_file, logObject=logObject)
	logObject.info("Converting Genbanks from Expansion listing file into FASTA per sample.")
	gw_fasta_dir = outdir + 'Sample_Assemblies_in_FASTA/'
	if not os.path.isdir(gw_fasta_dir): os.system('mkdir %s' % gw_fasta_dir)
	gw_fasta_listing_file = outdir + 'Genome_FASTA_Listings.txt'
	Pan_Object.convertGenbanksIntoFastas(gw_fasta_dir, gw_fasta_listing_file)
	logObject.info("Successfully performed conversions.")

	logObject.info("Running MASH Analysis Between Genomes.")
	gw_pairwise_differences = util.calculateMashPairwiseDifferences(gw_fasta_listing_file, outdir, 'genome_wide', 10000, cores, logObject, prune_set=sample_retention_set)
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

	lineage_phylogeny_from_mash = False
	if not lineage_phylogeny_file:
		# Run MASH Analysis Between Genomic Assemblies
		logObject.info("Using MASH estimated distances between genomes to infer neighbor-joining tree.")
		mash_nj_tree = outdir + 'MASH_NeighborJoining_Tree.nwk'

		# create neighbor-joining tree
		cmd = ['Rscript', RSCRIPT_FOR_NJTREECONSTRUCTION, mash_matrix_file, mash_nj_tree]
		try:
			util.run_cmd(cmd, logObject)
			assert(util.is_newick(mash_nj_tree))
			lineage_phylogeny_file = mash_nj_tree
			lineage_phylogeny_from_mash = True
		except Exception as e:
			logObject.error("Had issues with creating neighbor joining tree and defining populations using treestructure.")
			raise RuntimeError("Had issues with creating neighbor joining tree and defining populations using treestructure.")

	elif lineage_phylogeny_file and sample_retention_set != None:
		# Pruning lineage phylogeny provided
		logObject.info("Pruning lineage phylogeny to retain only samples of interest.")
		update_lineage_phylogeny_file = outdir + 'Lineage_Phylogeny.Pruned.nwk'

		t = Tree(lineage_phylogeny_file)
		R = t.get_midpoint_outgroup()
		t.set_outgroup(R)
		t.prune(sample_retention_set)
		#t.resolve_polytomy(recursive=True)
		t.write(format=5, outfile=update_lineage_phylogeny_file)
		lineage_phylogeny_file = update_lineage_phylogeny_file
		logObject.info("Successfully refined lineage phylogeny for sample set of interest.")

	if sample_retention_set != None:
		# Editing GCF listing files and input listing file to keep only samples of interest
		update_input_listing_file = outdir + 'Sample_Annotation_Files.Pruned.txt'
		update_input_listing_handle = open(update_input_listing_file, 'w')
		with open(input_listing_file) as oilf:
			for line in oilf:
				line = line.strip()
				ls = line.split('\t')
				if ls[0] in sample_retention_set:
					update_input_listing_handle.write(line + '\n')
		update_input_listing_handle.close()

		update_gcf_listing_dir = outdir + 'GCF_Listings_Pruned/'
		if not os.path.isdir(update_gcf_listing_dir):
			os.system('mkdir %s' % update_gcf_listing_dir)

		for g in os.listdir(gcf_listing_dir):
			update_gcf_listing_file = update_gcf_listing_dir + g
			update_gcf_listing_handle = open(update_gcf_listing_file, 'w')
			with open(gcf_listing_dir + g) as ogf:
				for line in ogf:
					line = line.strip()
					ls = line.split('\t')
					if ls[0] in sample_retention_set:
						update_gcf_listing_handle.write(line + '\n')
			update_gcf_listing_handle.close()

		input_listing_file = update_input_listing_file
		gcf_listing_dir = update_gcf_listing_dir

	all_samples = set([])
	with open(input_listing_file) as oilf:
		for line in oilf:
			all_samples.add(line.split('\t')[0])

	population_listing_file = None
	if population_analysis:
		population_listing_file = outdir + 'Populations_Defined.txt'
		populations_on_nj_tree_pdf = outdir + 'Populations_on_Lineage_Tree.pdf'
		# create neighbor-joining tree

		cmd = ['Rscript', RSCRIPT_FOR_DEFINECLADES_FROM_PHYLO, lineage_phylogeny_file, str(num_populations), population_listing_file, populations_on_nj_tree_pdf]
		try:
			util.run_cmd(cmd, logObject)
		except Exception as e:
			logObject.error(
				"Had issues with creating neighbor joining tree and defining populations using cutree.")
			raise RuntimeError(
				"Had issues with creating neighbor joining tree and defining populations using cutree.")

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
				   '-i', gcf_id, '-c', str(cores), '-pi', gw_fasta_listing_file, '-pr', outdir + 'genome_wide.out']
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
				   '-i', gcf_id, '-a',	gcf_pop_outdir + 'Codon_Alignments_Listings.txt', '-c', str(cores),
				   '-pi', gw_fasta_listing_file, '-pr', outdir + 'genome_wide.out']
			try:
				util.run_cmd(cmd, logObject)
			except Exception as e:
				logObject.warning("lsaBGC-Divergence.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-Divergence.py was unsuccessful for GCF %s\n" % gcf_id)

		# 4. Run lsaBGC-DiscoVary.py
		if discovary_analysis_id and discovary_input_listing:
			gcf_dis_outdir = dis_outdir + gcf_id + '/'
			if True: #not os.path.isdir(gcf_dis_outdir):
				os.system('mkdir %s' % gcf_dis_outdir)
				cmd = ['lsaBGC-DiscoVary.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-o',
					   gcf_dis_outdir, '-i', gcf_id, '-c', str(cores), '-p', discovary_input_listing, '-a',
					   gcf_pop_outdir + 'Codon_Alignments_Listings.txt', '-l', input_listing_file]
				try:
					util.run_cmd(cmd, logObject, stderr=sys.stderr, stdout=sys.stdout)
				except Exception as e:
					logObject.warning("lsaBGC-DiscoVary.py was unsuccessful for GCF %s" % gcf_id)
					sys.stderr.write("Warning: lsaBGC-DiscoVary.py was unsuccessful for GCF %s\n" % gcf_id)

	combined_gene_plotting_input_file = outdir + 'GCF_Gene_Plotting_Input.txt'
	combined_consensus_similarity_file = outdir + 'GCF_Ortholog_Group_Consensus_Sequence_Similarity.txt'
	combined_divergence_results_file = outdir + 'GCF_Divergences.txt'

	combined_gene_plotting_input_handle = open(combined_gene_plotting_input_file, 'w')
	combined_consensus_similarity_handle = open(combined_consensus_similarity_file, 'w')
	combined_orthoresults_unrefined_handle = open(outdir + 'GCF_Ortholog_Group_Information.txt', 'w')
	combined_divergence_results_handle = open(combined_divergence_results_file, 'w')

	combined_consensus_similarity_handle.write('\t'.join(['GCF', 'GCF_Order', 'Homolog_Group', 'Homolog_Group_Order', 'label', 'Difference_to_Consensus_Sequence']) + '\n')
	for i, g in enumerate(os.listdir(gcf_listing_dir)):
		gcf_id = g.split('.txt')[0]
		gcf_pop_outdir = pop_outdir + gcf_id + '/'
		gcf_div_outdir = div_outdir + gcf_id + '/'

		gcf_pop_results_unrefined = gcf_pop_outdir + 'Ortholog_Group_Information.txt'
		gcf_div_results = gcf_div_outdir + 'Relative_Divergence_Report.txt'

		include_header = False
		if i == 0: include_header = True
		if os.path.isfile(gcf_pop_results_unrefined): writeToOpenHandle(gcf_pop_results_unrefined, combined_orthoresults_unrefined_handle, include_header)
		if os.path.isfile(gcf_div_results): writeToOpenHandle(gcf_div_results, combined_divergence_results_handle, include_header)

		data = []
		with open(gcf_pop_results_unrefined) as ogpru:
			for j, line in enumerate(ogpru):
				line = line.strip('\n')
				ls = line.split('\t')
				if j == 0 and not include_header: continue
				elif j == 0 and include_header:
					if population_analysis:
						combined_gene_plotting_input_handle.write('\t'.join(ls[:2] + ls[3:6] + ['gene_start', 'gene_stop'] + ls[6:-5]) + '\n')
					else:
						combined_gene_plotting_input_handle.write('\t'.join(ls[:2] + ls[3:6] + ['gene_start', 'gene_stop'] + ls[6:-1]) + '\n')
				elif ls[3] != 'NA':
					data.append([int(ls[3]), ls])

		previous_end = 1
		for tupls in sorted(data, key=itemgetter(0)):
			ls = tupls[1]
			if population_analysis:
				combined_gene_plotting_input_handle.write('\t'.join(ls[:2] + ls[3:6] + [str(previous_end), str(previous_end + int(float(ls[6])))] + ls[6:15] + [ls[15].split(' [')[0].replace('Conserved', 'NA').replace('Infinite', 'NA').strip()] + ls[16:-5]) + '\n')
			else:
				combined_gene_plotting_input_handle.write('\t'.join(ls[:2] + ls[3:6] + [str(previous_end), str(previous_end + int(float(ls[6])))] + ls[6:-2] + [ls[-2].split(' [')[0].replace('Conserved', 'NA').replace('Infinite', 'NA').strip()]) + '\n')
			previous_end = previous_end + int(float(ls[6])) + 1

		hg_ordering = defaultdict(lambda: 'NA')
		with open(gcf_pop_results_unrefined) as ogpru:
			for i, line in enumerate(ogpru):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				hg_ordering[ls[1]] = ls[3]

		gcf_pop_stats_outdir = gcf_pop_outdir + 'Codon_PopGen_Analyses/'
		for f in os.listdir(gcf_pop_stats_outdir):
			if f.endswith("_sim_to_consensus.txt"):
				samps_accounted = set([])
				hg = None
				with open(gcf_pop_stats_outdir + f) as ogpso:
					for line in ogpso:
						line = line.strip()
						hg, samp, diff = line.split('\t')
						samps_accounted.add(samp)
						combined_consensus_similarity_handle.write(gcf_id + '\t' + gcf_id.split('_')[1] + '\t' + hg + '\t' + hg_ordering[hg] + '\t' + samp + '\t' + str(diff) + '\n')
				for samp in all_samples:
					if not samp in samps_accounted:
						combined_consensus_similarity_handle.write(gcf_id + '\t' + gcf_id.split('_')[1] + '\t' + hg + '\t' + hg_ordering[hg] + '\t' + samp + '\t' + str(1.0) + '\n')

	combined_gene_plotting_input_handle.close()
	combined_consensus_similarity_handle.close()
	combined_orthoresults_unrefined_handle.close()
	combined_divergence_results_handle.close()

	# Create Final R plots:

	# create big-picture heatmap of presence/sequence-similarity to consensus sequence of homolog groups from each gcf
	big_picture_heatmap_pdf_file = outdir + 'Consensus_Sequence_Similarity_of_Homolog_Groups.pdf'
	cmd = ['Rscript', RSCRIPT_FOR_BIGPICTUREHEATMAP, lineage_phylogeny_file, combined_consensus_similarity_file,
		   population_listing_file, big_picture_heatmap_pdf_file]
	try:
		util.run_cmd(cmd, logObject)
	except Exception as e:
		logObject.error("Had issues with creating big picture heatmap.")
		raise RuntimeError("Had issues with creating big picture heatmap.")

	# create gcf pop gen stats and conservation plots
	gcf_gene_views_pdf_file = outdir + 'GCF_Conservation_and_PopStats_Views.pdf'
	cmd = ['Rscript', RSCRIPT_FOR_GCFGENEPLOTS, combined_gene_plotting_input_file, gcf_gene_views_pdf_file]
	try:
		util.run_cmd(cmd, logObject)
	except Exception as e:
		logObject.error("Had issues with creating GCF gene conservation and population stats views.")
		raise RuntimeError("Had issues with creating GCF gene conservation and population stats views.")

	# create divergence plot
	gcf_divergence_pdf_file = outdir + 'GCF_Divergence.pdf'
	cmd = ['Rscript', RSCRIPT_FOR_DIVERGENCEPLOT, combined_divergence_results_file, gcf_divergence_pdf_file]
	try:
		util.run_cmd(cmd, logObject)
	except Exception as e:
		logObject.error("Had issues with creating GCF vs. genome-wide divergence plots.")
		raise RuntimeError("Had issues with creating GCF vs. genome-wide divergence plots.")


	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	lsaBGC_AutoAnalyze()
