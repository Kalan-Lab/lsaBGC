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
import traceback
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
RSCRIPT_FOR_POPULATIONS_ON_PHYLO = lsaBGC_main_directory + '/lsaBGC/Rscripts/GeneRatePhylogeny.R'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-AutoAnalyze.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Wrapper program to automate running lsaBGC analytical programs for each GCF. 
	
	Iteratively runs lsaBGC-See.py, lsaBGC-PopGene.py, lsaBGC-Divergence.py, and optionally lsaBGC-DiscoVary.py for
	each GCF in a GCF listings directory, produced by lsaBGC-Ready, lsaBGC-Cluster, or lsaBGC-AutoExpansion.py.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-o', '--output_directory', help="Parent output/workspace directory.", required=True)
	parser.add_argument('-l', '--input_listing', help="Path to tab delimited file listing: (1) sample name\n(2) path to whole-genome Genbank and (3) path to whole-genome predicted proteome\n(an output of lsaBGC-Ready.py or lsaBGC-AutoExpansion.py).", required=False, default=None)
	parser.add_argument('-g', '--gcf_listing_dir', help='Directory with GCF listing files.', required=True)
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder homolog group by sample matrix.", required=True)
	parser.add_argument('-k', '--sample_set', help="Sample set to keep in analysis. Should be file with one sample id per line.", required=False)
	parser.add_argument('-s', '--species_phylogeny', help="Path to species phylogeny. If not provided a FastANI based neighborjoining tree will be constructed and used.", default=None, required=False)
	parser.add_argument('-w', '--genome_wide_distances', help="Path to file listing genome-wide distances between genomes/samples. This is the Genome-Wide_Estimates.txt file produced by the computeGenomeWideEstimates.py script")
	parser.add_argument('-r', '--aai', action='store_true', help='AAI was used to compute genome wise distances instead of ANI. E.g. if CompareM was used.')
	parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is antiSMASH.', default='antiSMASH', required=False)
	parser.add_argument('-u', '--populations', help='Path to user defined populations/groupings file. Tab delimited with 2 columns: (1) sample name and (2) group identifier.', required=False, default=None)
	parser.add_argument('-i', '--discovary_input_listing', help="Sequencing readsets for DiscoVary analysis. Tab delimited file listing: (1) sample name, (2) forward readset, (3) reverse readset for metagenomic/isolate sequencing data.", required=False, default=None)
	parser.add_argument('-n', '--discovary_analysis_name', help="Identifier/name for DiscoVary. Not providing this parameter will avoid running lsaBGC-DiscoVary step.", required=False, default=None)
	parser.add_argument('-c', '--cores', type=int, help="Total number of cores to use.", required=False, default=1)

	args = parser.parse_args()
	return args

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
		raise RuntimeError('Input directory with GCF listings does not exist, sample annotation file, or the OrthoFinder matrix does not exist. Exiting now ...')

	all_samples = set([])
	try:
		with open(input_listing_file) as oilf:
			for line in oilf:
				line = line.strip()
				sample, sample_gbk, sample_faa = line.split('\t')
				assert(util.is_genbank(sample_gbk))
				assert(util.is_fasta(sample_faa))
				all_samples.add(line.split('\t')[0])
	except:
		raise RuntimeError("Issue with validating input listing file with sample genome-wide genbanks and predicted proteomes is in correct format.")

	try:
		gcf_listings = 0
		for f in os.listdir(gcf_listing_dir):
			if f.startswith('GCF_') and f.endswith('.txt'):
				with open(gcf_listing_dir + f) as ogf:
					for line in ogf:
						line = line.strip()
						samp, bgc_path = line.split('\t')
						assert(samp in all_samples)
						assert(util.is_genbank(bgc_path))
				gcf_listings += 1
		assert(gcf_listings > 0)
	except:
		raise RuntimeError("Issue with validating at least single GCF listing file in the GCF listings directory.")

	try:
		with open(original_orthofinder_matrix_file) as omf:
			for i, line in enumerate(omf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					for samp in ls[1:]:
						if samp != '':
							assert(samp in all_samples)
				else:
					assert(ls[0].startswith("OG"))
	except:
		raise RuntimeError("Issue with validating OrthoFinder matrix is in correct format.")


	"""
	PARSE OPTIONAL INPUTS
	"""

	sample_set_file = myargs.sample_set
	bgc_prediction_software = myargs.bgc_prediction_software.upper()
	species_phylogeny_file = myargs.species_phylogeny
	genomewide_distances_file = myargs.genome_wide_distances
	aai_flag = myargs.aai
	population_listing_file = myargs.populations
	discovary_analysis_id = myargs.discovary_analysis_name
	discovary_input_listing = myargs.discovary_input_listing
	cores = myargs.cores

	try:
		assert (bgc_prediction_software in set('ANTISMASH', 'DEEPBGC', 'GECCO'))
	except:
		raise RuntimeError('BGC prediction software option is not a valid option.')

	if species_phylogeny_file != None:
		try:
			assert(util.is_newick(species_phylogeny_file))
		except:
			raise RuntimeError('Species phylogeny provided either does not exist or is not in the proper format. Exiting now ...')

	if genomewide_distances_file != None:
		try:
			assert(os.path.isfile(genomewide_distances_file))
			with open(genomewide_distances_file) as ogdf:
				for line in ogdf:
					line = line.strip()
					sample_1, sample_2, gen_est = line.split('\t')
					assert(sample_1 in all_samples and sample_2 in all_samples and util.is_numeric(gen_est))
		except:
			raise RuntimeError('GenomeWide estimates file provided does not exist or does not meet expected formatting. Exiting now ...')

	if population_listing_file != None:
		try:
			with open(population_listing_file) as oplf:
				for line in oplf:
					line = line.strip()
					sample, pop = line.split('\t')
					assert(sample in all_samples)
		except:
			raise RuntimeError('Population listings file provided does meet expected formatting.')

	if discovary_input_listing != None:
		try:
			with open(discovary_input_listing) as odilf:
				for line in odilf:
					line = line.strip()
					sample, frw_pe_file, rev_pe_file = line.split('\t')
					assert(os.path.isfastq(frw_pe_file))
					assert(os.path.isfastq(rev_pe_file))
					assert(sample in all_samples)
		except:
			raise RuntimeError('DiscoVary paired-end sequencing data listings file does not meet expected formatting.')

	"""
	START WORKFLOW
	"""

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Continuing in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for easier determination of results' provenance in the future.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [gcf_listing_dir, input_listing_file, original_orthofinder_matrix_file, outdir,
						species_phylogeny_file, genomewide_distances_file, aai_flag, population_listing_file,
						discovary_analysis_id, discovary_input_listing, bgc_prediction_software, sample_set_file, cores]
	parameter_names = ["GCF Listings Directory", "Listing File of Sample Annotation Files for Initial Set of Samples",
					   "OrthoFinder Homolog Matrix", "Output Directory", "Species Phylogeny File in Newick Format",
					   "File with GenomeWide Distance Estimations", "CompareM AAI Was Used for GenomeWide Distance Estimations?",
					   "Clade/Population Listings File", "DiscoVary Analysis ID",
					   "DiscoVary Sequencing Data Location Specification File", "BGC Prediction Software",
					   "Sample Retention Set", "Cores"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 0: (Optional) Parse sample set retention specifications file, if provided by the user.
	sample_retention_set = util.getSampleRetentionSet(sample_set_file)

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
					sample, bgc_path = line.split('\t')
					if sample in sample_retention_set:
						update_gcf_listing_handle.write(line + '\n')

			update_gcf_listing_handle.close()

		input_listing_file = update_input_listing_file
		gcf_listing_dir = update_gcf_listing_dir

		if genomewide_distances_file != None:
			update_genomewide_distances_file = outdir + 'GenomeWide_Distance_Estimates.Pruned.txt'
			update_genomewide_distances_handle = open(update_genomewide_distances_file, 'w')
			with open(genomewide_distances_file) as ogdf:
				for line in ogdf:
					line = line.strip()
					sample_1, sample_2, gw_est = line.split('\t')
					if sample_1 in sample_retention_set and sample_2 in sample_retention_set:
						update_genomewide_distances_handle.write(line + '\n')
			update_genomewide_distances_handle.close()
			genomewide_distances_file = update_genomewide_distances_file

	if population_listing_file and species_phylogeny_file:
		populations_on_tree_pdf = outdir + 'Populations_on_Species_Tree.pdf'
		cmd = ['Rscript', RSCRIPT_FOR_POPULATIONS_ON_PHYLO, species_phylogeny_file, population_listing_file, populations_on_tree_pdf]
		try:
			util.run_cmd(cmd, logObject)
		except Exception as e:
			logObject.error("Had issues showcasing manually defined population labels on phylogeny/tree.")
			raise RuntimeError("Had issues showcasing manually defined population labels on phylogeny/tree.")

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
		lsabgc_see_checkpoint = gcf_see_outdir + 'CHECKPOINT.txt'
		if not os.path.isfile(lsabgc_see_checkpoint):
			os.system('rm -rf %s' % gcf_see_outdir)
			os.system('mkdir %s' % gcf_see_outdir)
			cmd = ['lsaBGC-See.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-o', gcf_see_outdir,
				   '-i', gcf_id, '-s', species_phylogeny_file, '-y', '-p', bgc_prediction_software, '-c', str(cores)]
			try:
				util.run_cmd(cmd, logObject)
				assert(os.path.isfile())
			except Exception as e:
				logObject.warning("lsaBGC-See.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-See.py was unsuccessful for GCF %s\n" % gcf_id)

		# 2. Run lsaBGC-PopGene.py
		gcf_pop_outdir = pop_outdir + gcf_id + '/'
		lsabgc_popgene_checkpoint = gcf_pop_outdir + 'CHECKPOINT.txt'
		if not os.path.isfile(lsabgc_popgene_checkpoint):
			os.system('rm -rf %s' % gcf_pop_outdir)
			os.system('mkdir %s' % gcf_pop_outdir)
			cmd = ['lsaBGC-PopGene.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-o', gcf_pop_outdir,
				   '-i', gcf_id, '-p', bgc_prediction_software, '-c', str(cores)]
			if genomewide_distances_file != None:
				cmd += ['-f', genomewide_distances_file]
			if population_listing_file != None:
				cmd += ['-u', population_listing_file]
			if aai_flag:
				cmd += ['-cm']
			try:
				util.run_cmd(cmd, logObject, stderr=sys.stderr, stdout=sys.stdout)
			except Exception as e:
				logObject.warning("lsaBGC-PopGene.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-PopGene.py was unsuccessful for GCF %s\n" % gcf_id)

		# 3. Run lsaBGC-Divergence.py
		gcf_div_outdir = div_outdir + gcf_id + '/'
		lsabgc_divergence_checkpoint = gcf_div_outdir + 'CHECKPOINT.txt'
		if not os.path.isfile(lsabgc_divergence_checkpoint):
			os.system('rm -rf %s' % gcf_div_outdir)
			os.system('mkdir %s' % gcf_div_outdir)
			cmd = ['lsaBGC-Divergence.py', '-g', gcf_listing_file, '-l', input_listing_file, '-o', gcf_div_outdir,
				   '-i', gcf_id, '-a',	gcf_pop_outdir + 'Codon_Alignments_Listings.txt', '-c', str(cores)]
			if genomewide_distances_file != None:
				cmd += ['-f', genomewide_distances_file]
			if aai_flag:
				cmd += ['-cm']
			try:
				util.run_cmd(cmd, logObject)
			except Exception as e:
				logObject.warning("lsaBGC-Divergence.py was unsuccessful for GCF %s" % gcf_id)
				sys.stderr.write("Warning: lsaBGC-Divergence.py was unsuccessful for GCF %s\n" % gcf_id)

		# 4. Run lsaBGC-DiscoVary.py
		if discovary_analysis_id and discovary_input_listing:
			gcf_dis_outdir = dis_outdir + gcf_id + '/'
			lsabgc_discovary_checkpoint = gcf_dis_outdir + 'CHECKPOINT.txt'
			if not os.path.isfile(lsabgc_discovary_checkpoint):
				os.system('rm -rf %s' % gcf_dis_outdir)
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
	combined_consensus_similarity_file = outdir + 'GCF_Homolog_Group_Consensus_Sequence_Similarity.txt'
	combined_divergence_results_file = outdir + 'GCF_Divergences.txt'

	combined_gene_plotting_input_handle = open(combined_gene_plotting_input_file, 'w')
	combined_consensus_similarity_handle = open(combined_consensus_similarity_file, 'w')
	combined_orthoresults_unrefined_handle = open(outdir + 'GCF_Homolog_Group_Information.txt', 'w')
	combined_divergence_results_handle = open(combined_divergence_results_file, 'w')

	if sample_retention_set == None:
		sample_retention_set = all_samples

	combined_consensus_similarity_handle.write('\t'.join(['GCF', 'GCF_Order', 'Homolog_Group', 'Homolog_Group_Order', 'label', 'Difference_to_Consensus_Sequence']) + '\n')
	for i, g in enumerate(os.listdir(gcf_listing_dir)):
		gcf_id = g.split('.txt')[0]
		gcf_pop_outdir = pop_outdir + gcf_id + '/'
		gcf_div_outdir = div_outdir + gcf_id + '/'

		gcf_pop_results_unrefined = gcf_pop_outdir + 'Homolog_Group_Information.txt'
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
					if population_listing_file != None:
						combined_gene_plotting_input_handle.write('\t'.join(ls[:3] + ls[4:7] + ['gene_start', 'gene_stop'] + ls[7:-2]) + '\n')
					else:
						combined_gene_plotting_input_handle.write('\t'.join(ls[:3] + ls[4:7] + ['gene_start', 'gene_stop'] + ls[7:-1]) + '\n')
				elif ls[4] != 'NA':
					data.append([int(ls[4]), ls])

		previous_end = 1
		for tupls in sorted(data, key=itemgetter(0)):
			ls = tupls[1]
			if population_listing_file != None:
				combined_gene_plotting_input_handle.write('\t'.join(ls[:3] + ls[4:7] + [str(previous_end), str(previous_end + int(float(ls[7])))] + ls[7:-2]) + '\n')
			else:
				combined_gene_plotting_input_handle.write('\t'.join(ls[:3] + ls[4:7] + [str(previous_end), str(previous_end + int(float(ls[7])))] + ls[7:-1]) + '\n')
			previous_end = previous_end + int(float(ls[7])) + 1

		hg_ordering = defaultdict(lambda: 'NA')
		with open(gcf_pop_results_unrefined) as ogpru:
			for i, line in enumerate(ogpru):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				hg_ordering[ls[2]] = ls[4]

		gcf_pop_stats_outdir = gcf_pop_outdir + 'Codon_PopGen_Analyses/'
		gcf_consensus_sim_plot_lines = []
		total_samples_accounted = set([])
		for f in os.listdir(gcf_pop_stats_outdir):
			if f.endswith("_sim_to_consensus.txt"):
				samps_accounted = set([])
				hg = None
				with open(gcf_pop_stats_outdir + f) as ogpso:
					for line in ogpso:
						line = line.strip()
						hg, samp, diff = line.split('\t')
						samps_accounted.add(samp)
						total_samples_accounted.add(samp)
						gcf_consensus_sim_plot_lines.append(gcf_id + '\t' + gcf_id.split('_')[1] + '\t' + hg + '\t' + hg_ordering[hg] + '\t' + samp + '\t' + str(diff))
				for samp in all_samples:
					if not samp in samps_accounted:
						gcf_consensus_sim_plot_lines.append(gcf_id + '\t' + gcf_id.split('_')[1] + '\t' + hg + '\t' + hg_ordering[hg] + '\t' + samp + '\t' + "NA")

		if float(len(total_samples_accounted))/float(len(sample_retention_set)) >= 0.1:
			combined_consensus_similarity_handle.write('\n'.join(gcf_consensus_sim_plot_lines) + '\n')
	combined_gene_plotting_input_handle.close()
	combined_consensus_similarity_handle.close()
	combined_orthoresults_unrefined_handle.close()
	combined_divergence_results_handle.close()

	# Create Final R plots:

	# create big-picture heatmap of presence/sequence-similarity to consensus sequence of homolog groups from each gcf
	big_picture_heatmap_pdf_file = outdir + 'Consensus_Sequence_Similarity_of_Homolog_Groups.pdf'
	cmd = ['Rscript', RSCRIPT_FOR_BIGPICTUREHEATMAP, species_phylogeny_file, combined_consensus_similarity_file,
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

def writeToOpenHandle(gcf_results_file, combined_results_handle, include_header):
	with open(gcf_results_file) as ogrf:
		for i, line in enumerate(ogrf):
			if i == 0 and include_header:
				combined_results_handle.write(line)
			elif i != 0:
				combined_results_handle.write(line)

if __name__ == '__main__':
	lsaBGC_AutoAnalyze()