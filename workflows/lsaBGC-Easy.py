#!/usr/bin/env python

### Program: lsaBGC-Easy.py
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
from lsaBGC import util
import subprocess
import traceback
import multiprocessing
import math
import itertools
from operator import itemgetter
from collections import defaultdict
os.environ['OMP_NUM_THREADS'] = '4'
lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'


def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Easy.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	
	*******************************************************************************************************************
	CONSIDERATIONS:	
	* If you have the computational resources available, please consider running with OrthoFinder mode set to 
	Genome_Wide (not the default). Also consider using antiSMASH, as these were the methods we used when 
	benchmarking lsaBGC in our manuscript.

	* Check out the considerations wiki page: https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background-&-Considerations  	

	*******************************************************************************************************************
	DESCRIPTION:	
	Downloads genomes for a taxa from NCBI using ncbi-genome-download, performs BGC predictions 
	using GECCO (very-light weight dependency, see tutorial on Wiki on how to use antiSMASH or DeepBGC instead), and 
	then runs the full lsaBGC suite to generate a "quick" (not guaranteed, depends on the number of cores you have) and 
	"easy" (I think its pretty easy - but obviously there are considerations here and evolutionary statistics will need 
	to be interpreted with caution depending on the population structure of the dataset!). Automatic dereplication
	is attempted at 0.999 amino acid identity along single-copy core genes used for phylogeny construction by GToTree.
	
	Specifically all genomes for the taxa from Genbank will be downloaded, then use Treemmer to select 30 distinctive 
	representatives (if > 30 genomes exist). GECCO will then be run on these genomes and the BGCs predicted within these
	genomes will be grouped together using lsaBGC-Cluster.py into GCFs. GCFs will be searched in additional genomes 
	(not one of the 30 reps) using lsaBGC-AutoExpansion.py. Instead of GECCO you could specify antiSMASH or DeepBGC,
	however this won't be automatic because they are not part of the lsaBGC conda environment. You will need to run the 
	command task file produced by lsaBGC-Easy.py separately and then once they have finished restart lsaBGC-Easy.py with
	the same command. It should pick up the new results and the completed steps and continue onward. Besides DeepBGC,
	GECCO and antiSMASH predictions will be attempted to be filtered to retain only complete BGC instances in primary
	genomes to improve clustering quality at the expense of sensitivity. For better clustering of fragmented BGCs we
	recommend using lsaBGC programs independently and incorporating BiG-SCAPE clustering - which has a local synteny 
	similarity measure. 

	Currently only works for bacteria, but lsaBGC can handle fungi now - just not lsaBGC-Easy ... yet

	*******************************************************************************************************************
	OVERVIEW OF STEPS TAKEN:	
		- Check number of genomes for taxa is not too crazy (<500) and greater than 10
		- Step 1: Download all genomes in Genbank format and extract proteins to folder
		- Step 2: Run MASH to Quickly Assess Sample Similarity and Select Primary Genomes Through a Rough Dereplication
		- Step 3: Create primary and additional genomes listing
		- Step 4: Run GECCO based annotation of BGCs and crete BGC listing
		- Step 5: Run lsaBGC-Ready.py with lsaBGC-Cluster.py & lsaBGC-Expansion.py
		- Step 6: Dereplicate and Group Samples into Populations/Clades
		- Step 7: Run lsaBGC-AutoAnalyze.py

	*******************************************************************************************************************
	AVAILBLE SCG MODELS IN GTOTREE:
           Actinobacteria                    (138 genes)
           Alphaproteobacteria               (117 genes)
           Archaea                            (76 genes)
           Bacteria                           (74 genes)
           Bacteria_and_Archaea               (25 genes)
           Bacteroidetes                      (90 genes)
           Betaproteobacteria                (203 genes)
           Chlamydiae                        (286 genes)
           Cyanobacteria                     (251 genes)
           Epsilonproteobacteria             (260 genes)
           Firmicutes                        (119 genes)
           Gammaproteobacteria               (172 genes)
           Proteobacteria                    (119 genes)
           Tenericutes                        (99 genes)
           Universal_Hug_et_al                (16 genes)
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-n', '--taxa_name', help='Name of the taxa of interest. Currently restricted to bacteria - but lsaBGC\nitself- when using antiSMASH BGC predictions can run on\nEuks - check out the tutorial!! If there is a space in the\nname, please surround by quotes.', required=True)
	parser.add_argument('-o', '--output_directory', help='Parent output/workspace directory.', required=True)
	parser.add_argument('-g', '--user_genomes_directory', help='A directory with additional genomes, e.g. those recently sequenced by the\nuser, belonging to the taxa. Accepted formats include FASTA.\nAccepted suffices include: .fna, .fa, .fasta.', required=False, default=None)
	parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is GECCO - which will be automatic. For the other two,\nyou will need to install them and run manually and then rerun lsaBGC-Easy.py with\nthe same command.', default='gecco', required=False)
	parser.add_argument('-m', '--max_primary_genomes', type=int, help="Maximum number of primary genomes to use for initial BGC finding. Default is 30.", default=30, required=False)
	parser.add_argument('-om', '--orthofinder_mode', help="Method for running OrthoFinder2. (Options: Genome_Wide, BGC_Only).\nDefault is BGC_Only (slightly experimental but much faster and should work\nespecially well for taxa with many BGCs).", default='BGC_Only', required=False)
	parser.add_argument('-sd', '--skip_dereplication', action='store_true', help="Whether to skip dereplication based on secondary GToTree analysis - not recommended.", default=False, required=False)
	parser.add_argument('-gtm', '--gtotree_model', help="SCG model for secondary GToTree analysis and what would be used for dereplication.", default='Bacteria', required=False)
	parser.add_argument('-iib', '--include_incomplete_bgcs', help="Whether to account for incomplete BGCs in primary genomes prior to clustering - not recommended.", default=False, required=False)
	parser.add_argument('-c', '--cores', type=int, help="Total number of cores/threads. Note, this is the total number of\nthreads to use. BGC prediction commands (via GECCO, antiSMASH, and DeepBGC) will\neach be set to use 4 cores by default.", required=False, default=4)
	parser.add_argument('-dt', '--dereplicate_threshold', type=float, help="Amino acid similarity threshold of SCGs for considering\ntwo genomes as redundant.", default=0.999, required=False)
	parser.add_argument('-pt', '--population_threshold', type=float, help="Amino acid similarity threshold of SCGs for considering\ntwo genomes as belonging to the same population.", default=0.99, required=False)
	args = parser.parse_args()
	return args

def lsaBGC_Easy():
	myargs = create_parser()

	taxa_name = myargs.taxa_name.strip('"')
	outdir = os.path.abspath(myargs.output_directory) + '/'
	cores = myargs.cores
	max_primary_genomes = myargs.max_primary_genomes
	user_genomes_directory = myargs.user_genomes_directory
	bgc_prediction_software = myargs.bgc_prediction_software.upper()
	gtotree_model = myargs.gtotree_model
	skip_dereplication_flag = myargs.skip_dereplication
	include_incomplete_bgcs_flag = myargs.include_incomplete_bgcs
	dereplicate_threshold = myargs.dereplicate_threshold
	population_threshold = myargs.population_threshold
	orthofinder_mode = myargs.orthofinder_mode.upper()

	try:
		assert (orthofinder_mode in set(['GENOME_WIDE', 'BGC_ONLY']))
	except:
		raise RuntimeError('BGC prediction software option is not a valid option.')

	try:
		assert (bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
	except:
		raise RuntimeError('BGC prediction software option is not a valid option.')

	if user_genomes_directory:
		try:
			user_genomes_directory = os.path.abspath(user_genomes_directory) + '/'
			assert (os.path.isdir(user_genomes_directory))
		except:
			raise RuntimeError('User defined genomes directory could not be validated to exist!')

	checkModelIsValid(gtotree_model)

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory already exists! Overwriting in 5 seconds, but only where needed will use checkpoint files to skip certain steps...\n")
	else:
		util.setupReadyDirectory([outdir])

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)
	logObject.info("Saving parameters for future records.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [taxa_name, outdir, cores]
	parameter_names = ["Taxa", "Output Directory", "Number of Cores"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	parallel_jobs_4cpu = max(math.floor(cores / 4), 1)

	genome_listing_file = outdir + 'NCBI_Genomes_from_Genbank_for_Taxa.txt'
	if not os.path.isfile(genome_listing_file):
		# Check number of genomes for taxa is not too crazy (<500) and greater than 10
		ngd_dry_cmd = ['ncbi-genome-download', '--dry-run', '--section', 'genbank', '--genera',
					   '"' + taxa_name + '"', 'bacteria', '>', genome_listing_file]
		runCmdViaSubprocess(ngd_dry_cmd, logObject, check_files=[genome_listing_file])

	oglf = open(genome_listing_file)
	genome_count = len(oglf.readlines())
	oglf.close()
	if genome_count >= 500 or genome_count < 5:
		logObject.error("Currently not recommended to use lsaBGC-Easy for taxa with > 500 genomes or < 5 genomes. In the future this will likely be updated to be a warning, but would rather not risk harming your server!")
		raise RuntimeError("Currently not recommended to use lsaBGC-Easy for taxa with > 500 genomes or < 5 genomes. In the future this will likely be updated to be a warning, but would rather not risk harming your computer/server!")
	checkUlimitSettings(genome_count, logObject)

	# Step 1: Download all genomes in FASTA formata
	genomes_directory = outdir + 'genbank/'
	all_genomes_listing_file = outdir + 'All_Genomes_Listing.txt'
	uncompress_dir = outdir + 'Uncompressed_Genomes/'
	if not os.path.isfile(all_genomes_listing_file):
		ngd_real_cmd = ['ncbi-genome-download', '--formats', 'fasta', '--retries', '2', '--section',
						'genbank', '--genera', '"' + taxa_name + '"', 'bacteria', '-o', outdir]
		runCmdViaSubprocess(ngd_real_cmd, logObject, check_directories=[genomes_directory])
		list_all_genomes_cmd = ['listAllGenomesInDirectory.py', '-i', genomes_directory, '-u', uncompress_dir, '-z',
								'-d', genome_listing_file, '>', all_genomes_listing_file]
		runCmdViaSubprocess(list_all_genomes_cmd, logObject, check_files=[all_genomes_listing_file])

		if user_genomes_directory:
			list_all_user_genomes_cmd = ['listAllGenomesInDirectory.py', '-i', user_genomes_directory, '>>',
										 all_genomes_listing_file]
			runCmdViaSubprocess(list_all_user_genomes_cmd, logObject)

	# Step 2: Run MASH to Quickly Assess Sample Similarity and Select Primary Genomes Through a Rough Dereplication
	mash_input_file = outdir + 'MASH_Input.txt'
	mash_input_handle = open(mash_input_file, 'w')
	all_genome_listings_fna = {}
	mash_to_genome_name_mapping = {}
	with open(all_genomes_listing_file) as oglf:
		for line in oglf:
			line = line.strip()
			sample, fna_path = line.split('\t')
			fna_to_genome_name_mapping[fna_path] = sample
			all_genome_listings_fna[sample] = fna_path
			mash_input_handle.write()

	# Step 2: Construct GToTree phylogeny
	gtotree_outdir = outdir + 'GToTree_output/'
	guiding_tree_file = gtotree_outdir + 'GToTree_output.tre'
	if not os.path.isfile(guiding_tree_file):
		os.system('rm -rf %s' % gtotree_outdir)
		gtotree_cmd = ['GToTree', '-A', all_proteomes_listing_file, '-H', 'Bacteria', '-n', '4', '-j',
					   str(parallel_jobs_4cpu), '-M', '4', '-o', gtotree_outdir]
		runCmdViaSubprocess(gtotree_cmd, logObject, check_files=[guiding_tree_file])

	# Step 3: Select 30 most diverse genome representatives if needed using Treemmer
	need_for_additional_genomes = True
	primary_samples = set([])
	if genome_count <= max_primary_genomes:
		need_for_additional_genomes = False
		primary_samples = set(all_genome_listings_gbk.keys())
	else:
		nr_listing_file = guiding_tree_file + '_trimmed_list_X_' + str(max_primary_genomes)
		if not os.path.isfile(nr_listing_file):
			treemmer_cmd = ['Treemmer_v0.3.py', guiding_tree_file, '-X', str(max_primary_genomes), '-c', str(cores)]
			runCmdViaSubprocess(treemmer_cmd, logObject)
			with open(nr_listing_file) as onlf:
				for line in onlf:
					sample = line.strip()
					primary_samples.add(sample)

	# Step 4: Create primary and additional genomes listing
	primary_genomes_listing_file = outdir + 'Primary_Genomes.txt'
	primary_bgc_pred_directory = outdir + 'Primary_Genomes_BGC_Predictions/'
	additional_genomes_listing_file = outdir + 'Additional_Genomes.txt'
	if not (os.path.isfile(primary_genomes_listing_file)):
		primary_genomes_listing_handle = open(primary_genomes_listing_file, 'w')
		if need_for_additional_genomes:
			additional_genomes_listing_handle = open(additional_genomes_listing_file, 'w')

		util.setupReadyDirectory([primary_bgc_pred_directory])

		primary_bgc_pred_cmds = []
		for s in all_genome_listings_gbk:
			if s in primary_samples:
				primary_genomes_listing_handle.write(s + '\t' + all_genome_listings_gbk[s] + '\n')
				if bgc_prediction_software == 'GECCO':
					gecco_cmd = ['gecco', 'run', '-j', '4', '-o', primary_bgc_pred_directory + s + '/', '-g',
								 all_genome_listings_fna[s]]
					if not include_incomplete_bgcs_flag:
						gecco_cmd += ['-E', '10']
					gecco_cmd += [logObject]
					primary_bgc_pred_cmds.append(gecco_cmd)
				if bgc_prediction_software == 'ANTISMASH':
					antismash_cmd = ['antismash', '--output-dir', primary_bgc_pred_directory + s + '/', '-c', '4',
									 '--genefinding-tool', 'none', '--output-basename', s, all_genome_listings_gbk[s]]
					primary_bgc_pred_cmds.append(antismash_cmd)
				elif bgc_prediction_software == 'DEEPBGC':
					deepbgc_cmd = ['deepbgc', 'pipeline', '--output', primary_bgc_pred_directory + s + '/',
								   all_genome_listings_fna[s]]
					primary_bgc_pred_cmds.append(deepbgc_cmd)

			elif need_for_additional_genomes:
				additional_genomes_listing_handle.write(s + '\t' + all_genome_listings_gbk[s] + '\n')
		primary_genomes_listing_handle.close()
		if need_for_additional_genomes:
			additional_genomes_listing_handle.close()

		# Step 5: Run GECCO or print commands for others and exit based annotation of BGCs and crete BGC listing
		if bgc_prediction_software == 'GECCO':
			p = multiprocessing.Pool(parallel_jobs_4cpu)
			p.map(util.multiProcess, primary_bgc_pred_cmds)
			p.close()
		else:
			cmd_file = outdir + 'bgc_prediction.cmds'
			cmd_handle = open(cmd_file, 'w')
			for cmd in primary_bgc_pred_cmds:
				cmd_handle.write(' '.join(cmd) + '\n')
			cmd_handle.close()
			print('%s BGC prediction commands written to %s. Please run these and afterwards restart lsaBGC-Easy.py with the same command afterwards.' % (bgc_prediction_software, cmd_file))
			sys.exit(0)

	primary_bgcs_listing_file = outdir + 'BGCs_in_Primary_Genomes.txt'
	if not os.path.isfile(primary_bgcs_listing_file):
		list_primary_bgc_cmd = ['listAllBGCGenbanksInDirectory.py', '-i', primary_bgc_pred_directory, '-f', '-p',
								bgc_prediction_software, '>', primary_bgcs_listing_file]
		runCmdViaSubprocess(list_primary_bgc_cmd, logObject, check_files=[primary_bgcs_listing_file])

	# Step 6: Run lsaBGC-Ready.py with lsaBGC-Cluster.py & lsaBGC-Expansion.py
	lsabgc_ready_directory = outdir + 'lsaBGC_Ready_Results/'
	lsabgc_ready_results_directory = lsabgc_ready_directory + 'Final_Results/'
	annotation_listing_file = lsabgc_ready_results_directory + 'Expanded_Sample_Annotation_Files.txt'
	gcf_listing_dir = lsabgc_ready_results_directory + 'Expanded_GCF_Listings/'
	orthogroups_matrix_file = lsabgc_ready_results_directory + 'Expanded_Orthogroups.tsv'
	species_tree_file = lsabgc_ready_results_directory + 'GToTree_output.tre'
	expected_similarities_file = lsabgc_ready_results_directory + 'GToTree_Expected_Similarities.txt'
	samples_to_keep_file = lsabgc_ready_results_directory + 'Samples_in_GToTree_Tree.txt'
	if not os.path.isdir(lsabgc_ready_results_directory) or len([f for f in os.listdir(lsabgc_ready_results_directory)])<=2:
		os.system('rm -rf %s' % lsabgc_ready_directory)
		lsabgc_ready_cmd = ['lsaBGC-Ready.py', '-i', primary_genomes_listing_file, '-l', primary_bgcs_listing_file,
							'-o', lsabgc_ready_directory, '-lc', '-c', str(cores), '-p', bgc_prediction_software, '-t',
							'-gtm', gtotree_model, '-m', orthofinder_mode]
		db_listing_file = lsaBGC_main_directory + '/db/database_location_paths.txt'
		if os.path.isfile(db_listing_file):
			lsabgc_ready_cmd += ['-a']
		if need_for_additional_genomes:
			lsabgc_ready_cmd += ['-d', additional_genomes_listing_file, '-le']
		runCmdViaSubprocess(lsabgc_ready_cmd, logObject, check_directories=[lsabgc_ready_results_directory])

	if not os.path.isfile(annotation_listing_file):
		try:
			annotation_listing_file = lsabgc_ready_results_directory + 'Primary_Sample_Annotation_Files.txt'
			assert(os.path.isfile(annotation_listing_file))
		except:
			raise RuntimeError("Issues validating the presence of the sample annotation files listing expected!")

	if not os.path.isdir(gcf_listing_dir):
		try:
			gcf_listing_dir = lsabgc_ready_results_directory + 'GCF_Listings/'
			assert(os.path.isdir(gcf_listing_dir))
		except:
			raise RuntimeError("Issues validating the presence of a directory with GCF listings!")

	if not os.path.isfile(orthogroups_matrix_file):
		try:
			orthogroups_matrix_file = lsabgc_ready_results_directory + 'Orthogroups.tsv'
		except:
			raise RuntimeError("Issue validating the presence of an OrthoFinder homolog group by sample matrix file!")

	try:
		assert(os.path.isfile(species_tree_file) and os.path.isfile(expected_similarities_file) and os.path.isfile(samples_to_keep_file))
	except:
		raise RuntimeError("Issue validating GToTree construction in lsaBGC-Ready was successful.")

	# Step 7: Perform dereplication and population clustering if needed
	population_file = outdir + 'Genomes_to_Populations_Mapping.txt'
	if not os.path.isfile(population_file):
		all_samples = set([])
		population_pairs = []
		samples_in_population_pairs = set([])
		dereplication_pairs = []
		samples_in_dereplication_pairs = set([])
		with open(expected_similarities_file) as oedf:
			for line in oedf:
				sample1, sample2, exsim = line.strip().split('\t')
				exsim = float(exsim)
				all_samples.add(sample1)
				all_samples.add(sample2)
				if exsim >= population_threshold:
					population_pairs.append([sample1, sample2])
					samples_in_population_pairs.add(sample1)
					samples_in_population_pairs.add(sample2)
				if exsim >= dereplicate_threshold:
					dereplication_pairs.append([sample1, sample2])
					samples_in_dereplication_pairs.add(sample1)
					samples_in_dereplication_pairs.add(sample2)

		if not skip_dereplication_flag:
			n50_cmds = []
			n50_directory = outdir + 'N50_Calculations/'
			util.setupReadyDirectory([n50_directory])
			for s, fna in all_genome_listings_fna.items():
				cmd = ['n50', fna, '>', n50_directory + s + '.txt', logObject]
				n50_cmds.append(cmd)

			p = multiprocessing.Pool(cores)
			p.map(util.multiProcess, n50_cmds)
			p.close()

			sample_n50s = defaultdict(lambda: 0.0)
			for f in os.listdir(n50_directory):
				sample = '.txt'.join(f.split('.txt')[:-1])
				with open(n50_directory + f) as onf:
					for line in onf:
						n50 = float(line.strip())
						sample_n50s[sample] = n50

			"""	
			Solution for single-linkage clustering taken from mimomu's response in the stackoverflow page:
			https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements?lq=1
			"""
			L = list(dereplication_pairs)
			LL = set(itertools.chain.from_iterable(L))
			for each in LL:
				components = [x for x in L if each in x]
				for i in components:
					L.remove(i)
				L += [list(set(itertools.chain.from_iterable(components)))]

			for s in all_samples:
				if not s in samples_in_dereplication_pairs:
					L.append([s])

			samples_to_keep_file = outdir + 'Dereplicated_Set_of_Genomes.txt'
			samples_to_keep_handle = open(samples_to_keep_file, 'w')
			for sc in L:
				cluster_sample_n50 = []
				for s in sc:
					cluster_sample_n50.append([s, n50])
				for i, s in enumerate(sorted(cluster_sample_n50, key=itemgetter(1), reverse=True)):
					if i == 0:
						samples_to_keep_handle.write(s[0] + '\n')
						break
			samples_to_keep_handle.close()

		"""	
		Solution for single-linkage clustering taken from mimomu's response in the stackoverflow page:
		https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements?lq=1
		"""
		L = list(population_pairs)
		LL = set(itertools.chain.from_iterable(L))
		for each in LL:
			components = [x for x in L if each in x]
			for i in components:
				L.remove(i)
			L += [list(set(itertools.chain.from_iterable(components)))]

		for s in all_samples:
			if not s in samples_in_population_pairs:
				L.append([s])

		pop_spec_handle = open(population_file, 'w')
		for i, sc in enumerate(L):
			pop_id = 'Clade_' + str(i)
			for s in sc:
				pop_spec_handle.write(s + '\t' + pop_id + '\n')
		pop_spec_handle.close()


	# Step 8: Run lsaBGC-AutoAnalyze.py
	lsabgc_autoanalyze_dir = outdir + 'lsaBGC_AutoAnalyze_Results/'
	lsabgc_autoanalyze_results_dir = lsabgc_autoanalyze_dir + 'Final_Results/'
	if not os.path.isdir(lsabgc_autoanalyze_results_dir):
		os.system('rm -rf %s' % lsabgc_autoanalyze_dir)
		lsabgc_autoanalyze_cmd = ['lsaBGC-AutoAnalyze.py', '-i', annotation_listing_file, '-g', gcf_listing_dir, '-m',
								  orthogroups_matrix_file, '-s', species_tree_file, '-w', expected_similarities_file, '-k',
								  samples_to_keep_file, '-c', str(cores), '-o', lsabgc_autoanalyze_dir, '-p',
								  bgc_prediction_software, '-u', population_file]
		runCmdViaSubprocess(lsabgc_autoanalyze_cmd, logObject, check_directories=[lsabgc_autoanalyze_results_dir])

	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

def runCmdViaSubprocess(cmd, logObject, check_files=[], check_directories=[]):
	logObject.info('Running %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		for cf in check_files:
			assert(os.path.isfile(cf))
		for cd in check_directories:
			assert(os.path.isdir(cd))
		logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))

def checkUlimitSettings(num_genomes, logObject):
	# Check ulimit settings!
	try:
		num_files = num_genomes * num_genomes
		uH = subprocess.check_output('ulimit -Hn', shell=True)
		uS = subprocess.check_output('ulimit -Sn', shell=True)
		uH = uH.decode('utf-8').strip()
		uS = uS.decode('utf-8').strip()
		if uH != 'unlimited':
			uH = int(uH)
		else:
			uH = 1e15
		if uS != 'unlimited':
			uS = int(uS)
		else:
			uS = 1e15
		if num_files > uS:
			logObject.error("Too many files will be produced and need to be read at once by OrthoFinder2. Luckily, you can resolve this quite easily via: ulimit -n 250000 . After running the command rerun lsaBGC-Easy.py and you should not get stuck here again.")
			sys.stderr.write("Too many files will be produced and need to be read at once by OrthoFinder2. Luckily, you can resolve this quite easily via: ulimit -n 250000 . After running the command rerun lsaBGC-Easy.py and you should not get stuck here again.")
			exit(1)
		if num_files > uH:
			logObject.error("Too many files will be produced and need to be read at once by OrthoFinder2. Your system requires root privs to change this, which I do not recommend. See the following OrthoFinder2 Github issue for more details: https://github.com/davidemms/OrthoFinder/issues/384")
			sys.stderr.write("Too many files will be produced and need to be read at once by OrthoFinder2. Your system requires root privs to change this, which I do not recommend. See the following OrthoFinder2 Github issue for more details: https://github.com/davidemms/OrthoFinder/issues/384")
			exit(1)

		if num_files < uS and num_files < uH:
			logObject.info(
				"Great news! Your ulimit settings look good and we believe OrthoFinder2 should be able to run smoothly!")
	except:
		logObject.error(
			"Difficulties validating ulimit settings are properly set to allow for successful OrthoFinder2 run.")
		raise RuntimeError(
			"Difficulties validating ulimit settings are properly set to allow for successful OrthoFinder2 run.")

def checkModelIsValid(gtotree_model):
	try:
		assert (gtotree_model in set(['Actinobacteria', 'Alphaproteobacteria', 'Bacteria', 'Archaea',
									  'Bacteria_and_Archaea', 'Bacteroidetes', 'Betaproteobacteria', 'Chlamydiae',
									  'Cyanobacteria', 'Epsilonproteobacteria', 'Firmicutes', 'Gammaproteobacteria',
									  'Proteobacteria', 'Tenericutes', 'Universal_Hug_et_al']))
	except:
		raise RuntimeError('Model for GToTree specified is not a valid model!')

if __name__ == '__main__':
    lsaBGC_Easy()