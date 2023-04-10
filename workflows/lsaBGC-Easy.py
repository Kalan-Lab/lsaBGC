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
import subprocess
import traceback
import multiprocessing
import math
import itertools
from operator import itemgetter
from collections import defaultdict
import gzip
import resource
from time import sleep
from lsaBGC import util
from ete3 import Tree

os.environ['OMP_NUM_THREADS'] = '4'
lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	   __              ___   _____  _____
	  / /  ___ ___ _  / _ ) / ___/ / ___/
	 / /  (_-</ _ `/ / _  |/ (_ / / /__  
	/_/  /___/\_,_/ /____/ \___/  \___/  
	*******************************************************************************************************************
	Program: lsaBGC-Easy.py (Version 1.31)
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	
	QUICK DESCRIPTION: 
	Workflow to run the majority of lsaBGC functionalities for a specific taxa or manually defined set of user-provided
	genomes or both. For information on the workflow and examples on how to run please take a look at the Wiki:
	https://github.com/Kalan-Lab/lsaBGC/wiki/14.-lsaBGC-Easy-Tutorial:-Combining-lsaBGC-with-ncbi-genome-download
		
	Currently only works for bacteria, but lsaBGC can handle fungi now - just not lsaBGC-Easy ... yet
	*******************************************************************************************************************
	CONSIDERATIONS:	
	* Check out the considerations wiki page: https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background-&-Considerations  	
	* If interested in running just a directory of genomes (you don't want to downoad genomes from Genbank for some taxa)
	  then you can just set the required -n argument to "None".
	*******************************************************************************************************************
	OVERVIEW OF STEPS TAKEN:	
		- Check number of genomes for taxa is not too crazy (<1000) or too small (>10)
		- Step 1: Get set of Genbank assembly accessions from recent GTDB release matching taxa.
		- Step 2: Download all genomes in FASTA format using ncbi-genome-download and perform gene-calling with prodigal.
		- Step 3: Run GToTree, Dereplicate, Group Samples into Populations/Clades, and Create genomes listing
		- Step 4: Run GECCO based annotation of BGCs and crete BGC listing or Create Task File with antiSMASH/DeepBGC coomands
		- Step 5: Run lsaBGC-Ready.py with lsaBGC-Cluster or BiG-SCAPE
		- Step 6: Run lsaBGC-AutoExpansion.py polishing to find GCF instances fragmented on multiple scaffolds
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

	parser.add_argument('-n', '--taxa_name',
						help='Name of the taxa of interest as listed in GTDB. Currently restricted to bacteria - but lsaBGC\nitself- when using antiSMASH BGC predictions can run on\nEuks - check out the tutorial!! If there is a space in the\nname, please surround by quotes.',
						required=True)
	parser.add_argument('-o', '--output_directory', help='Parent output/workspace directory.', required=True)
	parser.add_argument('-g', '--user_genomes_directory',
						help='A directory with additional genomes, e.g. those recently sequenced by the\nuser, belonging to the taxa. Accepted formats include FASTA.\nAccepted suffices include: .fna, .fa, .fasta.',
						required=False, default=None)
	parser.add_argument('-p', '--bgc_prediction_software',
						help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is GECCO - which will be automatic. For the other two,\nlsaBGC-Easy will produce a task-file which you will need to run in a seperate environment\nwith antiSMASH or DeepBGC installed, then rerun lsaBGC-Easy.py using the original command..',
						default='gecco', required=False)
	parser.add_argument('-om', '--orthofinder_mode',
						help="Method for running OrthoFinder2. (Options: Genome_Wide, BGC_Only).\nDefault is Genome_Wide (BGC_Only is slightly experimental but much faster and should work\nespecially well for taxa with many BGCs).",
						default='Genome_Wide', required=False)
	parser.add_argument('-sd', '--skip_dereplication', action='store_true',
						help="Whether to skip dereplication based on GToTree alignments of SCGs - not\nrecommended and can cause issues if there are a lot of\ngenomes for the taxa of interest.",
						default=False, required=False)
	parser.add_argument('-gtm', '--gtotree_model',
						help="SCG model for secondary GToTree analysis and wdhat would be used for dereplication. Default is to use the \"Bacteria\" SCG set (n=74 genes).",
						default='Bacteria', required=False)
	parser.add_argument('-iib', '--include_incomplete_bgcs', action='store_true',
						help="Whether to account for incomplete BGCs prior to clustering - not recommended.",
						default=False, required=False)
	parser.add_argument('-mc', '--run_coarse_orthofinder', action='store_true',
						help='Use coarse clustering of homolog groups in OrthoFinder instead of more\nresolute hierarchical determined homolog groups. There are some advantages to coarse OGs, including their construction being deterministic.',
						required=False, default=False)
	parser.add_argument('-b', '--use_bigscape', action='store_true',
						help="Use BiG-SCAPE for BGC clustering into GCFs instead of lsaBGC-Cluster. Recommended if you want to include incomplete BGCs for clustering and are using antiSMASH.",
						required=False, default=False)
	parser.add_argument('-sae', '--skip_auto_expansion', action='store_true',
						help="Skip lsaBGC-AutoExpansion.py to find missing pieces of BGCs due to assembly fragmentation.",
						required=False, default=False)
	parser.add_argument('-c', '--cpus', type=int,
						help="Total number of cpus/threads. Note, this is the total number of\nthreads to use (Default is 4). BGC prediction commands (via GECCO, antiSMASH, and DeepBGC) will\neach be set to use 4 cpus by default.",
						required=False, default=4)
	parser.add_argument('-dt', '--dereplicate_threshold', type=float,
						help="Amino acid similarity threshold of SCGs for considering\ntwo genomes as redundant. (Default is 0.999)",
						default=0.999, required=False)
	parser.add_argument('-pt', '--population_threshold', type=float,
						help="Amino acid similarity threshold of SCGs for considering\ntwo genomes as belonging to the same population. (Default is 0.99)",
						default=0.99, required=False)
	parser.add_argument('-x', '--ignore_limits', action='store_true',
						help="Ignore limitations on number of genomes allowed.\nE.g. allow for analyses of taxa with more than 2000 genomes available and more than 100 genomes\nafter dereplication. Not recommend, be cautious!!! Also note,\nyou can always delete \"Dereplicated_Set_of_Genomes.txt\" in the results directory and redo\ndereplication with different threshold.")
	parser.add_argument('-py', '--use_pyrodigal', action='store_true', help='Use pyrodigal instead of prodigal.',
						required=False, default=False)

	args = parser.parse_args()
	return args


def lsaBGC_Easy():
	myargs = create_parser()

	taxa_name = myargs.taxa_name.strip('"').strip()
	outdir = os.path.abspath(myargs.output_directory) + '/'
	cpus = myargs.cpus
	user_genomes_directory = myargs.user_genomes_directory
	bgc_prediction_software = myargs.bgc_prediction_software.upper()
	gtotree_model = myargs.gtotree_model
	skip_dereplication_flag = myargs.skip_dereplication
	include_incomplete_bgcs_flag = myargs.include_incomplete_bgcs
	skip_auto_expansion_flag = myargs.skip_auto_expansion
	dereplicate_threshold = myargs.dereplicate_threshold
	population_threshold = myargs.population_threshold
	run_coarse_orthofinder = myargs.run_coarse_orthofinder
	orthofinder_mode = myargs.orthofinder_mode.upper()
	use_bigscape_flag = myargs.use_bigscape
	ignore_limits_flag = myargs.ignore_limits
	use_pyrodigal = myargs.use_pyrodigal

	try:
		assert (orthofinder_mode in set(['GENOME_WIDE', 'BGC_ONLY']))
	except:
		sys.stderr.write('BGC prediction software option is not a valid option.\n')
		sys.exit(1)
	try:
		assert (bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
	except:
		sys.stderr.write('BGC prediction software option is not a valid option.\n')
		sys.exit(1)
	if user_genomes_directory:
		try:
			user_genomes_directory = os.path.abspath(user_genomes_directory) + '/'
			assert (os.path.isdir(user_genomes_directory))
		except:
			sys.stderr.write('User defined genomes directory could not be validated to exist!\n')
			sys.exit(1)

	checkModelIsValid(gtotree_model)

	if os.path.isdir(outdir):
		sys.stderr.write(
			"Output directory already exists! Overwriting in 5 seconds, but only where needed will use checkpoint files to skip certain steps...\n")
		sleep(5)
	else:
		util.setupReadyDirectory([outdir])

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)
	parameters_file = outdir + 'Command_Issued.txt'
	logo = """
   __              ___   _____  _____
  / /  ___ ___ _  / _ ) / ___/ / ___/
 / /  (_-</ _ `/ / _  |/ (_ / / /__  
/_/  /___/\_,_/ /____/ \___/  \___/  
*************************************"""
	version_string = util.parseVersionFromSetupPy()
	sys.stdout.write(logo + '\n\n')
	sys.stdout.write('Running version %s\n' % version_string)
	sys.stdout.write("Appending command issued for future records to: %s\n" % parameters_file)
	sys.stdout.write("Logging more details at: %s\n" % log_file)
	logObject.info("\nNEW RUN!!!\n**************************************")
	logObject.info('Running version %s' % version_string)
	logObject.info("Appending command issued for future records to: %s" % parameters_file)

	parameters_handle = open(parameters_file, 'a+')
	parameters_handle.write(' '.join(sys.argv) + '\n')
	parameters_handle.close()

	parallel_jobs_4cpu = max(math.floor(cpus / 4), 1)

	# Step 1: parse GTDB information file, get list of Genbank accessions, and perform dry-run with ncbi-genome-download
	genbank_accession_listing_file = outdir + 'NCBI_Genbank_Accession_Listing.txt'
	sys.stdout.write(
		"--------------------\nStep 1\n--------------------\nBeginning by assessing which genomic assemblies are available for the taxa %s in GTDB and NCBI's Genbank db\n" % taxa_name)
	logObject.info(
		"\n--------------------\nStep 1\n--------------------\nBeginning by assessing which genomic assemblies are available for the taxa %s in GTDB and NCBI's Genbank db" % taxa_name)
	if not os.path.isfile(genbank_accession_listing_file):
		genbank_accession_listing_handle = open(genbank_accession_listing_file, 'w')
		with gzip.open(lsaBGC_main_directory + 'db/GTDB_R207_Information.txt.gz', 'rt') as ogtdb:
			for line in ogtdb:
				line = line.strip('\n')
				ls = line.split('\t')
				if ls[0] == 'none': continue
				if len(taxa_name.split()) == 1:
					if ls[1] == taxa_name:
						genbank_accession_listing_handle.write(ls[0] + '\n')
				elif len(taxa_name.split()) == 2:
					if ls[2] == taxa_name:
						genbank_accession_listing_handle.write(ls[0] + '\n')
		genbank_accession_listing_handle.close()

	ogalf = open(genbank_accession_listing_file)
	accession_count = len(ogalf.readlines())
	ogalf.close()

	if (accession_count >= 2000 or accession_count < 10) and (not ignore_limits_flag):
		logObject.error(
			"Currently not recommended to use lsaBGC-Easy for taxa with > 2000 genomes or < 10 genomes. In the future this will likely be updated to be a warning, but would rather not risk harming your server!")
		sys.stderr.write(
			"Currently not recommended to use lsaBGC-Easy for taxa with > 2000 genomes or < 10 genomes.\nIn the future this will likely be updated to be a warning, but would rather not risk harming your computer/server!\n")
		sys.stderr.write(
			'Exiting now, but you can move past this if you feel comfortable and are willing to wait for longer processing times using the -x flag.\n')
		sys.exit(1)

	genome_listing_file = outdir + 'NCBI_Genomes_from_Genbank_for_Taxa.txt'
	if not os.path.isfile(genome_listing_file):
		if accession_count != 0:
			ngd_dry_cmd = ['ncbi-genome-download', '--dry-run', '--section', 'genbank', '-A',
						   genbank_accession_listing_file,
						   'bacteria', '>', genome_listing_file]
			runCmdViaSubprocess(ngd_dry_cmd, logObject, check_files=[genome_listing_file])
		else:
			of = open(genome_listing_file, 'w')
			of.close()

	oglf = open(genome_listing_file)
	genome_count = len(oglf.readlines())
	oglf.close()

	if (genome_count >= 2000 or genome_count < 10) and (not ignore_limits_flag):
		logObject.error(
			"Currently not recommended to use lsaBGC-Easy for taxa with > 2000 genomes or < 10 genomes. In the future this will likely be updated to be a warning, but would rather not risk harming your server!")
		sys.stderr.write(
			"Currently not recommended to use lsaBGC-Easy for taxa with > 2000 genomes or < 10 genomes.\nIn the future this will likely be updated to be a warning, but would rather not risk harming your computer/server!\n")
		sys.stderr.write(
			'Exiting now, but you can move past this if you feel comfortable and are willing to wait for longer processing times using the -x flag.\n')
		sys.exit(1)

	# Step 2: Download all genomes in FASTA format & prodigal gene calling
	genomes_directory = outdir + 'genbank/'
	all_genomes_listing_file = outdir + 'All_Genomes_Listing.txt'
	uncompress_dir = outdir + 'Uncompressed_Genomes/'
	if not os.path.isfile(all_genomes_listing_file):
		if genome_count != 0:
			ngd_real_cmd = ['ncbi-genome-download', '--formats', 'fasta', '--retries', '2', '--section',
							'genbank', '-A', genbank_accession_listing_file, 'bacteria', '-o', outdir]
			runCmdViaSubprocess(ngd_real_cmd, logObject, check_directories=[genomes_directory])
			list_all_genomes_cmd = ['listAllGenomesInDirectory.py', '-i', genomes_directory, '-u', uncompress_dir, '-z',
									'-d', genome_listing_file, '>', all_genomes_listing_file]
			runCmdViaSubprocess(list_all_genomes_cmd, logObject, check_files=[all_genomes_listing_file])

		if user_genomes_directory:
			list_all_user_genomes_cmd = ['listAllGenomesInDirectory.py', '-i', user_genomes_directory,
										 '-u', uncompress_dir, '-z', '>>', all_genomes_listing_file]
			runCmdViaSubprocess(list_all_user_genomes_cmd, logObject)

	if not os.path.isfile(all_genomes_listing_file):
		logObject.error('No genomes downloaded / provided. Exiting as more genomes are needed.')
		sys.stderr.write('No genomes downloaded / provided. Exiting ...\n')
		sys.exit(1)
	else:
		oaglf = open(all_genomes_listing_file)
		genome_included_count = len(oaglf.readlines())
		oaglf.close()
		if genome_included_count <= 2:
			logObject.error('Fewer than 3 genomes downloaded / provided. Exiting as more genomes are needed.')
			sys.stderr.write('Fewer than 3 genomes downloaded / provided. Exiting as more genomes are needed.\n')
			sys.exit(1)
		else:
			logObject.info('Great, we found or downloaded %d genomes belonging to the taxa!' % (genome_included_count))
			sys.stdout.write(
				'Great, we found or downloaded %d genomes belonging to the taxa!\n' % (genome_included_count))

	prodigal_results_directory = outdir + 'Prodigal_and_Processing_Results/'
	all_proteomes_listing_file = outdir + 'All_Proteomes_Listing.txt'
	logObject.info('\n--------------------\nStep 2\n--------------------\nBeginning gene calling and processing.')
	sys.stdout.write('--------------------\nStep 2\n--------------------\nBeginning gene calling and processing.\n')
	if not os.path.isfile(all_proteomes_listing_file):
		util.setupReadyDirectory(([prodigal_results_directory]))
		run_prodigal_cmds = []
		with open(all_genomes_listing_file) as oaglf:
			for line in oaglf:
				line = line.strip()
				sample, fna = line.split('\t')
				try:
					assert (os.path.isfile(fna))
				except:
					raise RuntimeError("Could not validate the genome for sample %s exists." % sample)
				rpc = ['runProdigalAndMakeProperGenbank.py', '-i', fna, '-s', sample, '-o', prodigal_results_directory,
					   logObject]
				run_prodigal_cmds.append(rpc)

		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, run_prodigal_cmds)
		p.close()

		all_proteomes_listing_handle = open(all_proteomes_listing_file, 'w')
		for f in os.listdir(prodigal_results_directory):
			if not f.endswith('.faa'): continue
			all_proteomes_listing_handle.write(prodigal_results_directory + f + '\n')
		all_proteomes_listing_handle.close()

		logObject.info('Gene calling and processing completed.')
		sys.stdout.write('Gene calling and processing completed.\n')
	else:
		logObject.info(
			'Listing of predicted proteomes (from prior gene-calling) already found in directory, will skip reprocessing.')
		sys.stdout.write(
			'Listing of predicted proteomes (from prior gene-calling) already found in directory, will skip reprocessing.\n')

	all_genome_listings_gbk = {}
	for f in os.listdir(prodigal_results_directory):
		if not f.endswith('.gbk'): continue
		sample = f.split('.gbk')[0]
		all_genome_listings_gbk[sample] = prodigal_results_directory + f

	all_genome_listings_fna = {}
	with open(all_genomes_listing_file) as oglf:
		for line in oglf:
			line = line.strip()
			sample, fna_path = line.split('\t')
			all_genome_listings_fna[sample] = fna_path

	# Step 3: Construct GToTree phylogeny - dereplication + grouping - creating listing
	gtotree_outdir = outdir + 'GToTree_output/'
	species_tree_file = gtotree_outdir + 'GToTree_output.tre'
	expected_similarities_file = outdir + 'GToTree_Expected_Similarities.txt'
	logObject.info(
		'\n--------------------\nStep 3\n--------------------\nBeginning construction of species tree using GToTree.')
	sys.stdout.write(
		'--------------------\nStep 3\n--------------------\nBeginning construction of species tree using GToTree.\n')
	if not os.path.isfile(species_tree_file) or not os.path.isfile(expected_similarities_file):
		os.system('rm -rf %s' % gtotree_outdir)
		gtotree_cmd = ['GToTree', '-A', all_proteomes_listing_file, '-H', gtotree_model, '-n', '4', '-j',
					   str(parallel_jobs_4cpu), '-M', '4', '-o', gtotree_outdir]
		runCmdViaSubprocess(gtotree_cmd, logObject, check_files=[species_tree_file])

		protein_msa_file = gtotree_outdir + 'Aligned_SCGs.faa'
		pair_seq_matching = util.determineSeqSimProteinAlignment(protein_msa_file)
		expected_sim_handle = open(expected_similarities_file, 'w')
		for s1 in pair_seq_matching:
			for s2 in pair_seq_matching:
				exp_sim = str(pair_seq_matching[s1][s2])
				if s1 == s2: exp_sim = '1.0'
				expected_sim_handle.write(s1 + '\t' + s2 + '\t' + exp_sim + '\n')
		expected_sim_handle.close()
		logObject.info('construction of species tree using GToTree completed.')
		sys.stdout.write('Construction of species tree using GToTree completed.\n')
	else:
		logObject.info('A prior species tree was found so skipping re-running GToTree.')
		sys.stdout.write('A prior species tree was found so skipping re-running GToTree.\n')

	population_file = outdir + 'Genomes_to_Populations_Mapping.txt'
	samples_to_keep_file = outdir + 'Dereplicated_Set_of_Genomes.txt'
	if not os.path.isfile(population_file) or not os.path.isfile(samples_to_keep_file):
		logObject.info(
			'Beginning determination of populations/clades and dereplication (assuming dereplication not requested to be skipped via -sd flag).')
		sys.stdout.write(
			'Beginning determination of populations/clades and dereplication (assuming dereplication not requested to be skipped via -sd flag).\n')

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

		n50_cmds = []
		n50_directory = outdir + 'N50_Calculations/'
		util.setupReadyDirectory([n50_directory])
		for s, fna in all_genome_listings_fna.items():
			cmd = ['n50', fna, '>', n50_directory + s + '.txt', logObject]
			n50_cmds.append(cmd)

		p = multiprocessing.Pool(cpus)
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

		samples_to_keep_handle = open(samples_to_keep_file, 'w')
		for sc in L:
			cluster_sample_n50 = []
			for s in sc:
				cluster_sample_n50.append([s, n50])
			for i, s in enumerate(sorted(cluster_sample_n50, key=itemgetter(1, 0), reverse=True)):
				if i == 0:
					samples_to_keep_handle.write(s[0] + '\n')
				elif skip_dereplication_flag:
					samples_to_keep_handle.write(s[0] + '\n')
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

		assert (os.path.isfile(samples_to_keep_file))
		assert (os.path.isfile(population_file))
		logObject.info('Completed population/clade inference and dereplication.')
		logObject.info('Samples kept for analysis can be found at: %s.' % samples_to_keep_file)
		logObject.info('Population/clading info can be found at: %s.' % population_file)
		sys.stdout.write('Completed population/clade inference and dereplication.\n')
		sys.stdout.write('Samples kept for analysis can be found at: %s.\n' % samples_to_keep_file)
		sys.stdout.write('Population/clading info can be found at: %s.\n' % population_file)
	else:
		logObject.info(
			'Found prior determined population/clade inference and sample listing files. Using them instead of reassessing them.')
		sys.stdout.write(
			'Found prior determined population/clade inference and sample listing files. Using them instead of reassessing them.\n')

	count_of_dereplicated_sample_set = 0
	with open(samples_to_keep_file) as oskf:
		for line in oskf:
			count_of_dereplicated_sample_set += 1
	if (count_of_dereplicated_sample_set > 100) and (not ignore_limits_flag):
		logObject.error(
			"Currently not recommended to use lsaBGC-Easy for taxa with > 100 genomes after dereplication. In the future this will likely be updated to be a warning, but would rather not risk harming your computer/server!")
		sys.stderr.write(
			"Currently not recommended to use lsaBGC-Easy for taxa with > 100 genomes\nafter dereplication. In the future this will likely be updated\nto be a warning, but would rather not risk harming\nyour computer/server!")
		sys.stderr.write(
			'Exiting now, but you can move past this if you feel comfortable and are willing to wait for longer processing times using the -x flag.\n')
	else:
		logObject.info(
			"After dereplication (if performed) there are %d samples/genomes being considered." % count_of_dereplicated_sample_set)
		sys.stdout.write(
			"After dereplication (if performed) there are %d samples/genomes being considered.\n" % count_of_dereplicated_sample_set)

	checkUlimitSettings(count_of_dereplicated_sample_set, logObject)

	if bgc_prediction_software != 'GECCO':
		logObject.info("\n--------------------\nStep 4\n--------------------\nStarting creation of %s commands for de novo BGC prediction." % bgc_prediction_software)
		sys.stdout.write("--------------------\nStep 4\n--------------------\nStarting creation of %s commands for de novo BGC prediction.\n" % bgc_prediction_software)
	else:
		logObject.info("\n--------------------\nStep 4\n--------------------\nRunning GECCO for de novo BGC predictions.")
		sys.stdout.write("--------------------\nStep 4\n--------------------\nRunning GECCO for de novo BGC predictions.\n")

	primary_genomes_listing_file = outdir + 'Primary_Genomes.txt'
	primary_bgc_pred_directory = outdir + 'Primary_Genomes_BGC_Predictions/'
	cmd_file = outdir + 'bgc_prediction.cmds'
	if not (os.path.isfile(primary_genomes_listing_file)):
		primary_genomes_listing_handle = open(primary_genomes_listing_file, 'w')
		util.setupReadyDirectory([primary_bgc_pred_directory])
		primary_bgc_pred_cmds = []
		with open(samples_to_keep_file) as oskf:
			for s in oskf:
				s = s.strip()
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
		primary_genomes_listing_handle.close()

		# Step 4: Run GECCO or print commands for others and exit based annotation of BGCs and crete BGC listing
		if bgc_prediction_software == 'GECCO':
			p = multiprocessing.Pool(parallel_jobs_4cpu)
			p.map(util.multiProcess, primary_bgc_pred_cmds)
			p.close()
		else:
			cmd_handle = open(cmd_file, 'w')
			for cmd in primary_bgc_pred_cmds:
				cmd_handle.write(' '.join(cmd) + '\n')
			cmd_handle.close()
			logObject.info(
				'%s BGC prediction commands written to %s.Please run these and afterwards restart lsaBGC-Easy.py with the same command afterwards.' % (
				bgc_prediction_software, cmd_file))
			sys.stdout.write(
				'************************\nPlease run the following BGC prediction commands using a different conda envirnoment with the BGC prediction\nsoftware (%s) installed and afterwards restart lsaBGC-Easy.py with the same command as used initially.\nExiting now, see you back soon, here is the task file with the %s commands:\n%s\n' % (
				bgc_prediction_software, bgc_prediction_software, cmd_file))
			sys.exit(0)

	bgc_prediciton_count = 0
	for f in os.listdir(primary_bgc_pred_directory):
		bgc_prediciton_count += 1

	if bgc_prediciton_count == 0 and (bgc_prediction_software == 'ANTISMASH' or bgc_prediction_software == 'DEEPBGC'):
		logObject.error(
			"************************\nLooks like you didn't run the BGC prediction commands (or at least successfully),\nplease run them in a different environment where %s is present in the environment.\nOnly after they have been successfully run and written to the directory:\n%s\ncan you restart the lsaBGC command. Here are the %s commands again:\n%s" % (
			bgc_prediction_software, primary_bgc_pred_directory, bgc_prediction_software, cmd_file))
		sys.stderr.write(
			"************************\nLooks like you didn't run the BGC prediction commands (or at least successfully),\nplease run them in a different environment where %s is present in the environment.\nOnly after they have been successfully run and written to the directory:\n%s\ncan you restart the lsaBGC command. Here are the %s commands again:\n%s\n" % (
			bgc_prediction_software, primary_bgc_pred_directory, bgc_prediction_software, cmd_file))
		sys.exit(1)

	primary_bgcs_listing_file = outdir + 'BGCs_in_Primary_Genomes.txt'
	if not os.path.isfile(primary_bgcs_listing_file):
		list_primary_bgc_cmd = ['listAllBGCGenbanksInDirectory.py', '-i', primary_bgc_pred_directory, '-p',
								bgc_prediction_software, '>', primary_bgcs_listing_file]
		if not include_incomplete_bgcs_flag:
			list_primary_bgc_cmd += ['-f']
		runCmdViaSubprocess(list_primary_bgc_cmd, logObject, check_files=[primary_bgcs_listing_file])

	logObject.info("Primary BGC predictions appear successful and are listed for samples/genomes at:\n%s" % primary_bgcs_listing_file)
	sys.stdout.write("Primary BGC predictions appear successful and are listed for samples/genomes at:\n%s\n" % primary_bgcs_listing_file)

	# Step 5: Run lsaBGC-Ready.py with lsaBGC-Cluster or BiG-SCAPE
	logObject.info('\n--------------------\nStep 5\n--------------------\nBeginning clustering of BGCs using either lsaBGC-Cluster (default) or BiG-SCAPE (can be requested).')
	sys.stdout.write('--------------------\nStep 5\n--------------------\nBeginning clustering of BGCs using either lsaBGC-Cluster (default) or BiG-SCAPE (can be requested).\n')
	bigscape_listing_file = lsaBGC_main_directory + 'external_tools/bigscape_location.txt'
	bigscape_prog_location, pfam_directory = [None] * 2
	if os.path.isfile(bigscape_listing_file):
		bigscape_prog_location, pfam_directory = open(bigscape_listing_file).readlines()[0].strip().split('\t')
	bigscape_results_dir = outdir + 'BiG_SCAPE_Clustering_Results/'
	if use_bigscape_flag and not os.path.isdir(bigscape_results_dir) and os.path.isfile(bigscape_prog_location):
		bigscape_cmd = ['python', bigscape_prog_location, '-i', primary_bgc_pred_directory, '-o', bigscape_results_dir,
						'-c', str(cpus), '--include_singletons', '--pfam_dir', pfam_directory]
		runCmdViaSubprocess(bigscape_cmd, logObject, check_directories=[bigscape_results_dir])
		try:
			assert (os.path.isdir(bigscape_results_dir + 'network_files/'))
		except:
			logObject.error("Something unexpected in BiG-SCAPE processing occurred.")
			sys.stderr.write("Something unexpected in BiG-SCAPE processing occurred.\n")
			exit(1)

	lsabgc_ready_directory = outdir + 'lsaBGC_Ready_Results/'
	lsabgc_ready_results_directory = lsabgc_ready_directory + 'Final_Results/'
	annotation_listing_file = lsabgc_ready_results_directory + 'Primary_Sample_Annotation_Files.txt'
	gcf_listing_dir = lsabgc_ready_results_directory + 'GCF_Listings/'
	orthogroups_matrix_file = lsabgc_ready_results_directory + 'Orthogroups.tsv'
	if not os.path.isdir(lsabgc_ready_results_directory) or len(
			[f for f in os.listdir(lsabgc_ready_results_directory)]) <= 2:
		os.system('rm -rf %s' % lsabgc_ready_directory)
		lsabgc_ready_cmd = ['lsaBGC-Ready.py', '-i', primary_genomes_listing_file, '-l', primary_bgcs_listing_file,
							'-o', lsabgc_ready_directory, '-c', str(cpus), '-p', bgc_prediction_software,
							'-m', orthofinder_mode]

		if use_bigscape_flag and os.path.isdir(bigscape_results_dir):
			lsabgc_ready_cmd += ['-b', bigscape_results_dir]
		else:
			lsabgc_ready_cmd += ['-lc']
		if use_pyrodigal:
			lsabgc_ready_cmd += ['-py']
		if run_coarse_orthofinder:
			lsabgc_ready_cmd += ['-mc']

		db_listing_file = lsaBGC_main_directory + '/db/database_location_paths.txt'
		if os.path.isfile(db_listing_file):
			lsabgc_ready_cmd += ['-a']
		runCmdViaSubprocess(lsabgc_ready_cmd, logObject,
							check_directories=[lsabgc_ready_results_directory, gcf_listing_dir],
							check_files=[annotation_listing_file, orthogroups_matrix_file])

	# Step 6: Run lsaBGC-AutoExpansion.py polishing to find GCF instances fragmented on multiple scaffolds
	if not skip_auto_expansion_flag:
		logObject.info('\n--------------------\nStep 6\n--------------------\nBeginning identification of fragmented BGC instances from completed instances using lsaBGC-AutoExpansion.')
		sys.stdout.write('--------------------\nStep 6\n--------------------\nBeginning identification of fragmented BGC instances from completed instances using lsaBGC-AutoExpansion.\n')

	checkLsaBGCInputsExist(annotation_listing_file, gcf_listing_dir, orthogroups_matrix_file)

	lsabgc_autoexpansion_results_directory = outdir + 'lsaBGC_AutoExpansion_Results/'
	exp_annotation_listing_file = lsabgc_autoexpansion_results_directory + 'Sample_Annotation_Files.txt'
	exp_orthogroups_matrix_file = lsabgc_autoexpansion_results_directory + 'Orthogroups.expanded.tsv'
	exp_gcf_listing_dir = lsabgc_autoexpansion_results_directory + 'Updated_GCF_Listings/'

	if not os.path.isdir(exp_gcf_listing_dir) and (not skip_auto_expansion_flag):
		lsabgc_expansion_cmd = ['lsaBGC-AutoExpansion.py', '-g', gcf_listing_dir, '-m', orthogroups_matrix_file, '-l',
								annotation_listing_file, '-e', annotation_listing_file, '-q', '-c', str(cpus),
								'-o', lsabgc_autoexpansion_results_directory, '-p', bgc_prediction_software]
		runCmdViaSubprocess(lsabgc_expansion_cmd, logObject,
							check_directories=[lsabgc_autoexpansion_results_directory, exp_gcf_listing_dir],
							check_files=[exp_orthogroups_matrix_file, exp_orthogroups_matrix_file])

	if skip_auto_expansion_flag:
		exp_annotation_listing_file = annotation_listing_file
		exp_orthogroups_matrix_file = orthogroups_matrix_file
		exp_gcf_listing_dir = gcf_listing_dir

	# Step 7: Run lsaBGC-AutoAnalyze.py
	logObject.info('\n--------------------\nStep 7\n--------------------\nBeginning lsaBGC-AutoExpansion (which runs lsaBGC-See.py, lsaBGC-PopGene.py, and lsaBGC-Divergence.py + summarizes results across GCFs at the end.\n')
	sys.stdout.write('--------------------\nStep 7\n--------------------\nBeginning lsaBGC-AutoExpansion (which runs lsaBGC-See.py, lsaBGC-PopGene.py,\nand lsaBGC-Divergence.py + summarizes results across GCFs at the end.\n')

	checkLsaBGCInputsExist(exp_annotation_listing_file, exp_gcf_listing_dir, exp_orthogroups_matrix_file)
	try:
		assert (os.path.isfile(species_tree_file) and os.path.isfile(expected_similarities_file) and os.path.isfile(
			samples_to_keep_file))
	except:
		raise RuntimeError("Issue validating GToTree construction in lsaBGC-Ready was successful.")

	lsabgc_autoanalyze_dir = outdir + 'lsaBGC_AutoAnalyze_Results/'
	lsabgc_autoanalyze_results_dir = lsabgc_autoanalyze_dir + 'Final_Results/'
	if not os.path.isdir(lsabgc_autoanalyze_results_dir):
		os.system('rm -rf %s' % lsabgc_autoanalyze_dir)
		lsabgc_autoanalyze_cmd = ['lsaBGC-AutoAnalyze.py', '-i', exp_annotation_listing_file, '-g', exp_gcf_listing_dir,
								  '-m',
								  exp_orthogroups_matrix_file, '-s', species_tree_file, '-w',
								  expected_similarities_file, '-k',
								  samples_to_keep_file, '-c', str(cpus), '-o', lsabgc_autoanalyze_dir, '-p',
								  bgc_prediction_software, '-u', population_file]
		runCmdViaSubprocess(lsabgc_autoanalyze_cmd, logObject, check_directories=[lsabgc_autoanalyze_results_dir])


	gseef_results_dir = outdir + 'GSeeF_Results/'
	gseef_final_results_dir = gseef_results_dir + 'Final_Results/'
	if not os.path.isdir(gseef_results_dir):
		os.system('rm -rf %s' % gseef_results_dir)

		pruned_species_tree_file = outdir + 'Species_Tree.pruned.tre'
		try:
			samples_to_keep_list = []
			with open(samples_to_keep_file) as ostkf:
				for line in ostkf:
					line = line.strip()
					samples_to_keep_list.append(line.strip())
			t = Tree(species_tree_file)
			t.prune(samples_to_keep_list)
			t.write(format=1, outfile=pruned_species_tree_file)
		except Exception as e:
			sys.stderr.write('Issues with creating pruned species tree for GSeeF analysis!\n')
			sys.stderr.write(e + '\n')
			logObject.error('Issues with creating pruned species tree for GSeeF analysis!')
			sys.exit(1)

		gseef_cmd = ['GSeeF.py', '-l', exp_annotation_listing_file, '-g', exp_gcf_listing_dir,
					 '-s', pruned_species_tree_file,  '-c', str(cpus),
					 '-o', gseef_results_dir, '-p', bgc_prediction_software]
		runCmdViaSubprocess(gseef_cmd, logObject, check_directories=[gseef_final_results_dir])

	# Close logging object and exit
	logObject.info('lsaBGC-Easy completed! Check out the major results in the folder: %s' % lsabgc_autoanalyze_results_dir)
	sys.stdout.write('lsaBGC-Easy completed! Check out the major results (from lsaBGC-AutoAnalyze) in the folder:\n%s\n' % lsabgc_autoanalyze_results_dir)
	sys.stdout.write('And for a coarse view of GCFs across a species phylogeny check out the GSeeF output at:\n%s\n' % gseef_final_results_dir)
	util.closeLoggerObject(logObject)
	sys.exit(0)


def runCmdViaSubprocess(cmd, logObject, check_files=[], check_directories=[]):
	logObject.info('Running %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		for cf in check_files:
			assert (os.path.isfile(cf))
		for cd in check_directories:
			assert (os.path.isdir(cd))
		logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))


def checkUlimitSettings(num_genomes, logObject):
	# Check ulimit settings!
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
		logObject.warning(
			"Because many files will be produced and need to be read at once by OrthoFinder2, we are increasing the current shell's limits! Five second pause to allow you to exit the program if you do not which to continue")
		sleep(5)
		resource.setrlimit(resource.RLIMIT_NOFILE, (250000, 250000))
	if num_files > uH:
		logObject.error(
			"Too many files will be produced and need to be read at once by OrthoFinder2. Your system requires root privs to change this, which I do not recommend. See the following OrthoFinder2 Github issue for more details: https://github.com/davidemms/OrthoFinder/issues/384")
		sys.stderr.write(
			"Too many files will be produced and need to be read at once by OrthoFinder2. Your system requires root privs to change this, which I do not recommend. See the following OrthoFinder2 Github issue for more details: https://github.com/davidemms/OrthoFinder/issues/384")
		sys.exit()
	if num_files < uS and num_files < uH:
		logObject.info(
			"Great news! Your ulimit settings look good and we believe OrthoFinder2 should be able to run smoothly!")


def checkModelIsValid(gtotree_model):
	try:
		assert (gtotree_model in set(['Actinobacteria', 'Alphaproteobacteria', 'Bacteria', 'Archaea',
									  'Bacteria_and_Archaea', 'Bacteroidetes', 'Betaproteobacteria', 'Chlamydiae',
									  'Cyanobacteria', 'Epsilonproteobacteria', 'Firmicutes', 'Gammaproteobacteria',
									  'Proteobacteria', 'Tenericutes', 'Universal_Hug_et_al']))
	except:
		raise RuntimeError('Model for GToTree specified is not a valid model!')


def checkLsaBGCInputsExist(annotation_listing_file, gcf_listing_dir, orthogroups_matrix_file):
	try:
		assert (os.path.isfile(annotation_listing_file))
	except:
		raise RuntimeError("Issues validating the presence of the sample annotation files listing expected!")

	try:
		assert (os.path.isfile(orthogroups_matrix_file))
	except:
		raise RuntimeError("Issue validating the presence of an OrthoFinder homolog group by sample matrix file!")

	try:
		assert (os.path.isdir(gcf_listing_dir))
	except:
		raise RuntimeError("Issues validating the presence of a directory with GCF listings!")


if __name__ == '__main__':
	lsaBGC_Easy()
