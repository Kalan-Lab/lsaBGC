#!/usr/bin/env python

### Program: lsaBGC-Process.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Microbiology and Immunology

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
from lsaBGC import processing, util

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Process.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Microbiology and Immunology
	
	This program will automatically run or create task files for running Prokka (gene calling and annotation), 
	antiSMASH (biosynthetic gene cluster annotation), and OrthoFinder (de novo ortholog group construction).
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-a', '--assembly_listing', type=str,
						help="Tab delimited text file. First column is the sample name and the second is the path to its assembly in FASTA format. Please remove troublesome characters in the sample name.",
						required=True)
	parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
	parser.add_argument('-cp', '--conda_path', type=str, help="Path to anaconda/miniconda installation directory itself.", required=True)
	parser.add_argument('-pe', '--prokka_env_path', type=str, help="Path to conda environment for Prokka.", required=True)
	parser.add_argument('-oe', '--orthofinder_env_path', type=str, help="Path to conda environment for OrthoFinder. Optional, if not used, locus tags will be 3 characters insteado just 2.",
						required=False, default=None)
	parser.add_argument('-ae', '--antiSMASH_env_path', type=str,
						help="Path to conda environment for antiSMASH. Database should automatically configured for antiSMASH loaded by the environment.",
						required=False)
	parser.add_argument('-g', '--genus', type=str, help="The genus under investigation. The lineage of interest could be species, but for this, just use the genus.", required=False,
						default="Genus")
	parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=8)
	parser.add_argument('-d', '--dry_run', action='store_true',
						help="Just create task files with commands for running prodigal, antiSMASH, and OrthoFinder. Useful for parallelizing across an HPC.",
						required=False, default=False)
	parser.add_argument('-s', '--append_singleton_hgs', action='store_true', help="Append homolog groups with only one protein representative to the Orthogroups.csv homolog group matrix. This enables more reliable detection of homologous rare/singleton BGCs downstream in the pipeline.", required=False, default=False)
	parser.add_argument('-q', '--fast_annotation', action='store_true', help="Skip basic/standard annotation in Prokka.", required=False, default=False)
	parser.add_argument('-p', '--only_run_prokka', action='store_true', help="Only run Prokka for gene annotation and Genbank creation. Skip the rest.", required=False, default=False)
	parser.add_argument('-f', '--refined_orthofinder', action='store_true', help="Only run OrthoFinder on proteins from antiSMASH proteomes only. This has implications downstream on being able to identify multi-copy genes across the genome.", required=False, default=False)
	args = parser.parse_args()
	return args

def lsaBGC_Process():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	assembly_listing_file = os.path.abspath(myargs.assembly_listing)
	outdir = os.path.abspath(myargs.output_directory) + '/'
	prokka_env_path = os.path.abspath(myargs.prokka_env_path) + '/'
	conda_path = myargs.conda_path

	prokka_load_code = '. %s/etc/profile.d/conda.sh && conda activate %s &&' % (conda_path, prokka_env_path)

	try:
		assert (os.path.isfile(assembly_listing_file))
	except:
		raise RuntimeError(
			'Input file listing the location of assemblies for samples is corrupt or an incorrect path was provided. Exiting now ...')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	"""
	PARSE OPTIONAL INPUTS
	"""

	cores = myargs.cores
	lineage = myargs.genus
	dry_run_flag = myargs.dry_run
	fast_annotation_flag = myargs.fast_annotation
	only_run_prokka = myargs.only_run_prokka
	refined_orthofinder = myargs.refined_orthofinder
	append_singleton_hgs_flag = myargs.append_singleton_hgs

	orthofinder_env_path = None
	orthofinder_load_code = None
	antiSMASH_env_path = None
	antiSMASH_load_code = None
	locus_tag_length = 4

	if not only_run_prokka:
		locus_tag_length = 3
		try:
			orthofinder_env_path = os.path.abspath(myargs.orthofinder_env_path) + '/'
			orthofinder_load_code = '. %s/etc/profile.d/conda.sh && conda activate %s &&' % (conda_path, orthofinder_env_path)

			antiSMASH_env_path = os.path.abspath(myargs.antiSMASH_env_path) + '/'
			antiSMASH_load_code = '. %s/etc/profile.d/conda.sh && conda activate %s &&' % (conda_path, antiSMASH_env_path)
		except Exception as e:
			raise RuntimeError("Had issues validating the conda environments for antiSMASH and OrthoFinder are valid. Exiting ...")

	"""
	START WORKFLOW
	"""
	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for future provenance.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [assembly_listing_file, outdir, cores, dry_run_flag, fast_annotation_flag, prokka_env_path,
						antiSMASH_env_path, orthofinder_env_path, only_run_prokka, refined_orthofinder,
						append_singleton_hgs_flag]
	parameter_names = ["Assembly Listing File", "Output Directory", "Cores", "Dry Run Flagged",
					   "Fast Prokka Annotation Requested?", "Prokka Env Path", "AntiSMASH Env Path",
					   "OrthoFinder Env Path", "Only Prokka Annotations to be Run?",
					   "Run OrthoFinder with Only BGC Proteins?",
					   "Append Singleton HGs to the Homolog Group by Sample Matrix?"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Parse sample assemblies
	logObject.info("Beginning to parse sample assemblies from listing file.")
	sample_assemblies = processing.readInAssemblyListing(assembly_listing_file, logObject)
	logObject.info("Successfully parsed sample assemblies.")

	# Step 2: Run Prokka for gene calling
	prokka_outdir = outdir + 'Prokka_Results/'
	prokka_proteomes_dir = prokka_outdir + 'Prokka_Proteomes/'
	prokka_genbanks_dir = prokka_outdir + 'Prokka_Genbanks/'
	try:
		if not os.path.isdir(prokka_outdir): os.system('mkdir %s' % prokka_outdir)
		if not os.path.isdir(prokka_proteomes_dir): os.system('mkdir %s' % prokka_proteomes_dir)
		if not os.path.isdir(prokka_genbanks_dir): os.system('mkdir %s' % prokka_genbanks_dir)
	except:
		logObject.error("Can't create Prokka results directories. Exiting now ...")
		raise RuntimeError("Can't create Prokka results directories. Exiting now ...")

	logObject.info("Running/setting-up Prokka for all samples!")
	processing.runProkka(sample_assemblies, prokka_outdir, prokka_proteomes_dir, prokka_genbanks_dir, prokka_load_code,
					 lineage, cores, locus_tag_length, logObject, dry_run_flag=dry_run_flag, skip_annotation_flag=fast_annotation_flag)

	prokka_results_listing_file = outdir + 'Sample_Annotation_Files.txt'
	prlf_handle = open(prokka_results_listing_file, 'w')
	for s in sample_assemblies:
		try:
			prokka_genbank = prokka_genbanks_dir + s + '.gbk'
			prokka_proteome = prokka_proteomes_dir + s + '.faa'
			assert(os.path.isfile(prokka_genbank) + os.path.isfile(prokka_proteome))
			prlf_handle.write(s + '\t' + prokka_genbank + '\t' + prokka_proteome + '\n')
		except Exception as e:
			logObject.warning("Unable to validate Prokka ran successfully for sample %s, skipping ..." % s)
	prlf_handle.close()
	logObject.info("Successfully ran/set-up Prokka.")

	if not only_run_prokka:
		# Step 3: Run AntiSMASH for identification of BGCs
		antismash_outdir = outdir + 'AntiSMASH_Results/'
		try:
			os.system('mkdir %s' % (antismash_outdir))
		except Exception as e:
			logObject.error("Can't create AntiSMASH results directories. Exiting now ...")
			raise RuntimeError("Can't create AntiSMASH results directories. Exiting now ...")

		logObject.info("Running/setting-up AntiSMASH for all samples!")
		processing.runAntiSMASH(prokka_genbanks_dir, antismash_outdir, antiSMASH_load_code, cores, logObject,
								dry_run_flag=dry_run_flag)
		logObject.info("Successfully ran/set-up AntiSMASH.")

		refined_proteomes_outdir = outdir + 'AntiSMASH_Sample_Proteomes/'
		if refined_orthofinder:
			os.system('mkdir %s' % refined_proteomes_outdir)

		# Write resulting list of BGC Genbanks (to be used as input for lsaBGC-Cluster)
		antismash_bgc_listing_file = outdir + 'All_AntiSMASH_BGCs.txt'
		antismash_bgc_listing_handle = open(antismash_bgc_listing_file, 'w')
		for s in os.listdir(antismash_outdir):
			sample_antismash_outdir = antismash_outdir + s + '/'
			if not os.path.isdir(sample_antismash_outdir): continue
			sample_bgcs = set([])
			for f in os.listdir(sample_antismash_outdir):
				if '.region' in f and f.endswith('.gbk'):
					bgc_genbank = sample_antismash_outdir + f
					sample_bgcs.add(bgc_genbank)
					antismash_bgc_listing_handle.write(s + '\t' + bgc_genbank + '\n')
				else:
					os.system('rm -fr %s' % sample_antismash_outdir + f)
			if refined_orthofinder:
				util.writeRefinedProteomes(s, sample_bgcs, refined_proteomes_outdir, logObject)
		antismash_bgc_listing_handle.close()

		# Step 4: Run OrthoFinder for de novo ortholog construction
		orthofinder_outdir = outdir + 'OrthoFinder_Results/'
		logObject.info("Running/setting-up OrthoFinder!")
		if not refined_orthofinder:
			processing.runOrthoFinder(prokka_proteomes_dir, orthofinder_outdir, orthofinder_load_code, cores, logObject, dry_run_flag=dry_run_flag)
		else:
			processing.runOrthoFinder(refined_proteomes_outdir, orthofinder_outdir, orthofinder_load_code, cores, logObject, dry_run_flag=dry_run_flag)
		logObject.info("Successfully ran/set-up OrthoFinder.")

		logObject.info("Organizing results directory.")

		# Move select result files from OrthoFinder to main directory to make more easy to access/find
		orthofinder_homolog_matrix = orthofinder_outdir + 'Orthogroups/Orthogroups.tsv'
		flag = False
		if os.path.isfile(orthofinder_homolog_matrix) and append_singleton_hgs_flag:
			unassigned_orthofinder_homolog_matrix = orthofinder_outdir + 'Orthogroups/Orthogroups_UnassignedGenes.tsv'
			result_file = outdir + 'Orthogroups.tsv'
			processing.appendSingletonHGsToPresenceMatrix(orthofinder_homolog_matrix, unassigned_orthofinder_homolog_matrix, result_file, logObject)
			flag = True
		elif os.path.isfile(orthofinder_homolog_matrix):
			os.system('mv %s %s' % (orthofinder_homolog_matrix, outdir))
			flag = True

		if not flag:
			raise RuntimeWarning("Orthofinder did not run successfully, please check logs in output directory.")
		orthofinder_species_tree = orthofinder_outdir + 'Species_Tree/SpeciesTree_rooted.txt'
		if os.path.isfile(orthofinder_species_tree):
			os.system('mv %s %s' % (orthofinder_species_tree, outdir))

		logObject.info('Main Resulting Files:')
		logObject.info('Homolog Group Presence Absence/Matrix OrthoFinder2: ' + orthofinder_homolog_matrix)
		logObject.info('Listing of BGCs predicted by AntiSMASH: ' + antismash_bgc_listing_file)
		logObject.info('Listing of Sample Proteome and Genbanks from Prokka: ' + prokka_results_listing_file)

	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	lsaBGC_Process()