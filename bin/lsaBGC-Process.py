#!/usr/bin/env python

### Program: lsaBGC-Process.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Microbiology and Immunology

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

	parser.add_argument('-i', '--assembly_listing', type=str,
						help="Tab delimited text file. First column is the sample name and the second is the path to its assembly in FASTA format. Please remove troublesome characters in the sample name.",
						required=True)
	parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
	parser.add_argument('-ae', '--antiSMASH_env_path', type=str,
						help="Path to conda environment for antiSMASH. Database should automatically configured for antiSMASH loaded by the environment.",
						required=True)
	parser.add_argument('-pe', '--prokka_env_path', type=str, help="Path to conda environment for Prokka.",
						required=True)
	parser.add_argument('-cp', '--conda_path', type=str,
						help="Path to anaconda/miniconda installation directory itself.", required=True)
	parser.add_argument('-oe', '--orthofinder_env_path', type=str, help="Path to conda environment for OrthoFinder. Optional, if not used, locus tags will be 3 characters insteado just 2.",
						required=False, default=None)
	parser.add_argument('-l', '--lineage', type=str, help="The lineage under investigation.", required=False,
						default="Lineage")
	parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=8)
	parser.add_argument('-d', '--dry_run', action='store_true',
						help="Just create task files with commands for running prodigal, antiSMASH, and OrthoFinder. Useful for parallelizing across an HPC.",
						required=False, default=False)
	parser.add_argument('-s', '--skip_annotation', action='store_true', help="Skip basic/standard annotation in Prokka", required=False, default=False)
	parser.add_argument('-e', '--extract_protein_fasta', action='store_true', help="Whether to extract protein fasta and append as third column path to such FASTA in file listing BGCs Genbanks.", required=False, default=False)
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
	antiSMASH_env_path = os.path.abspath(myargs.antiSMASH_env_path) + '/'
	prokka_env_path = os.path.abspath(myargs.prokka_env_path) + '/'
	conda_path = myargs.conda_path

	antismash_load_code = '. %s/etc/profile.d/conda.sh && conda activate %s &&' % (conda_path, antiSMASH_env_path)
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
	lineage = myargs.lineage
	dry_run_flag = myargs.dry_run
	skip_annotation_flag = myargs.skip_annotation
	orthofinder_env_path = None
	orthofinder_load_code = None
	extract_protein_fasta = myargs.extract_protein_fasta

	locus_tag_length = 4
	if myargs.orthofinder_env_path:
		locus_tag_length = 3
		orthofinder_env_path = os.path.abspath(myargs.orthofinder_env_path) + '/'
		orthofinder_load_code = '. %s/etc/profile.d/conda.sh && conda activate %s &&' % (conda_path, orthofinder_env_path)

	"""
	START WORKFLOW
	"""
	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for future provedance.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [assembly_listing_file, outdir, cores, dry_run_flag, skip_annotation_flag, antiSMASH_env_path, prokka_env_path,
						orthofinder_env_path, extract_protein_fasta]
	parameter_names = ["Assembly Listing File", "Output Directory", "Cores", "Dry Run Flagged", "Skip Prokka Annotation",
					   "AntiSMASH Env Path", "Prokka Env Path", "OrthoFinder Env Path", "Extract BGC Proteins as FASTA"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Parse sample assemblies
	logObject.info("Parse sample assemblies from listing file.")
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
					 lineage, cores, locus_tag_length, logObject, dry_run_flag=dry_run_flag, skip_annotation_flag=skip_annotation_flag)
	logObject.info("Successfully ran/set-up Prokka.")

	# Step 3: Run AntiSMASH for identification of BGCs
	antismash_outdir = outdir + 'AntiSMASH_Results/'
	bgc_proteomes_outdir = outdir + 'BGC_Proteomes/'
	try:
		if extract_protein_fasta:
			os.system('mkdir %s' % (bgc_proteomes_outdir))
		os.system('mkdir %s' % (antismash_outdir))
	except:
		logObject.error("Can't create AntiSMASH results directories. Exiting now ...")
		raise RuntimeError("Can't create AntiSMASH results directories. Exiting now ...")

	logObject.info("Running/setting-up AntiSMASH for all samples!")
	processing.runAntiSMASH(prokka_genbanks_dir, antismash_outdir, antismash_load_code, cores, logObject, dry_run_flag=dry_run_flag)
	logObject.info("Successfully ran/set-up AntiSMASH.")

	# Step 4: Run OrthoFinder for de novo ortholog construction
	if orthofinder_env_path:
		orthofinder_outdir = outdir + 'OrthoFinder_Results/'
		logObject.info("Running/setting-up OrthoFinder!")
		processing.runOrthoFinder(prokka_proteomes_dir, orthofinder_outdir, orthofinder_load_code, cores, logObject, dry_run_flag=dry_run_flag)
		logObject.info("Successfully ran/set-up OrthoFinder.")

	# Write resulting list of BGC Genbanks (to be used as input for lsaBGC-Cluster)
	antismash_bgc_listing_file = outdir + 'All_AntiSMASH_BGCs.txt'
	antismash_bgc_listing_handle = open(antismash_bgc_listing_file, 'w')
	for s in os.listdir(antismash_outdir):
		sample_antismash_outdir = antismash_outdir + s + '/'
		if not os.path.isdir(sample_antismash_outdir): continue
		for f in os.listdir(sample_antismash_outdir):
			if '.region' in f and f.endswith('.gbk'):
				bgc_genbank = sample_antismash_outdir + f
				if extract_protein_fasta:
					bgc_proteome = processing.extractBGCProteomes(s, bgc_genbank, bgc_proteomes_outdir, logObject)
					antismash_bgc_listing_handle.write(s + '\t' + bgc_genbank + '\t' + bgc_proteome + '\n')
				else:
					antismash_bgc_listing_handle.write(s + '\t' + bgc_genbank + '\n')
			else:
				os.system('rm -f %s' % sample_antismash_outdir + f)
	antismash_bgc_listing_handle.close()

	if orthofinder_env_path:
		# Move select result files from OrthoFinder to main directory to make more easy to access/find
		orthofinder_homolog_matrix = orthofinder_outdir + 'Orthogroups.csv'
		orthofinder_species_tree = [orthofinder_outdir + od for od in os.listdir(orthofinder_outdir) if od.startswith("Orthologues_")][0] + '/SpeciesTree_rooted.txt'
		if os.path.isfile(orthofinder_species_tree):
			os.system('mv %s %s' % (orthofinder_species_tree, outdir))
		if os.path.isfile(orthofinder_homolog_matrix):
			os.system('mv %s %s' % (orthofinder_homolog_matrix, outdir))

	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	lsaBGC_Process()