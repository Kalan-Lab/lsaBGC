#!/usr/bin/env python

### Program: lsaBGC-HMMExpander.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Microbiology and Immunology

import os
import sys
from time import sleep
import argparse
from lsaBGC import lsaBGC


def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-HMMExpander.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Microbiology and Immunology

	This program will automatically run or create task files for running prokka (gene calling and annotation), 
	antiSMASH (biosynthetic gene cluster annotation), and OrthoFinder (de novo ortholog group construction).
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--assemblies', type=str,
						help="Tab delimited text file. First column is the sample name and the second is the path to its assembly in FASTA format. Please remove troublesome characters in the sample name.",
						required=True)
	parser.add_argument('-o', '--outdir', type=str, help="The resulting output directory.", required=True)
	parser.add_argument('-ae', '--antiSMASH_env_path', type=str,
						help="Path to conda enviroment for antiSMASH. Database should automatically configured for antiSMASH loaded by the environment.",
						required=True)
	parser.add_argument('-pe', '--prokka_env_path', type=str, help="Path to conda enviroment for Prokka.",
						required=True)
	parser.add_argument('-oe', '--orthofinder_env_path', type=str, help="Path to conda enviroment for OrthoFinder.",
						required=True)
	parser.add_argument('-cp', '--conda_path', type=str,
						help="Path to anaconda/miniconda installation directory itself.", required=True)
	parser.add_argument('-l', '--lineage', type=str, help="The lineage under investigation.", required=False,
						default="Lineage")
	parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=8)
	parser.add_argument('-d', '--dry_run', action='store_true',
						help="Just create task files with commands for running prodigal, antiSMASH, and OrthoFinder. Useful for parallelizing across an HPC.",
						required=False, default=False)
	args = parser.parse_args()
	return args


def createInputsForLsaBGC():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	assembly_listing_file = os.path.abspath(myargs.assemblies)
	outdir = os.path.abspath(myargs.outdir) + '/'
	antiSMASH_env_path = os.path.abspath(myargs.antiSMASH_env_path) + '/'
	prokka_env_path = os.path.abspath(myargs.prokka_env_path) + '/'
	orthofinder_env_path = os.path.abspath(myargs.orthofinder_env_path) + '/'
	conda_path = myargs.conda_path

	antismash_load_code = '. %s/etc/profile.d/conda.sh && %s activate %s &&' % (conda_path, antiSMASH_env_path)
	prokka_load_code = '. %s/etc/profile.d/conda.sh && %s activate %s &&' % (conda_path, prokka_env_path)
	orthofinder_load_code = '. %s/etc/profile.d/conda.sh && %s activate %s &&' % (conda_path, orthofinder_env_path)

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	try:
		assert (os.path.isfile(assembly_listing_file))
	except:
		raise RuntimeError(
			'Input file listing the location of assemblies for samples is corrupt or an incorrect path was provided. Exiting now ...')

	"""
	PARSE OPTIONAL INPUTS
	"""

	cores = myargs.cores
	lineage = myargs.lineage
	dry_run_flag = myargs.dry_run

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = lsaBGC.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for future provedance.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [assembly_listing_file, outdir, cores, dry_run_flag, antiSMASH_env_path, prokka_env_path,
						orthofinder_env_path]
	parameter_names = ["Assembly Listing File", "Output Directory", "Cores", "Dry Run Flagged", "AntiSMASH Env Path",
					   "Prokka Env Path", "OrthoFinder Env Path"]
	lsaBGC.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Parse sample assemblies
	logObject.info("Parse sample assemblies from listing file.")
	sample_assemblies = lsaBGC.readInAssemblyListing(assembly_listing_file, logObject)
	logObject.info("Successfully parsed sample assemblies.")

	# Step 2: Run Prokka for gene calling
	prokka_outdir = outdir + 'Prokka_Results/'
	prokka_proteomes_dir = prokka_outdir + 'Prokka_Proteomes/'
	prokka_genbanks_dir = prokka_outdir + 'Prokka_Genbanks/'
	try:
		os.system('mkdir %s %s %s' % (prokka_outdir, prokka_proteomes_dir, prokka_genbanks_dir))
	except:
		logObject.error("Can't create Prokka results directories. Exiting now ...")
		raise RuntimeError("Can't create Prokka results directories. Exiting now ...")

	logObject.info("Running/setting-up Prokka for all samples!")
	lsaBGC.runProkka(sample_assemblies, prokka_outdir, prokka_proteomes_dir, prokka_genbanks_dir, prokka_load_code,
					 dry_run_flag, lineage, cores, logObject)
	logObject.info("Successfully ran/set-up Prokka.")

	# Step 3: Run AntiSMASH for identification of BGCs
	antismash_outdir = outdir + 'AntiSMASH_Results/'
	try:
		os.system('mkdir %s' % (antismash_outdir))
	except:
		logObject.error("Can't create AntiSMASH results directories. Exiting now ...")
		raise RuntimeError("Can't create AntiSMASH results directories. Exiting now ...")

	logObject.info("Running/setting-up AntiSMASH for all samples!")
	lsaBGC.runAntiSMASH(prokka_genbanks_dir, antismash_outdir, antismash_load_code, dry_run_flag, cores, logObject)
	logObject.info("Successfully ran/set-up AntiSMASH.")

	# Step 4: Run OrthoFinder for de novo ortholog construction
	orthofinder_outdir = outdir + 'OrthoFinder_Results/'
	logObject.info("Running/setting-up OrthoFinder!")
	lsaBGC.runOrthoFinder(prokka_proteomes_dir, orthofinder_outdir, orthofinder_load_code, dry_run_flag, cores,
						  logObject)
	logObject.info("Successfully ran/set-up OrthoFinder.")

	# Close logging object, write resulting list of BGC Genbanks, and exit!
	lsaBGC.closeLoggerObject(logObject)
	sys.exit(0)


if __name__ == '__main__':
	createInputsForLsaBGC()