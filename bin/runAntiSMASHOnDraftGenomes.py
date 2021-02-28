#!/usr/bin/env python

### Program: runAntiSMASHOnDraftGenomes.py
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
	Program: runAntiSMASHOnDraftGenomes.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Microbiology and Immunology
	
	This program will automatically run barebones AntiSMASH analysis to annotate BGCs on genomic assemblies.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--assemblies', type=str,
						help="Tab delimited text file. First column is the sample name and the second is the path to its assembly in FASTA format. Please remove troublesome characters in the sample name.",
						required=True)
	parser.add_argument('-o', '--outdir', type=str, help="The resulting output directory.", required=True)
	parser.add_argument('-ae', '--antiSMASH_env_path', type=str,
						help="Path to conda enviroment for antiSMASH. Database should automatically configured for antiSMASH loaded by the environment.",
						required=True)
	parser.add_argument('-cp', '--conda_path', type=str,
						help="Path to anaconda/miniconda installation directory itself.", required=True)
	parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)
	parser.add_argument('-d', '--dry_run', action='store_true',
						help="Just create task files with commands for running prodigal, antiSMASH, and OrthoFinder. Useful for parallelizing across an HPC.",
						required=False, default=False)
	args = parser.parse_args()
	return args

def runAntiSMASHOnDraftGenomes():
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
	conda_path = myargs.conda_path
	antismash_load_code = '. %s/etc/profile.d/conda.sh && conda activate %s &&' % (conda_path, antiSMASH_env_path)

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
	parameter_values = [assembly_listing_file, outdir, cores, dry_run_flag, antiSMASH_env_path]
	parameter_names = ["Assembly Listing File", "Output Directory", "Cores", "Dry Run Flagged", "AntiSMASH Env Path"]
	lsaBGC.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Parse sample assemblies
	logObject.info("Parse sample assemblies from listing file.")
	sample_assemblies = lsaBGC.readInAssemblyListing(assembly_listing_file, logObject)
	logObject.info("Successfully parsed sample assemblies.")

	# Step 2: Run AntiSMASH for identification of BGCs
	antismash_outdir = outdir + 'AntiSMASH_Results/'
	bgc_proteomes_outdir = outdir + 'BGC_Proteomes/'
	try:
		os.system('mkdir %s %s' % (antismash_outdir, bgc_proteomes_outdir))
	except:
		logObject.error("Can't create AntiSMASH results directories. Exiting now ...")
		raise RuntimeError("Can't create AntiSMASH results directories. Exiting now ...")

	logObject.info("Running/setting-up AntiSMASH for all samples!")
	lsaBGC.runAntiSMASHBareBones(sample_assemblies, antismash_outdir, antismash_load_code, dry_run_flag, cores, logObject)
	logObject.info("Successfully ran/set-up AntiSMASH.")

	# Step 3: Extract proteomes from BGC genbanks and write resulting list of BGC Genbanks
	antismash_bgc_gbk_listing_file = outdir + 'All_AntiSMASH_BGCs.txt'
	antismash_bgc_pro_listing_file = outdir + 'All_AntiSMASH_BGC_Proteomes.txt'
	abg_handle = open(antismash_bgc_gbk_listing_file, 'w')
	abp_handle = open(antismash_bgc_pro_listing_file, 'w')
	for s in os.listdir(antismash_outdir):
		sample_antismash_outdir = antismash_outdir + s + '/'
		if not os.path.isdir(sample_antismash_outdir): continue
		for f in os.listdir(sample_antismash_outdir):
			if '.region' in f and f.endswith('.gbk'):
				bgc_genbank = sample_antismash_outdir + f
				bgc_proteome = lsaBGC.extractBGCProteomes(s, bgc_genbank, bgc_proteomes_outdir, logObject)
				abg_handle.write(s + '\t' + bgc_genbank + '\n')
				abp_handle.write(s + '\t' + bgc_proteome + '\n')
			else:
				os.system('rm -rf %s' % sample_antismash_outdir + f)
	abg_handle.close()
	abp_handle.close()

	# Close logging object and exit
	lsaBGC.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	runAntiSMASHOnDraftGenomes()