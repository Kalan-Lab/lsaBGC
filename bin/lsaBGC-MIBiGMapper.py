#!/usr/bin/env python

### Program: lsaBGC-MIBiGMapper.py
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
from lsaBGC import util
from lsaBGC.classes.GCF import GCF
import subprocess
import traceback
from collections import defaultdict
from operator import itemgetter
from scipy import stats

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-MIBiGMapper.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will assess whether MIBiG BGCs can be grouped in with a single GCF as determined by lsaBGC-Cluster
	or BiG-SCAPE. 
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name\nwhile the 2nd column is the path to a BGC prediction\nin Genbank format.', required=True)
	parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder homolog by sample matrix.", required=True)
	parser.add_argument('-o', '--output_directory', help="Output directory.", required=True)
	parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO)\n[Default is antiSMASH].', default='antiSMASH', required=False)
	parser.add_argument('-id', '--identity', type=float, help="Minimal identity of MIBiG proteins to match GCF proteins to be\nassigned same homolog group [Default is 60.0].", required=False, default=60.0)
	parser.add_argument('-cov', '--coverage', type=float, help="Minimal coverage of GCF proteins needed to assign MIBiG protein\nto same homolog group [Default is 70.0].", required=False, default=70.0)
	parser.add_argument('-sh', '--shared', type=int, help="Minimal number of homolog groups of representative BGC found in MIBiG\nBGC to consider potential match [Default is 5].", required=False, default=5)
	parser.add_argument('-pr', '--syntenic_correlation_cutoff', type=float, help="Minimum absolute correlation coefficient of\nMIBiG BGC with representative BGC from GCF. Is not regarded if '--draft_mode' specified. [Default is 0.0].", required=False, default=0.0)
	parser.add_argument('-d', '--draft_mode', action='store_true', help='BGCs in GCF listing are fragmented. Will assess MIBiG\nBGCs to all instances of GCF in sample - assuming multiple BGCs are due to assembly\nfragmentation.', default=False, required=False)
	parser.add_argument('-c', '--cpus', type=int, help="Number of CPUs to use [Default is 1].", required=False, default=1)
	args = parser.parse_args()
	return args

def mapToMIBiG():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	gcf_listings_file = os.path.abspath(myargs.gcf_listing)
	gcf_id = myargs.gcf_id
	orthofinder_matrix_file = os.path.abspath(myargs.orthofinder_matrix)
	outdir = os.path.abspath(myargs.output_directory) + '/'

	### vet input files quickly
	try:
		assert(os.path.isfile(orthofinder_matrix_file))
		assert(os.path.isfile(gcf_listings_file))
	except:
		raise RuntimeError('One or more of the input files provided, does not exist. Exiting now ...')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	mibig_info_file = None
	mibig_db_file = None
	try:
		db_locations = lsaBGC_main_directory + 'db/database_location_paths.txt'
		assert (os.path.isfile(db_locations))
		with open(db_locations) as odls:
				for line in odls:
					line = line.strip()
					ls = line.split('\t')
					if ls[0] == 'mibig':
						mibig_info_file, mibig_db_file = ls[1:]
						assert (os.path.isfile(mibig_info_file) and os.path.isfile(mibig_db_file))
	except:
		raise RuntimeError("It appears the MIBiG database was not setup or setup successfully. Please run/rerun the script setup_annotation_dbs.py and report to Github issues if issue persists.")

	"""
	PARSE OPTIONAL INPUTS
	"""

	bgc_prediction_software = myargs.bgc_prediction_software.upper()
	identity_cutoff = float(myargs.identity)
	coverage_cutoff = float(myargs.coverage)
	shared_cutoff = float(myargs.shared)
	correlation_cutoff = float(myargs.syntenic_correlation_cutoff)
	draft_mode = myargs.draft_mode
	cpus = myargs.cpus

	try:
		assert (bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
	except:
		raise RuntimeError('BGC prediction software option is not a valid option.')

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)
	version_string = util.parseVersionFromSetupPy()
	logObject.info('Running lsaBGC version %s' % version_string)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for future records.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [gcf_listings_file, gcf_id, orthofinder_matrix_file, outdir, bgc_prediction_software,
						identity_cutoff, coverage_cutoff, shared_cutoff, correlation_cutoff, draft_mode, cpus]
	parameter_names = ["GCF Listing File", "GCF ID", "OrthoFinder Orthogroups.tsv File", "Output Directory",
					   "BGC Prediction Method", "MIBiG Protein to BGC Protein Identity Cutoff",
					   "MIBiG Protein to BGC Protein Coverage Cutoff", "MIBiG BGC to GCF BGC shared HGs Cutoff",
					   "MIBiG BGC to GCF BGC Syntenic Correlation Cutoff", "Run in Draft Mode?", "CPUs"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Extract Proteins from BGCs belonging to GCF.
	# Create GCF object
	GCF_Object = GCF(gcf_listings_file, gcf_id=gcf_id, logObject=logObject)

	# Process GenBanks
	logObject.info("Processing BGC GenBanks from GCF listing file.")
	GCF_Object.readInBGCGenbanks(comprehensive_parsing=True, prediction_method=bgc_prediction_software)
	logObject.info("Successfully parsed BGC GenBanks and associated with unique IDs.")

	# extract proteins from BGC GenBanks to FASTA file
	gcf_prots_file = outdir + 'GCF_proteins.faa'
	GCF_Object.aggregateProteins(gcf_prots_file, draft_mode=draft_mode)

	# Step 2: Perform DIAMOND alignment of protein.
	diamond_result_file = outdir + 'DIAMOND_Results.txt'
	diamond_cmd = ['diamond', 'blastp', '--threads', str(cpus), '--db', mibig_db_file, '--very-sensitive', '--query',
				   gcf_prots_file, '-k', '0', '--out', diamond_result_file, '--outfmt', '6', 'qseqid', 'sseqid', 'pident',
				   'evalue', 'bitscore', 'qcovhsp']
	try:
		subprocess.call(' '.join(diamond_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		assert (os.path.isfile(diamond_result_file))
		logObject.info('Successfully ran: %s' % ' '.join(diamond_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(diamond_cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(diamond_cmd))

	# Step 3: Parse OrthoFinder Homolog vs Sample Matrix and associate each homolog group with a color
	logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
	gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy = util.parseOrthoFinderMatrix(orthofinder_matrix_file, GCF_Object.pan_genes)
	GCF_Object.inputHomologyInformation(gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy)
	logObject.info("Successfully parsed homolog matrix.")

	# Step 4: Read in MIBiG Database - connect IDs to  (un-needed step - but maybe might be needed in the future
	# depending on formatting changes).
	mibig_id_to_desc = {}
	with open(mibig_info_file) as omif:
		for line in omif:
			line = line.strip()
			ls = line.split('\t')
			mibig_id_to_desc[ls[0]] = ls[1]

	# Step 5: Parse DIAMOND Results
	best_mibig_prot_hits = defaultdict(lambda: [[], [], 0.0])
	with open(diamond_result_file) as odrf:
		for line in odrf:
			line = line.strip()
			qid, sid = line.split('\t')[:2]
			pid, eval, bitscore, qcovhsp = [float(x) for x in line.split('\t')[2:]]
			if pid >= identity_cutoff and qcovhsp >= coverage_cutoff:
				mibig_prot = mibig_id_to_desc[sid]
				qid_hg = gene_to_hg[qid.split('|')[1]]
				if bitscore > best_mibig_prot_hits[mibig_prot][2]:
					best_mibig_prot_hits[mibig_prot] = [[qid], [qid_hg], bitscore]
				elif bitscore == best_mibig_prot_hits[mibig_prot][2]:
					best_mibig_prot_hits[mibig_prot][0].append(qid)
					best_mibig_prot_hits[mibig_prot][1].append(qid_hg)


	mibig_bgc_prots = defaultdict(set)
	mibig_prot_to_hg = {}
	for mbp in best_mibig_prot_hits:
		if len(set(best_mibig_prot_hits[mbp][1])) == 1:
			mibig_prot_to_hg[mbp] = best_mibig_prot_hits[mbp][1][0]
			mb = mbp.split('|')[0]
			mibig_bgc_prots[mb].add(mbp)
		else:
			logObject.warning('Note, MIBiG protein %s matched multiple homolog-groups, avoiding inclusion in downstream analyses.' % (mbp))

	sample_hgs = defaultdict(set)
	bgc_hgs = defaultdict(set)
	bgc_hg_ranked_orders = defaultdict(dict)
	for bgc in GCF_Object.bgc_hgs:
		sc_hg_coords = []
		for hg in GCF_Object.bgc_hgs[bgc]:
			sample = GCF_Object.bgc_sample[bgc]
			sample_hgs[sample].add(hg)
			bgc_hgs[bgc].add(hg)

			gene_counts = 0
			hg_midpoint = None
			hg_strand = None
			for g in GCF_Object.hg_genes[hg]:
				if GCF_Object.comp_gene_info[g]['bgc_name'] == bgc:
					hg_midpoint = (GCF_Object.comp_gene_info[g]['start'] + GCF_Object.comp_gene_info[g]['end'])/2.0
					hg_strand = GCF_Object.comp_gene_info[g]['direction']
					gene_counts += 1
			if gene_counts == 1:
				sc_hg_coords.append([hg, hg_midpoint, hg_strand])

		for hg_ord, hg_its in enumerate(sorted(sc_hg_coords, key=itemgetter(1))):
			bgc_hg_ranked_orders[bgc][hg_its[0]] = [hg_ord, hg_its[1], hg_its[2]]

	accounted_for_tups = set([])
	result_file = outdir + 'GCF_to_MIBiG_Relations.txt'
	rf_handle = open(result_file, 'w')
	rf_handle.write('gcf_id\tmibig_bgc_id\tgcf_hg\tmibig_protein\n')
	for mb in mibig_bgc_prots:
		mbp_hg_counts = defaultdict(int)
		for mbp in mibig_bgc_prots[mb]:
			mbp_hg = mibig_prot_to_hg[mbp]
			mbp_hg_counts[mbp_hg] += 1
		mhg_info = {}
		for mbp in mibig_bgc_prots[mb]:
			mbp_hg = mibig_prot_to_hg[mbp]
			if mbp_hg_counts[mbp_hg] != 1: continue
			mhg_midpoint = sum([float(x) for x in mbp.split('|')[2].split('-')])/2.0
			mhg_strand = mbp.split('|')[3]
			mhg_info[mbp_hg] = [mhg_midpoint, mhg_strand, mbp]
		if draft_mode:
			for samp in sample_hgs:
				gcf_samp_shared_prop = (len(sample_hgs[samp].intersection(set(mhg_info.keys())))/float(len(sample_hgs[samp])))*100.0
				gcf_samp_shared = len(sample_hgs[samp].intersection(set(mhg_info.keys())))
				if gcf_samp_shared >= shared_cutoff:
					mibig_prot_hg_tup = tuple([mhg_info[hg][2], hg])
					if not mibig_prot_hg_tup in accounted_for_tups:
						rf_handle.write('\t'.join([gcf_id, mb, hg, mhg_info[hg][2]]) + '\n')
					accounted_for_tups.add(mibig_prot_hg_tup)
		else:
			for bgc in bgc_hgs:
				sample = GCF_Object.bgc_sample[bgc]
				bgc_shared_prop = ((len(bgc_hgs[bgc].intersection(set(mhg_info.keys()))))/float(len(bgc_hgs[bgc])))*100.0
				bgc_shared = len(bgc_hgs[bgc].intersection(set(mhg_info.keys())))
				if bgc_shared >= shared_cutoff:
					best_corr = None
					mib_hg_order = []
					mib_hg_direction = []
					ref_hg_order = []
					ref_hg_direction = []

					for hg in mhg_info:
						mhg_midpoint, mhg_strand, mbp = mhg_info[hg]
						mib_hg_order.append(mhg_midpoint)
						mib_hg_direction.append(mhg_strand)
						lt_matching = []
						for lt in GCF_Object.hg_genes[hg]:
							if lt in GCF_Object.bgc_genes[bgc]:
								lt_matching.append(lt)
						if len(lt_matching) == 1:
							ref_lt_midpoint = (GCF_Object.comp_gene_info[lt]['start'] + GCF_Object.comp_gene_info[lt]['end'])/2.0
							ref_lt_direction = (GCF_Object.comp_gene_info[lt]['direction'])
							ref_hg_order.append(ref_lt_midpoint)
							ref_hg_direction.append(ref_lt_direction)
						else:
							ref_hg_order.append(None)
							ref_hg_direction.append(None)

					best_corr = None
					try:
						list1_same_dir = []
						list2_same_dir = []
						list1_comp_dir = []
						list2_comp_dir = []
						for iter, hgval1 in enumerate(mib_hg_order):
							hgdir1 = mib_hg_direction[iter]
							hgval2 = ref_hg_order[iter]
							hgdir2 = ref_hg_direction[iter]
							if hgval1 == None or hgval2 == None: continue
							if hgdir1 == None or hgdir2 == None: continue
							if hgdir1 == hgdir2:
								list1_same_dir.append(hgval1)
								list2_same_dir.append(hgval2)
							else:
								list1_comp_dir.append(hgval1)
								list2_comp_dir.append(hgval2)

						if len(list1_same_dir) >= 3:
							corr, pval = stats.pearsonr(list1_same_dir, list2_same_dir)
							corr = abs(corr)
							if (pval < 0.1) and ((best_corr and best_corr < corr) or (not best_corr)):
								best_corr = corr
						if len(list1_comp_dir) >= 3:
							corr, pval = stats.pearsonr(list1_comp_dir, list2_comp_dir)
							corr = abs(corr)
							if (pval < 0.1) and ((best_corr and best_corr < corr) or (not best_corr)):
								best_corr = corr
					except:
						raise RuntimeError(traceback.format_exc())

					if best_corr != None and best_corr >= correlation_cutoff:
						for hg in mhg_info:
							mibig_prot_hg_tup = tuple([mhg_info[hg][2], hg])
							if not mibig_prot_hg_tup in accounted_for_tups:
								rf_handle.write('\t'.join([gcf_id, mb, hg, mhg_info[hg][2]]) + '\n')
							accounted_for_tups.add(mibig_prot_hg_tup)

	rf_handle.close()

	# Write checkpoint file for lsaBGC-AutoAnalyze.py
	checkpoint_file = outdir + 'CHECKPOINT.txt'
	checkpoint_handle = open(checkpoint_file, 'w')
	checkpoint_handle.write('lsaBGC-MIBiGMapper completed successfully!')
	checkpoint_handle.close()

	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	mapToMIBiG()
