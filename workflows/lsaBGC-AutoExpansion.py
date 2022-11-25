#!/usr/bin/env python

### Program: lsaBGC-AutoExpansion.py
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
from lsaBGC.classes.BGC import BGC
from lsaBGC import util
from decimal import Decimal
import math
from collections import defaultdict
import traceback

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_TREESTRUCTURE = lsaBGC_main_directory + '/lsaBGC/Rscripts/createNJTreeAndDefineClades.R'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-AutoExpansion.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Program to run lsaBGC-Expansion on each GCF, resolve conflicts across GCFs, and then consolidate result files when
	overlapping BGC predictions for different GCFs exist.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-g', '--gcf_listing_dir', help='Directory with GCF listing files.', required=True)
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder homolog group by sample matrix.", required=True)
	parser.add_argument('-l', '--initial_listing', type=str, help="Path to tab delimited text file for samples with three columns: (1) sample name (2) Prokka generated Genbank file (*.gbk), and (3) Prokka generated predicted-proteome file (*.faa). Please remove troublesome characters in the sample name.", required=True)
	parser.add_argument('-e', '--expansion_listing', help="Path to tab delimited file listing: (1) sample name (2) path to Prokka Genbank and (3) path to Prokka predicted proteome. This file is produced by lsaBGC-AutoProcess.py.", required=True)
	parser.add_argument('-o', '--output_directory', help="Parent output/workspace directory.", required=True)
	parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is antiSMASH.', default='antiSMASH', required=False)
	parser.add_argument('-q', '--quick_mode', action='store_true', help='Whether to run lsaBGC-Expansion in quick mode?', required=False, default=False)
	parser.add_argument('-z', '--pickle_expansion_annotation_data', help="Pickle file with serialization of annotation data in the expansion listing file.", required=False, default=None)
	parser.add_argument('-c', '--cpus', type=int, help="Total number of cpus to use.", required=False, default=1)
	parser.add_argument('-ph', '--protocore_homologs', help="File with manual listings of proto-core homolog groups.\nThis should be provided as 2 column tab-delmited file: (1) GCF id and\n(2) space delmited listing of homolog groups.", required=False, default=None)

	args = parser.parse_args()
	return args


def lsaBGC_AutoExpansion():
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
	initial_listing_file = os.path.abspath(myargs.initial_listing)
	expansion_listing_file = os.path.abspath(myargs.expansion_listing)

	try:
		assert (os.path.isdir(gcf_listing_dir) and os.path.isfile(original_orthofinder_matrix_file) and os.path.isfile(initial_listing_file) and os.path.isfile(expansion_listing_file))
	except:
		raise RuntimeError('Input directory with GCF listings does not exist, the OrthoFinder sample by homolog matrix file, input file listing initial listings file, expansions listings file, or does not exist. Exiting now ...')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Continuing in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	"""
	PARSE OPTIONAL INPUTS
	"""

	cpus = myargs.cpus
	bgc_prediction_software = myargs.bgc_prediction_software.upper()
	quick_mode = myargs.quick_mode
	pickle_expansion_annotation_data_file = myargs.pickle_expansion_annotation_data
	protocore_homologs_file = myargs.protocore_homologs

	try:
		assert (bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
	except:
		raise RuntimeError('BGC prediction software option is not a valid option.')

	if pickle_expansion_annotation_data_file != None:
		try: assert(os.path.isfile(pickle_expansion_annotation_data_file))
		except: raise RuntimeError('Issue validating pickle expansion listing file exists.')

	if protocore_homologs_file != None:
		try: assert(os.path.isfile(protocore_homologs_file))
		except: raise RuntimeError('Issue validating proto-core homologs file exists.')

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)

	version_string = util.parseVersionFromSetupPy()
	logObject.info('Running version %s' % version_string)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for easier determination of results' provenance in the future.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [gcf_listing_dir, initial_listing_file,
						expansion_listing_file, original_orthofinder_matrix_file, outdir,
						pickle_expansion_annotation_data_file, quick_mode, bgc_prediction_software,
						protocore_homologs_file, cpus]
	parameter_names = ["GCF Listings Directory", "Listing File of Prokka Annotation Files for Initial Set of Samples",
					   "Listing File of Prokka Annotation Files for Expansion/Additional Set of Samples",
					   "OrthoFinder Homolog Matrix", "Output Directory",
					   "Pickle File with Annotation Data in Expansion Listing for Quick Loading",
					   "Run in Quick Mode?", "BGC Prediction Software", "ProtoCore-Like Homolog Group Specifications",
					   "cpus"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# parse manual proto-core homolog group specifications if provided
	protocore_hgs = None
	if protocore_homologs_file != None:
		protocore_hgs = {}
		with open(protocore_homologs_file) as ophf:
			for line in ophf:
				line = line.strip('\n')
				gcf, phgs = line.split('\t')
				protocore_hgs[gcf] = phgs

	exp_outdir = outdir + 'Expansion/'
	if not os.path.isdir(exp_outdir): os.system('mkdir %s' % exp_outdir)

	gcf_expansion_results = defaultdict(dict)
	bgc_lt_evals = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
	bgcs_with_func_core_lt = set([])
	bgc_lts = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
	bgc_lt_to_hg = defaultdict(dict)
	original_gcfs = set([])
	original_gcf_samples = defaultdict(set)
	for g in os.listdir(gcf_listing_dir):
		gcf_id = g.split('.txt')[0]
		gcf_listing_file = gcf_listing_dir + g
		orthofinder_matrix_file = original_orthofinder_matrix_file
		logObject.info("Beginning expansion for GCF %s" % gcf_id)
		sys.stderr.write("Beginning expansion for GCF %s\n" % gcf_id)

		with open(gcf_listing_file) as oglf:
			for line in oglf:
				line = line.strip()
				sample, bgc_gbk_path = line.split('\t')
				expansion_flag = False
				if '_Expansion_BGC' in bgc_gbk_path:
					expansion_flag = True
				BGC_Object = BGC(bgc_gbk_path, bgc_gbk_path, expansion_flag, prediction_method=bgc_prediction_software)
				BGC_Object.parseGenbanks(comprehensive_parsing=False)
				curr_bgc_lts = set(BGC_Object.gene_information.keys())
				for lt in curr_bgc_lts:
					bgc_lt_evals[sample][gcf_id][bgc_gbk_path][lt] = -300
					if quick_mode:
						bgc_lt_evals[sample][gcf_id][bgc_gbk_path][lt] = -500
					bgc_lts[sample][gcf_id][bgc_gbk_path].add(lt)
				original_gcfs.add(line.split('\t')[1])
				original_gcf_samples[g].add(line.split('\t')[0])

		# Run lsaBGC-Expansion.py for GCF
		### TODO add additional options in expansion to auto-expansion
		gcf_exp_outdir = exp_outdir + gcf_id + '/'
		if not os.path.isdir(gcf_exp_outdir):
			cmd = ['lsaBGC-Expansion.py', '-g', gcf_listing_file, '-m', orthofinder_matrix_file, '-l',
				   initial_listing_file, '-p', bgc_prediction_software, '-e', expansion_listing_file,
				   '-o', gcf_exp_outdir, '-i', gcf_id, '-c', str(cpus)]
			if quick_mode:
				cmd += ['-q']
			if pickle_expansion_annotation_data_file:
				cmd += ['-z', pickle_expansion_annotation_data_file]
			if protocore_hgs != None:
				cmd += ['-ph', '"' + protocore_hgs[gcf_id] + '"']
			try:
				util.run_cmd(cmd, logObject, stderr=sys.stderr, stdout=sys.stdout)
			except:
				logObject.warning("lsaBGC-Expansion.py was unsuccessful, skipping over GCF %s" % gcf_id)
				sys.stderr.write("lsaBGC-Expansion.py was unsuccessful, skipping over GCF %s\n" % gcf_id)
				continue
		try:
			updated_gcf_listing_file = gcf_exp_outdir + 'GCF_Expanded.txt'
			gcf_hmm_evalues_file = gcf_exp_outdir + 'GCF_NewInstances_HMMEvalues.txt'

			assert (os.path.isfile(updated_gcf_listing_file))
			assert (os.path.isfile(gcf_hmm_evalues_file))

			gcf_expansion_results[gcf_id] = updated_gcf_listing_file

			with open(gcf_hmm_evalues_file) as oghef:
				for line in oghef:
					line = line.strip()
					bgc_gbk_path, sample, lt, hg, eval, hg_is_functionally_core = line.split('\t')
					if sample in original_gcf_samples[g]: continue
					eval = Decimal(eval)
					d = Decimal(eval + Decimal(1e-300))
					bgc_lt_evals[sample][gcf_id][bgc_gbk_path][lt] = float(max([d.log10(), -300]))
					if quick_mode:
						d = Decimal(eval + Decimal(1e-500))
						bgc_lt_evals[sample][gcf_id][bgc_gbk_path][lt] = float(max([d.log10(), -500]))
					bgc_lts[sample][gcf_id][bgc_gbk_path].add(lt)
					if hg_is_functionally_core == 'True':
						bgcs_with_func_core_lt.add(bgc_gbk_path)
					bgc_lt_to_hg[bgc_gbk_path][lt] = hg

		except Exception as e:
			raise RuntimeException(traceback.format_exc())
			logObject.warning("Key output files from lsaBGC-Expansion.py appear to be missing, skipping over GCF %s" % gcf_id)
			sys.stderr.write("Key output files lsaBGC-Expansion.py appear to be missing, skipping over GCF %s\n" % gcf_id)
			continue

	updated_gcf_listing_dir = outdir + 'Updated_GCF_Listings/'
	if not os.path.isdir(updated_gcf_listing_dir): os.system('mkdir %s' % updated_gcf_listing_dir)

	bgcs_to_discard = set([])
	# refine and filter overlapping BGCs across GCF expansions
	for sample in bgc_lts:
		for i, gcf1 in enumerate(bgc_lts[sample]):
			for j, gcf2 in enumerate(bgc_lts[sample]):
				if i >= j: continue
				for bgc1 in bgc_lts[sample][gcf1]:
					bgc1_lts = bgc_lts[sample][gcf1][bgc1]
					for bgc2 in bgc_lts[sample][gcf2]:
						bgc2_lts = bgc_lts[sample][gcf2][bgc2]
						if len(bgc1_lts.intersection(bgc2_lts)) > 0 and ( (float(len(bgc1_lts.intersection(bgc2_lts)))/float(len(bgc1_lts)) >= 0.05) or (float(len(bgc1_lts.intersection(bgc2_lts)))/float(len(bgc2_lts)) >= 0.05)):
							address = None
							if bgc1 in original_gcfs and bgc2 in original_gcfs:
								address = "neither removed"
								pass
							elif bgc1 in original_gcfs and not bgc2 in original_gcfs:
								address = bgc2 + ' removed'
								bgcs_to_discard.add(bgc2)
							elif bgc2 in original_gcfs and not bgc1 in original_gcfs:
								address = bgc1 + ' removed'
								bgcs_to_discard.add(bgc1)
							if address != None:
								logObject.info("Overlap found between BGCs %s (%s) and %s (%s), %s." % (bgc1, gcf1, bgc2, gcf2, address))
							else:
								bgc1_score = 0.0
								bgc2_score = 0.0
								for lt in bgc1_lts:
									bgc1_score += bgc_lt_evals[sample][gcf1][bgc1][lt]
								for lt in bgc2_lts:
									bgc2_score += bgc_lt_evals[sample][gcf2][bgc2][lt]

								try:
									assert(bgc1_score != bgc2_score)
								except Exception as e:
									logObject.warning("Overlapping BGCs exist with equivalent scpus for different GCFs. This should not be possible. Both instances will be discarded.")
									logObject.warning(traceback.format_exc())

								address = None
								if bgc1_score < bgc2_score:
									address = bgc2 + ' removed'
									bgcs_to_discard.add(bgc2)
								elif bgc2_score < bgc1_score:
									address = bgc1 + ' removed'
									bgcs_to_discard.add(bgc1)
								else:
									address = 'both removed'
									bgcs_to_discard.add(bgc1)
									bgcs_to_discard.add(bgc2)

								logObject.info("Overlap found between BGCs %s (%s) and %s (%s), %s." % (bgc1, gcf1, bgc2, gcf2, address))

	# create updated general listings file
	updated_listings_file = outdir + 'Sample_Annotation_Files.txt'
	updated_listings_handle = open(updated_listings_file, 'w')
	all_samples = set([])
	with open(initial_listing_file) as oilf:
		for line in oilf:
			all_samples.add(util.cleanUpSampleName(line.strip().split('\t')[0]))
			updated_listings_handle.write(line)
	with open(expansion_listing_file) as oelf:
		for line in oelf:
			if not util.cleanUpSampleName(line.strip().split('\t')[0]) in all_samples:
				all_samples.add(util.cleanUpSampleName(line.strip().split('\t')[0]))
				updated_listings_handle.write(line)
	updated_listings_handle.close()

	# further filter out entire GCF presence in samples if needed
	sample_has_bgc_with_functional_core_lt = defaultdict(lambda: defaultdict(lambda: False))
	for gcf in gcf_expansion_results:
		expanded_gcf_listing_file = gcf_expansion_results[gcf]
		with open(expanded_gcf_listing_file) as oeglf:
			for line in oeglf:
				line = line.strip()
				sample, bgc_gbk_path = line.split('\t')
				if sample in original_gcf_samples[gcf]: continue
				if (not bgc_gbk_path in bgcs_to_discard) and (bgc_gbk_path in bgcs_with_func_core_lt):
					sample_has_bgc_with_functional_core_lt[sample][gcf] = True

	# update OrthoFinder homolog group vs. sample matrix
	sample_hg_lts = defaultdict(lambda: defaultdict(set))
	for gcf in gcf_expansion_results:
		expanded_gcf_listing_file = gcf_expansion_results[gcf]
		final_expanded_gcf_listing_file = updated_gcf_listing_dir + gcf + '.txt'
		final_expanded_gcf_listing_handle = open(final_expanded_gcf_listing_file, 'w')
		with open(expanded_gcf_listing_file) as oeglf:
			for line in oeglf:
				line = line.strip()
				sample, bgc_gbk_path = line.split('\t')
				if (sample_has_bgc_with_functional_core_lt[sample][gcf] == False) and (not bgc_gbk_path in original_gcfs):
					logObject.info("GCF %s presence in sample %s disregarded, because BGC instance with functionally core homolog group was removed due to overlap with BGC from another GCF." % (gcf, sample))
					continue
				if (bgc_gbk_path in bgcs_to_discard) and (not bgc_gbk_path in original_gcfs): continue
				final_expanded_gcf_listing_handle.write(line + '\n')
				if bgc_gbk_path in original_gcfs: continue
				expansion_flag = False
				if '_Expansion_BGC' in bgc_gbk_path:
					expansion_flag = True
				BGC_Object = BGC(bgc_gbk_path, bgc_gbk_path, expansion_flag, prediction_method=bgc_prediction_software)
				BGC_Object.parseGenbanks(comprehensive_parsing=False)
				curr_bgc_lts = set(BGC_Object.gene_information.keys())
				sample = util.cleanUpSampleName(sample)
				for lt in curr_bgc_lts:
					if lt in bgc_lt_to_hg[bgc_gbk_path]:
						hg = bgc_lt_to_hg[bgc_gbk_path][lt]
						sample_hg_lts[sample][hg].add(lt)
		final_expanded_gcf_listing_handle.close()

	original_samples = []
	all_hgs = set([])
	with open(orthofinder_matrix_file) as omf:
		for i, line in enumerate(omf):
			line = line.strip('\n')
			ls = line.split('\t')
			if i == 0:
				original_samples = [util.cleanUpSampleName(x) for x in ls[1:]]
				all_samples = all_samples.union(set(original_samples))
			else:
				hg = ls[0]
				all_hgs.add(hg)
				for j, prot in enumerate(ls[1:]):
					sample_hg_lts[original_samples[j]][hg] = sample_hg_lts[original_samples[j]][hg].union(set(prot.split(', ')))

	expanded_orthofinder_matrix_file = outdir + 'Orthogroups.expanded.tsv'
	expanded_orthofinder_matrix_handle = open(expanded_orthofinder_matrix_file, 'w')

	header = [''] + [s for s in sorted(all_samples)]
	expanded_orthofinder_matrix_handle.write('\t'.join(header) + '\n')
	for hg in sorted(all_hgs):
		printlist = [hg]
		for s in sorted(all_samples):
			printlist.append(', '.join([x.strip() for x in sample_hg_lts[s][hg] if x.strip() != '']))
		expanded_orthofinder_matrix_handle.write('\t'.join(printlist) + '\n')
	expanded_orthofinder_matrix_handle.close()

	# Close logging object and exit
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	lsaBGC_AutoExpansion()
