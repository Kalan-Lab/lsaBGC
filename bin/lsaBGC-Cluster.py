#!/usr/bin/env python

### Program: lsaBGC-Cluster.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Microbiology and Immunology

import os
import sys
from time import sleep
import argparse
from lsaBGC import lsaBGC
import statistics

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Cluster.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Microbiology and Immunology

	This program will cluster BGCs found by AntiSMASH using MCL based on similarity exhibited in ortholog group presence/
	absence data. Clustering uses MCL.""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-b', '--bgc_specs_file',
						help='BGC specifications file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
	parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
	parser.add_argument('-o', '--outdir', help="Output directory.", required=True)
	parser.add_argument('-c', '--cores', type=int, help="Number of cores to use for MCL step.", required=False,
					default=1)
	parser.add_argument('-i', '--mcl_inflation', type=float, help="Inflation parameter to be used for MCL.", required=False,
					default=1.4)
	parser.add_argument('-r', '--run_inflation_tests', action='store_true',
					help="Run tests for selecting best inflation parameter for MCL analysis and exit.", default=False,
					required=False)
	args = parser.parse_args()
	return args

def bgclust():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	bgc_specs_file = os.path.abspath(myargs.bgc_specs_file)
	orthofinder_matrix = os.path.abspath(myargs.orthofinder_matrix)
	outdir = os.path.abspath(myargs.outdir) + '/'

	### vet input files quickly
	try:
		assert(os.path.isfile(orthofinder_matrix))
		assert(os.path.isfile(bgc_specs_file))
	except:
		raise RuntimeError('One or more of the input files provided, does not exist. Exiting now ...')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	"""
	PARSE OPTIONAL INPUTS
	"""

	cores = myargs.cores
	mcl_inflation = myargs.mcl_inflation
	run_inflation_tests = myargs.run_inflation_tests


	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = lsaBGC.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.
	logObject.info("Saving parameters for future provedance.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [bgc_specs_file, orthofinder_matrix, outdir, cores, mcl_inflation, run_inflation_tests]
	parameter_names = ["BGC Listing File", "OrthoFinder Orthogroups.csv File", "Output Directory", "Cores",
					   "MCL Inflation Parameter", "Whether to Run Inflation Parameter Tests"]
	lsaBGC.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Parse BGCs from Listing File
	logObject.info("Starting to process sample assemblies from listing file.")
	bgc_sample, bgc_product, bgc_genes, all_genes = lsaBGC.readInBGCGenbanksComprehensive(bgc_specs_file, logObject)
	logObject.info("Successfully parsed sample assemblies.")

	# Step 2: Parse OrthoFinder Homolog vs Sample Matrix
	logObject.info("Starting to parse OrthoFinder Homolog vs Sample information.")
	gene_to_cog, avg_nz_cog_counts = lsaBGC.parseOrthoFinderMatrix(orthofinder_matrix, all_genes, logObject)
	logObject.info("Successfully parsed sample assemblies.")

	# Step 3: Calculate overlap in homolog profiles between pairs of BGCs
	logObject.info('Calculating overlap in Ortholog Groups between BGC GBKs ...\n')
	pair_relations_file = outdir + 'bgc_pair_relationships.txt'
	prf_handle = open(pair_relations_file, 'w')
	pairwise_relations = defaultdict(lambda: defaultdict(float))
	for i, bgc1 in enumerate(bgc_cogs):
		bgc1_cogs = set([x for x in bgc_cogs[bgc1] if avg_nz_cog_counts[x] < 2])
		for j, bgc2 in enumerate(bgc_cogs):
			if i < j:
				bgc2_cogs = set([x for x in bgc_cogs[bgc2] if avg_nz_cog_counts[x] < 2])
				overlap_metric = float(len(bgc1_cogs.intersection(bgc2_cogs)))/float(min([len(bgc1_cogs), len(bgc2_cogs)]))
				overlap_metric_scaled = 100.00*overlap_metric
				if overlap_metric_scaled > 0:
					pairwise_relations[bgc1][bgc2] = overlap_metric_scaled
					pairwise_relations[bgc2][bgc1] = overlap_metric_scaled
					prf_handle.write('%s\t%s\t%f\n' % (bgc1, bgc2, overlap_metric_scaled))
	prf_handle.close()

	mci_file = outdir + 'bgc_pair_relationships.mci'
	tab_file = outdir + 'bgc_pair_relationships.tab'
	mcxload_cmd = 'mcxload -abc %s --stream-mirror -write-tab %s -o %s' % (pair_relations_file, tab_file, mci_file)
	sys.stderr.write('Converting format of pair relationship file via mxcload ...\n')
	os.system(mcxload_cmd)

	stats_file = outdir + 'GCF_details.txt'
	sf_handle = open(stats_file, 'w')
	sf_handle.write('\t'.join(['inflation parameter', 'GCF id', 'number of BGCs',
							   'samples with multiple BGCs in GCF', 'size of the SCC', 'mean number of OGs',
							   'stdev for number of OGs', 'min difference', 'max difference', 'annotations']) + '\n')

	if run_inflation_tests:
		mcl_inflation_params = [0.8, 1.4, 2, 4, 8]
		for i in mcl_inflation_params:
			mcl_out = outdir + 'mcl.' + str(i).replace('.', '_') + '.out'
			mcxdump_out = outdir + 'final_mcl.' + str(i).replace('.', '_') + '.out'
			mcl_cmd = 'mcl %s -I %f -o %s -te %d' % (mci_file, i, mcl_out, threads)
			mcxdump_cmd = 'mcxdump -icl %s -tabr %s -o %s' % (mcl_out, tab_file, mcxdump_out)

			sys.stderr.write('Running MCL with inflation parameter set to %f ...\n' % i)
			#os.system(mcl_cmd)

			sys.stderr.write('Dumping results in human-readable format ...\n')
			#os.system(mcxdump_cmd)

			clustered_bgcs = set([])
			with open(mcxdump_out) as omo:
				for j, gcf in enumerate(omo):
					gcf = gcf.strip()
					gcf_mems = gcf.split()
					if len(gcf_mems) < 2: continue
					diffs = set([])
					samp_counts = defaultdict(int)
					samp_ogs = defaultdict(set)
					annots = set([])
					for a, bgc1 in enumerate(gcf_mems):
						samp_counts[bgc1.split('_')[0]] += 1
						annots.add(bgc_annot[bgc1])
						samp_ogs[bgc1.split('_')[0]] = samp_ogs[bgc1.split('_')[0]].union(bgc_cogs[bgc1])
						clustered_bgcs.add(bgc1)
						for b, bgc2 in enumerate(gcf_mems):
							if a < b:
								diffs.add(pairwise_relations[bgc1][bgc2])
					multi_same_sample = 0
					num_ogs = []
					for si, s in enumerate(samp_counts):
						if samp_counts[s] > 1:
							multi_same_sample += 1
						if si == 0:
							scc = samp_ogs[s]
						else:
							scc = scc.intersection(samp_ogs[s])
						num_ogs.append(len(samp_ogs[s]))
					gcf_stats = [i, 'gcf_' + str(j), len(gcf_mems), multi_same_sample, len(scc), statistics.mean(num_ogs), statistics.stdev(num_ogs), min(diffs), max(diffs), '; '.join(annots)]
					sf_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
			singleton_bgcs = set([])
			for bgc in bgc_cogs:
				if not bgc in clustered_bgcs: singleton_bgcs.add(bgc)
			sf_handle.write('\t'.join([str(x) for x in ([i, 'singletons', len(singleton_bgcs)] + ['NA']*7)]) + '\n')
	else:
		mcl_out = outdir + 'mcl.' + str(i).replace('.', '_') + '.out'
		mcxdump_out = outdir + 'final_mcl.' + str(i).replace('.', '_') + '.out'
		mcl_cmd = 'mcl %s -I %f -o %s -te %d' % (mci_file, i, mcl_out, threads)
		mcxdump_cmd = 'mcxdump -icl %s -tabr %s -o %s' % (mcl_out, tab_file, mcxdump_out)

		sys.stderr.write('Running MCL with inflation parameter set to %f ...\n' % i)
		os.system(mcl_cmd)

		sys.stderr.write('Dumping results in human-readable format ...\n')
		os.system(mcxdump_cmd)

		gcf_listing_dir = outdir + 'GCF_listings/'
		if not os.path.isdir(gcf_listing_dir): os.system('mkdir %s' % gcf_listing_dir)
		with open(mcxdump_out) as omo:
			for j, gcf in enumerate(omo):
				gcf = gcf.strip()
				gcf_mems = gcf.split()
				if len(gcf_mems) < 2: continue
				outf_list = open(gcf_listing_dir + 'gcf_' + str(j+1) + '.txt', 'w')
				for bgc in gcf_mems:
					sname = '_'.join(bgc.split('_')[:-1])
					gbkpath = bgc_gbks[bgc]
					outf_list.write('%s\t%s\n' % (gbkpath, sname))
				outf_list.close()

		clustered_bgcs = set([])
		with open(mcxdump_out) as omo:
			for j, gcf in enumerate(omo):
				gcf = gcf.strip()
				gcf_mems = gcf.split()
				if len(gcf_mems) < 2: continue
				diffs = set([])
				samp_counts = defaultdict(int)
				samp_ogs = defaultdict(set)
				annots = set([])
				for a, bgc1 in enumerate(gcf_mems):
					samp_counts[bgc1.split('_')[0]] += 1
					annots.add(bgc_annot[bgc1])
					samp_ogs[bgc1.split('_')[0]] = samp_ogs[bgc1.split('_')[0]].union(bgc_cogs[bgc1])
					clustered_bgcs.add(bgc1)
					for b, bgc2 in enumerate(gcf_mems):
						if a < b:
							diffs.add(pairwise_relations[bgc1][bgc2])
				multi_same_sample = 0
				num_ogs = []
				for si, s in enumerate(samp_counts):
					if samp_counts[s] > 1:
						multi_same_sample += 1
					if si == 0:
						scc = samp_ogs[s]
					else:
						scc = scc.intersection(samp_ogs[s])
					num_ogs.append(len(samp_ogs[s]))
				gcf_stats = [i, 'gcf_' + str(j), len(gcf_mems), multi_same_sample, len(scc), statistics.mean(num_ogs), statistics.stdev(num_ogs), min(diffs), max(diffs), '; '.join(annots)]
				sf_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
		singleton_bgcs = set([])
		for bgc in bgc_cogs:
			if not bgc in clustered_bgcs: singleton_bgcs.add(bgc)
		sf_handle.write('\t'.join([str(x) for x in ([i, 'singletons', len(singleton_bgcs)] + ['NA']*7)]) + '\n')
	sf_handle.close()

if __name__ == '__main__':
	#Parse arguments.
	bgclust()



def gatherAnnotation(gbk):
	product = "NA"
	with open(gbk) as ogbk:
		for rec in SeqIO.parse(ogbk, 'genbank'):
			for feature in rec.features:
				if feature.type == "region":
					product = feature.qualifiers.get('product')[0]
	return product

