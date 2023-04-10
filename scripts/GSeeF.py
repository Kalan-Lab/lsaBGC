#!/usr/bin/env python

### Program: GSeeF.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2022, Kalan-Lab
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
import shutil
import sys
import argparse
from lsaBGC import util
import subprocess
import traceback
import math
from operator import itemgetter
from collections import defaultdict
from Bio import SeqIO
from time import sleep
import random
from ete3 import Tree
from lsaBGC import util

os.environ['OMP_NUM_THREADS'] = '4'
lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
SEED = 12345

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: GSeeF.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	*******************************************************************************************************************
	GSeeF.py creates visualizations of GCF distributions across a species tree given BiG-SCAPE or lsaBGC-Clustering
	results. 
	
	E.g.: GSeeF.py -b BiG-SCAPE_Results/ -a AntiSMASH_Input_Dir_for_BiG-SCAPE/ -o GSeeF_Results/
	E.g.: GSeeF.py -g lsaBGC-Cluster_GCF_Listings_Results/ -l Sample_Annotation_Listings.txt -o GSeeF_Results/
	 
	*******************************************************************************************************************
	OVERVIEW OF STEPS TAKEN:
	    - Step 1: Run Gene Calling on Genomes. 
		- Step 2: Run GToTree to produce species tree.
		- Step 3: Create final visuals in PDF format (using ggtree) + annotation tracks for iTol
		
	*******************************************************************************************************************
	CITATION:
	- Please make sure to cite GToTree - Lee 2019 (this is how phylogenies are made)
	- lsaBGC - Salamzade et al. 2022
	- BiG_SCAPE - Navarro-Munoz & Selem-Mojica et al. 2020
	- iTol - Letunic & Bork 2021
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

	parser.add_argument('-b', '--bigscape_results_directory', help='BiG-SCAPE Results directory.', required=False, default=None)
	parser.add_argument('-a', '--antismash_results_directory', help='The directory with AntiSMASH results used as input for BiG-SCAPE. Should have folders matching individual genomes which contain AntiSMASH results for each genome.', required=False, default=None)
	parser.add_argument('-l', '--lsabgc_sample_annotation_file', help='lsaBGC-type sample annotation tab-delimited text file. First column is sample name, second is path to full genome GenBank, third is path to full genome proteome.', required=False, default=None)
	parser.add_argument('-g', '--lsabgc_gcf_listings_directory', help='lsaBGC-Cluster/Auto-Expansion GCF listings directory.')
	parser.add_argument('-o', '--output_directory', help='Parent output/workspace directory.', required=True)
	parser.add_argument('-s', '--species_tree', help='Provide species tree in Newick format. If not provided will run GToTree to generate species phylogeny.', required=False, default=None)
	parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is antiSMASH. Will only work for GECCO and DeepBGC if\n-l and -g arguments are used (a.k.a. in lsaBGC mode) to provide paths to genomes and BGC predictions.', default='antismash', required=False)
	parser.add_argument('-m', '--max_gcfs', type=int, help='The maximum number of GCFs to consider (will prioritize inclusion of more common GCFs). (Default is 50).', default=50, required=False)
	parser.add_argument('-gtm', '--gtotree_model', help="SCG model for secondary GToTree analysis and what would be used for dereplication. (Default is Bacteria)", default='Bacteria', required=False)
	parser.add_argument('-c', '--cpus', type=int, help="Total number of cpus/threads. Note, this is the total number of\nthreads to use. (Default is 1)", required=False, default=4)

	args = parser.parse_args()
	return args


def GSeeF():
	myargs = create_parser()

	bigscape_results_directory = myargs.bigscape_results_directory
	antismash_results_directory = myargs.antismash_results_directory
	lsabgc_sample_annotation_file = myargs.lsabgc_sample_annotation_file
	lsabgc_gcf_listings_directory = myargs.lsabgc_gcf_listings_directory
	species_tree_file = myargs.species_tree
	outdir = os.path.abspath(myargs.output_directory) + '/'
	cpus = myargs.cpus
	max_gcfs = myargs.max_gcfs
	gtotree_model = myargs.gtotree_model
	bgc_prediction_software = myargs.bgc_prediction_software.upper()

	checkModelIsValid(gtotree_model)

	try:
		assert (bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
	except:
		sys.stderr.write('BGC prediction software option is not a valid option.\n')
		sys.exit(1)

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory already exists! Overwriting in 5 seconds, but only where needed will use checkpoint files to skip certain steps...\n")
		sleep(5)
	else:
		util.setupReadyDirectory([outdir])

	final_results_dir = outdir + 'Final_Results/'
	if not os.path.isdir(final_results_dir):
		os.system('mkdir %s' % final_results_dir)

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)
	logObject.info("Saving command for future records.")
	parameters_file = outdir + 'Command_Issued.txt'
	parameters_handle = open(parameters_file, 'a+')
	parameters_handle.write(' '.join(sys.argv) + '\n')
	parameters_handle.close()
	sys.stdout.write("Appending command issued for future records to: %s\n" % parameters_file)
	sys.stdout.write("Logging more details at: %s\n" % parameters_file)
	parallel_jobs_4cpu = max(math.floor(cpus / 4), 1)

	mode = None
	if bigscape_results_directory != None and antismash_results_directory != None and \
			os.path.isdir(bigscape_results_directory) and os.path.isdir(antismash_results_directory):
		bigscape_results_directory = os.path.abspath(bigscape_results_directory) + '/'
		antismash_results_directory = os.path.abspath(antismash_results_directory) + '/'
		mode = 'BiG-SCAPE'
		logObject.info('\n---------------------\nRunning in BiG-SCAPE mode!\n---------------------\n')
		sys.stdout.write('---------------------\nRunning in BiG-SCAPE mode!\n---------------------\n')
	elif lsabgc_sample_annotation_file != None and lsabgc_gcf_listings_directory != None and \
			os.path.isfile(lsabgc_sample_annotation_file) and os.path.isdir(lsabgc_gcf_listings_directory):

		if not os.path.abspath(lsabgc_gcf_listings_directory).endswith('GCF_Listings'):
			sys.stderr.write('lsaBGC GCF Listings directory does not end with GCF_Listings/. Maybe you provided the parent directory instead?\n')
			logObject.error('lsaBGC GCF Listings directory does not end with GCF_Listings/. Maybe you provided the parent directory instead?')
			sys.exit(1)
		lsabgc_sample_annotation_file = os.path.abspath(lsabgc_sample_annotation_file)
		lsabgc_gcf_listings_directory = os.path.abspath(lsabgc_gcf_listings_directory) + '/'
		mode = 'lsaBGC'
		logObject.info('\n---------------------\nRunning in lsaBGC mode!\n---------------------\n')
		sys.stdout.write('---------------------\nRunning in lsaBGC mode!\n---------------------\n')
	else:
		sys.stderr.write('Difficulty figuring out whether you are trying to run with BiG-SCAPE or lsaBGC-Cluster results. Expected inputs were unable to be validated either provide arguments for -b/-a parameters for BiG-SCAPE mode or -l/-s for lsaBGC mode.\n')
		logObject.error('Difficulty figuring out whether you are trying to run with BiG-SCAPE or lsaBGC-Cluster results. Expected inputs were unable to be validated either provide arguments for -b/-a parameters for BiG-SCAPE mode or -l/-s for lsaBGC mode.')
		sys.exit(1)

	# Step 1: Get Genome-Wide Proteomes for GToTree Phylogeny Construction
	logObject.info('\n---------------------\nStep 1\n---------------------\nGathering predicted proteomes for all samples to use for phylogeny construction.')
	sys.stdout.write('---------------------\nStep 1\n---------------------\nGathering predicted prtoemes for all samples to use for phylogeny construction!\n')

	proteome_listing_file = outdir + 'Proteomes_Listing.txt'
	proteome_listing_handle = open(proteome_listing_file, 'w')
	bigscape_bgc_to_genome = {}
	bigscape_bgc_to_annotation = {}
	all_samples = set([])
	if mode == 'BiG-SCAPE':
		proteomes_dir = outdir + 'Proteomes/'
		if os.path.isdir(proteomes_dir):
			shutil.rmtree(proteomes_dir)
		os.system('mkdir %s' % proteomes_dir)

		for sd in os.listdir(antismash_results_directory):
			subdir = antismash_results_directory + sd + '/'
			if not os.path.isdir(subdir): continue
			full_genbank_file = None
			count_possible_full_genome_genbanks = 0
			for f in os.listdir(subdir):
				if f.endswith('.gbk') and not 'region' in f:
					full_genbank_file = subdir + f
					count_possible_full_genome_genbanks += 1
			if count_possible_full_genome_genbanks != 1 or full_genbank_file == None or not os.path.isfile(full_genbank_file): continue
			all_samples.add(sd)
			sample_proteome_file = proteomes_dir + sd + '.faa'
			sample_proteome_handle = open(sample_proteome_file, 'w')
			with open(full_genbank_file) as ofgf:
				for rec in SeqIO.parse(ofgf, 'genbank'):
					for feat in rec.features:
						if feat.type == 'CDS':
							prot_lt = feat.qualifiers.get('locus_tag')[0]
							prot_seq = str(feat.qualifiers.get('translation')[0]).replace('*', '')
							sample_proteome_handle.write('>' + prot_lt + '\n' + prot_seq + '\n')
			sample_proteome_handle.close()
			proteome_listing_handle.write(sample_proteome_file + '\n')

			for f in os.listdir(subdir):
				if f.endswith('.gbk') and 'region' in f:
					bgc_id = f.split('.gbk')[0]
					bigscape_bgc_to_genome[bgc_id] = sd
					bigscape_bgc_to_annotation[bgc_id] = parseAntiSMASHGBKForFunction(subdir + f, logObject)

	elif mode == 'lsaBGC':
		with open(lsabgc_sample_annotation_file) as olsaf:
			for line in olsaf:
				line = line.strip()
				sample, gbk, faa = line.split('\t')
				all_samples.add(sample)
				proteome_listing_handle.write(faa + '\n')
	proteome_listing_handle.close()

	# Step 2: Parse Clustering Results
	logObject.info('\n---------------------\nStep 2\n---------------------\nParsing clustering results and gathering annotations from antiSMASH BGC Genbanks.')
	sys.stdout.write('---------------------\nStep 2\n---------------------\nParsing clustering results and gathering annotations from antiSMASH BGC Genbanks.\n')

	gcf_samples = defaultdict(set)
	gcf_sample_annotations = defaultdict(lambda: defaultdict(set))
	if mode == 'BiG-SCAPE':
		for dirpath, dirnames, files in os.walk(bigscape_results_directory):
			for filename in files:
				if filename.endswith(".tsv") and "clustering" in filename:
					cluster_tsv = os.path.join(dirpath, filename)
					with open(cluster_tsv) as oct:
						for i, line in enumerate(oct):
							if i == 0: continue
							line = line.strip()
							bgc_id, gcf_id = line.split('\t')
							gcf_id = 'GCF_' + gcf_id
							if not bgc_id in bigscape_bgc_to_genome: continue
							bgc_annotation = bigscape_bgc_to_annotation[bgc_id]
							sample = bigscape_bgc_to_genome[bgc_id]
							gcf_samples[gcf_id].add(sample)
							gcf_sample_annotations[gcf_id][sample].add(bgc_annotation)

	elif mode == 'lsaBGC':
		for gcf_tsv in os.listdir(lsabgc_gcf_listings_directory):
			gcf_tsv_file = lsabgc_gcf_listings_directory + gcf_tsv
			gcf_id = gcf_tsv.split('.tsv')[0]
			try:
				with open(gcf_tsv_file) as ogtf:
					for line in ogtf:
						line = line.strip()
						sample, bgc_gbk = line.split('\t')
						gcf_samples[gcf_id].add(sample)
						if bgc_prediction_software == 'ANTISMASH':
							bgc_annotation = parseAntiSMASHGBKForFunction(bgc_gbk, logObject)
						elif bgc_prediction_software == 'DEEPBGC':
							bgc_annotation = parseDeepBGCGBKForFunction(bgc_gbk, logObject)
						elif bgc_prediction_software == 'GECCO':
							bgc_annotation = parseGECCOGBKForFunction(bgc_gbk, logObject)
						gcf_sample_annotations[gcf_id][sample].add(bgc_annotation)
			except:
				logObject.warning('Difficulty reading file %s in GCF Listings directory.' % gcf_tsv_file)

	# Step 3: Run GToTree Phylogeny Construction and Perform Midpoint Rooting
	logObject.info('\n---------------------\nStep 3\n---------------------\nGenerating species phylogeny using GToTree.')
	sys.stdout.write('---------------------\nStep 3\n---------------------\nGenerating species phylogeny using GToTree.\n')

	if species_tree_file != None:
		species_tree_file = os.path.abspath(species_tree_file)
		try:
			assert(os.path.isfile(species_tree_file) and util.is_newick(species_tree_file))
		except:
			sys.stderr.write('Error: Unable to find user provided species tree or uanble to validate it is in Newick format. Note, this isn\'t required, if you exclude providing this argument, GToTree will be run.\n')
			sys.exit(1)
	else:
		gtotree_outdir = outdir + 'GToTree_output/'
		species_tree_file = gtotree_outdir + 'GToTree_output.tre'
		if not os.path.isfile(species_tree_file):
			os.system('rm -rf %s' % gtotree_outdir)
			gtotree_cmd = ['GToTree', '-A', proteome_listing_file, '-H', gtotree_model, '-n', '4', '-j',
						   str(parallel_jobs_4cpu), '-M', '4', '-o', gtotree_outdir]
			runCmdViaSubprocess(gtotree_cmd, logObject, check_files=[species_tree_file])

	rooted_species_tree_file = final_results_dir + 'Species_Phylogeny.Midpoint_Rooted.tre'
	try:
		t = Tree(species_tree_file)
		R = t.get_midpoint_outgroup()
		t.set_outgroup(R)
		t.write(format=1, outfile=rooted_species_tree_file)
	except:
		sys.stderr.write('Error: Issue midpoint rooting the species tree.')
		sys.exit(1)

	# Determine single annotation category/class per sample-GCF pairing and brew colors palette
	all_annotation_classes = set(['unknown', 'multi-type'])
	gs_annots = defaultdict(dict)
	gcf_annots = {}
	for gcf in gcf_sample_annotations:
		all_gcf_annots = defaultdict(int)
		for sample in gcf_sample_annotations[gcf]:
			annotations = gcf_sample_annotations[gcf][sample]
			annot = 'unknown'
			if len(annotations) == 1:
				annot = list(annotations)[0]
			elif len(annotations) == 2 and 'NRPS' in annotations and 'NRPS-like' in annotations:
				annot = 'NRPS'
			elif len(annotations) >= 2:
				annot = 'multi-type'
			gs_annots[gcf][sample] = annot
			all_annotation_classes.add(annot)
			all_gcf_annots[annot] += 1

	colorBrewScript = lsaBGC_main_directory + 'lsaBGC/Rscripts/brewColors.R'
	colorsbrewed_file = outdir + 'colors_brewed.txt'
	colorbrew_cmd = ['Rscript', colorBrewScript, str(len(all_annotation_classes)-2), colorsbrewed_file]
	runCmdViaSubprocess(colorbrew_cmd, logObject, check_files=[colorsbrewed_file])

	colors = []
	with open(colorsbrewed_file) as ocf:
		colors = [x.strip() for x in ocf.readlines()]
	random.Random(SEED).shuffle(colors)

	annotation_colors = {}
	ci = 0
	for a in sorted(all_annotation_classes):
		if a == 'multi-type':
			annotation_colors[a] = '#000000'
		elif a == 'unknown':
			annotation_colors[a] = '#a2a4a6'
		else:
			annotation_colors[a] = colors[ci]
			ci += 1

	gcfs_order_index = {}
	for gcf in gcf_samples:
		gcfs_order_index[gcf] = len(gcf_samples[gcf])

	gcf_order = []
	#gcf_colors = []
	for gcf_tup in sorted(gcfs_order_index.items(), key=itemgetter(1), reverse=True):
		gcf_order.append(gcf_tup[0])
		#gcf_colors.append(annotation_colors[gcf_annots[gcf_tup[0]]])

	# Step 4: Create track files for Rscript and iTol
	logObject.info('\n---------------------\nStep 4\n---------------------\nCreating files for plotting.')
	sys.stdout.write('---------------------\nStep 4\n---------------------\nCreating files for plotting.\n')

	## get samples in tree to consider when creating tracks
	tree_samples = set([])
	t = Tree(rooted_species_tree_file)
	for node in t.traverse("postorder"):
		if node.is_leaf: tree_samples.add(node.name)

	itol_track_file = final_results_dir + 'GCF_Heatmap.iTol.txt'
	itol_track_handle = open(itol_track_file, 'w')
	gseef_track_file = outdir + 'gseef_track_input.txt'
	gseef_track_handle = open(gseef_track_file, 'w')

	itol_track_handle.write('DATASET_DOMAINS\n')
	itol_track_handle.write('SEPARATOR TAB\n')
	itol_track_handle.write('DATASET_LABEL\tGCFs\n')
	itol_track_handle.write('COLOR\t#000000\n')
	itol_track_handle.write('DATA\n')

	gseef_track_handle.write('gcf\tgcf_index\tlabel\tannotation\tcolors\n')

	tot_gcfs = len(gcf_order)
	for sample in sorted(all_samples):
		if not sample in tree_samples: continue
		printlist = [sample, str(tot_gcfs)]
		for i, gcf in enumerate(gcf_order):
			if i >= max_gcfs: continue
			col_array = []
			col_array.append('RE')
			col_array.append(str(i))
			col_array.append(str(i+1))
			if sample in gs_annots[gcf]:
				col_array.append(annotation_colors[gs_annots[gcf][sample]])
				gseef_track_handle.write('\t'.join([gcf, str(i), sample, gs_annots[gcf][sample], '"' + annotation_colors[gs_annots[gcf][sample]] + '"']) + '\n')
			else:
				gseef_track_handle.write('\t'.join([gcf, str(i), sample, 'absent', '"#FFFFFF"']) + '\n')
				col_array.append('#FFFFFF')
			col_array.append(gcf)
			printlist.append('|'.join(col_array))
		itol_track_handle.write('\t'.join(printlist) + '\n')
	itol_track_handle.close()
	gseef_track_handle.close()

	# Step 5: Run GSeeF.R Rscript
	logObject.info('\n---------------------\nStep 5\n---------------------\nGenerating final static PDF visualization and iTol track file.')
	sys.stdout.write('---------------------\nStep 5\n---------------------\nGenerating final static PDF visualization and iTol track file.\n')

	gSeefScript = lsaBGC_main_directory + 'lsaBGC/Rscripts/gSeef.R'
	resulting_png = final_results_dir + 'Phylogenetic_Heatmap.png'
	label_resulting_png = final_results_dir + 'Annotation_Legend.png'
	gseef_cmd = ['Rscript', gSeefScript, rooted_species_tree_file, gseef_track_file, resulting_png, label_resulting_png]
	runCmdViaSubprocess(gseef_cmd, logObject, check_files=[resulting_png, label_resulting_png])

	# Close logging object and exit
	logObject.info('GSeeF completed! Check out the major results in the folder: %s' % final_results_dir)
	sys.stdout.write('GSeeF completed! Check out the major results in the folder:\n%s\n' % final_results_dir)
	util.closeLoggerObject(logObject)
	sys.exit(0)

def parseGECCOGBKForFunction(bgc_gbk, logObject):
	try:
		rec = SeqIO.read(bgc_gbk, 'genbank')
		product = 'unknown'
		try:
			product = rec.annotations['structured_comment']['GECCO-Data']['biosyn_class']
		except:
			try:
				product = rec.annotations['structured_comment']['GECCO-Data']['cluster_type']
			except:
				pass
		if product == 'Unknown':
			product = 'unknown'
		return(product)
	except:
		logObject.error('Issues parsing BGC Genbank %s' % bgc_gbk)
		raise RuntimeError()

def parseDeepBGCGBKForFunction(bgc_gbk, logObject):
	product = 'unknown'
	try:
		products = set([])
		with open(bgc_gbk) as obg:
			for rec in SeqIO.parse(obg, 'genbank'):
				for feat in rec.features:
					if not feat.type == 'cluster': continue
					try:
						products.add(feat.qualifiers.get('product_class')[0])
					except:
						pass
		if len(products) == 1:
			product = list(products)[0]
		elif len(products) >= 2:
			product = 'multi-type'
	except:
		logObject.error('Issues parsing BGC Genbank %s' % bgc_gbk)
		raise RuntimeError()
	return (product)

def parseAntiSMASHGBKForFunction(bgc_gbk, logObject):
	product = 'unknown'
	try:
		products = set([])
		with open(bgc_gbk) as obg:
			for rec in SeqIO.parse(obg, 'genbank'):
				for feat in rec.features:
					if feat.type == 'protocluster':
						try:
							products.add(feat.qualifiers.get('product')[0])
						except:
							pass
		if len(products) == 1:
			product = list(products)[0]
		elif len(products) == 2 and 'NRPS-like' in products and 'NRPS' in products:
				product = 'NRPS'
		elif len(products) >= 2:
			product = 'multi-type'
	except:
		logObject.error('Issues parsing BGC Genbank %s' % bgc_gbk)
		raise RuntimeError()
	return(product)

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

def checkModelIsValid(gtotree_model):
	try:
		assert (gtotree_model in set(['Actinobacteria', 'Alphaproteobacteria', 'Bacteria', 'Archaea',
									  'Bacteria_and_Archaea', 'Bacteroidetes', 'Betaproteobacteria', 'Chlamydiae',
									  'Cyanobacteria', 'Epsilonproteobacteria', 'Firmicutes', 'Gammaproteobacteria',
									  'Proteobacteria', 'Tenericutes', 'Universal_Hug_et_al']))
	except:
		raise RuntimeError('Model for GToTree specified is not a valid model!')


if __name__ == '__main__':
	GSeeF()
