import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import logging
import subprocess
import statistics
from operator import itemgetter
from collections import defaultdict
import random

def createBGCGenbank(full_genbank_file, new_genbank_file, scaffold, start_coord, end_coord):
	"""
	Function to prune full genome-sized Genbank for only features in BGC of interest.
	"""
	try:
		ngf_handle = open(new_genbank_file, 'w')

		with open(full_genbank_file) as ogbk:
			for rec in recs:
				if not rec.id == scaffold: continue
				original_seq = str(rec.seq)
				filtered_seq = ""
				if end_coord == len(original_seq):
					filtered_seq = original_seq[start_coord - 1:]
				else:
					filtered_seq = original_seq[start_coord - 1:end_coord]

				new_seq_object = Seq(filtered_seq)

				updated_rec = copy.deepcopy(rec)
				updated_rec.seq = new_seq_object

				updated_features = []
				for feature in rec.features:
					start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])

					feature_coords = set(range(start, end + 1))
					if len(feature_coords.intersection(pruned_coords)) > 0:
						updated_start = start - start_coord + 1
						updated_end = end - start_coord + 1

						if end > end_coord: updated_end = end_coord - start_coord + 1
						if start < start_coord: updated_start = 1

						strand = 1
						if '(-)' in str(feature.location):
							strand = -1

						updated_location = FeatureLocation(updated_start - 1, updated_end, strand=strand)
						if feature.type == 'CDS':
							print(refined_genbank_file + '\t' + feature.qualifiers.get('locus_tag')[0] + '\t' + str(
								strand) + '\t' + str(updated_start) + '\t' + str(updated_end) + '\t' + str(
								len(filtered_seq[updated_start - 1:updated_end]) / 3.0))
						updated_feature = copy.deepcopy(feature)

						updated_feature.location = updated_location
						updated_features.append(updated_feature)
				updated_rec.features = updated_features
				SeqIO.write(updated_rec, ngf_handle, 'genbank')
		ngf_handle.close()
	except Exception as e:
		raise RuntimeError(traceback.format_exc())

def parseGenbankAndFindBoundaryGenes(sample_genbank, distance_to_scaffold_boundary=2000):
	"""
	Function to parse Genbanks from Prokka and return a dictionary of genes per scaffold, gene to scaffold, and a
	set of genes which lie on the boundary of scaffolds.

	:param sample_genbank: Prokka generated Genbank file.
	:param distance_to_scaffold_boundary: Distance to scaffold edge considered as boundary.

	:return gene_to_scaffold: Dictionary mapping each gene's locus tag to the scaffold it is found on.
	:return scaffold_genes: Dictionary with keys as scaffolds and values as a set of genes found on that scaffold.
	:return boundary_genes: Set of gene locus tag ids which are found within proximity to scaffold edges.
	"""

	gene_location = {}
	scaffold_genes = defaultdict(set)
	boundary_genes = set([])
	gene_id_to_order = defaultdict(dict)
	gene_order_to_id = defaultdict(dict)

	with open(sample_genbank) as osg:
		for rec in SeqIO.parse(osg, 'genbank'):
			scaffold = rec.id
			scaffold_length = len(str(rec.seq))
			boundary_ranges = set(range(1, distance_to_scaffold_boundary+1)).union(set(range(scaffold_length-distance_to_scaffold_boundary, scaffold_length + 1)))
			gene_starts = []
			for feature in rec.features:
				if not feature.type == 'CDS': continue
				locus_tag = feature.qualifiers.get('locus_tag')[0]
				start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
				end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])

				gene_location[locus_tag] = {'scaffold': scaffold, 'start': start, 'end:': end}
				scaffold_genes[scaffold].add(gene)

				gene_range = set(start, end+1)
				if len(gene_range.intersection(boundary_ranges)) > 0:
					bounary_genes.add(locus_tag)

				gene_starts.append([locus_tag, start])

			for i, g in enumerate(sorted(gene_starts, key=itemgetter(1))):
				gene_id_to_order[scaffold][g[0]] = i
				gene_order_to_id[scaffold][i] = g[0]

	return([gene_location, scaffold_genes, boundary_genes, gene_id_to_order, gene_order_to_id])

def calculateMashPairwiseDifferences(fasta_listing_file, outdir, name, sketch_size, cores, logObject):
	"""

	"""
	mash_db = outdir + name
	fastas = []
	fasta_to_name = {}
	try:
		with open(fasta_listing_file) as oflf:
			for line in oflf:
				line = line.strip()
				ls = line.split('\t')
				fastas.append(ls[1])
				fasta_to_name[ls[1]] = ls[0]
	except:
		error_message = "Had issues reading the FASTA listing file %s" % fasta_listing_file
		logObject.error(error_message)
		raise RuntimeError(error_message)

	# create mash database (using mash sketch)
	mash_sketch_cmd = ['mash', 'sketch', '-p', str(cores), '-s', str(sketch_size), '-o', mash_db] + fastas
	logObject.info('Running mash sketch with the following command: %s' % ' '.join(mash_sketch_cmd))
	try:
		subprocess.call(' '.join(mash_sketch_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(mash_sketch_cmd))
	except:
		error_message = 'Had an issue running: %s' % ' '.join(mash_sketch_cmd)
		logObject.error(error_message)
		raise RuntimeError(error_message)
	mash_db = mash_db + '.msh'

	try:
		assert (os.path.isfile(mash_db))
	except:
		error_message = "Had issue validating that MASH sketching worked properly, couldn't find: %s" % mash_db
		logObject.error(error_message)
		raise RuntimeError(error_message)

	# run mash distance estimation
	mash_dist_cmd = ['mash', 'dist', '-s', str(sketch_size), '-p', str(cores), mash_db, mash_db, '>',
					 outdir + name + '.out']
	logObject.info('Running mash dist with the following command: %s' % ' '.join(mash_dist_cmd))
	try:
		subprocess.call(' '.join(mash_dist_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(mash_dist_cmd))
	except:
		error_message = 'Had an issue running: %s' % ' '.join(mash_dist_cmd)
		logObject.error(error_message)
		raise RuntimeError(error_message)

	pairwise_distances = defaultdict(lambda: defaultdict(float))
	try:
		with open(outdir + name + '.out') as of:
			for line in of:
				line = line.strip()
				ls = line.split('\t')
				f1, f2, dist = ls[:3]
				dist = float(dist)
				n1 = fasta_to_name[f1]
				n2 = fasta_to_name[f2]
				pairwise_distances[n1][n2] = dist
	except:
		error_message = 'Had issues reading the output of MASH dist anlaysis in: %s' % outdir + name + '.out'
		logObject.error(error_message)
		raise RuntimeError(error_message)
	return pairwise_distances

def parseOrthoFinderMatrix(orthofinder_matrix_file, relevant_gene_lts):
	"""
	Function to parse and return information from OrthoFinderV2 de novo homolog group identification.

	:param orthofinder_matrix: OrthoFinderV2 matrix Orthogroups.csv file
	:param relevant_gene_lts: set of all the relevant gene locus tag identifiers found in BGC Genbanks

	:return gene_to_hg: dictionary mapping gene locus tags to homolog group
	:return hg_genes: dictionary with set of gene locus tags for each homolog group
	:return hg_median_gene_counts: median copy count for each homolog group
	:return hg_multicopy_proportion: proportion of samples with homolog group which have multiple (paralogous) genes in the homolog group.
	"""
	gene_to_hg = {}
	hg_genes = defaultdict(set)
	hg_multicopy_proportion = defaultdict(lambda: 'NA')
	hg_median_gene_counts = defaultdict(lambda: 'NA')
	with open(orthofinder_matrix_file) as ofm:
		for i, line in enumerate(ofm):
			if i == 0: continue
			line = line.strip('\n')
			ls = line.split('\t')
			hg = ls[0]
			flag_in_bgc = False
			for sgs in ls[1:]:
				for g in sgs.split(', '):
					if g in relevant_gene_lts:
						flag_in_bgc = True
						gene_to_hg[g] = hg
						hg_genes[hg].add(g)
			if flag_in_bgc:
				gene_counts = []
				for sgs in ls[1:]:
					gene_counts.append(len([x for x in sgs.split(', ') if len(x.split('_')[0]) == 3]))

				hg_multicopy_proportion[hg] = float(sum([1 for x in gene_counts if x > 1])) / sum([1 for x in gene_counts if x > 0])
				hg_median_gene_counts[hg] = statistics.median(gene_counts)

	return ([gene_to_hg, hg_genes, hg_median_gene_counts, hg_multicopy_proportion])

def multiProcess(input):
	"""
	Genralizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list), with last item in list corresponding to a logging object handle for logging
	progress.
	"""
	input_cmd = input[:-1]
	logObject = input[-1]
	logObject.info('Running the following command: %s' % ' '.join(input_cmd))
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(input_cmd))
	except:
		logObject.warning('Had an issue running: %s' % ' '.join(input_cmd))

def is_fasta(fasta):
	"""
	Function to validate if FASTA file is correctly formatted.
	"""
	try:
		with open(fasta) as of:
			SeqIO.parse(of, 'fasta')
		return True
	except:
		return False

def is_genbank(gbk):
	"""
	Function to check in Genbank file is correctly formatted.
	"""
	try:
		with open(gbk) as of:
			SeqIO.parse(of, 'genbank')
		return True
	except:
		return False

def createLoggerObject(log_file):
	"""
	Function which creates logger object.
	:param log_file: path to log file.
	:return: logging logger object.
	"""
	logger = logging.getLogger('task_logger')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	# create console handler with a higher log level
	ch = logging.StreamHandler()
	ch.setLevel(logging.ERROR)
	# create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)
	# add the handlers to the logger
	logger.addHandler(fh)
	logger.addHandler(ch)

	# q = multiprocessing.Queue()
	# queue_listner = QueueListener(q, ch)
	# queue_listner.start()

	# logger.handlers[0].stream = sys.stderr
	return logger

def closeLoggerObject(logObject):
	"""
	Function which closes/terminates loggerObject.
	:param logObject: logging logger object to close
	"""
	handlers = logObject.handlers[:]
	for handler in handlers:
		handler.close()
		logObject.removeHandler(handler)

def logParameters(parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to std.stderr
	"""
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		sys.stderr.write(pn + ': ' + str(pv) + '\n')

def logParametersToFile(parameter_file, parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to text file.
	"""
	parameter_handle = open(parameter_file, 'w')
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		parameter_handle.write(pn + ': ' + str(pv) + '\n')
	parameter_handle.close()

def logParametersToObject(logObject, parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to text file.
	"""
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		logObject.info(pn + ': ' + str(pv))




