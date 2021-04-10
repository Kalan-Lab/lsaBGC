import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from ete3 import Tree
import logging
import subprocess
import multiprocessing
from collections import defaultdict
from operator import itemgetter
import itertools
import statistics
import random
from scipy.stats import f_oneway
import pysam
import random

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_BGSEE = lsaBGC_main_directory + '/lsaBGC/bgSee.R'
RSCRIPT_FOR_PLOTTING = lsaBGC_main_directory + '/lsaBGC/plotCogConservation.R'
RSCRIPT_FOR_CLUSTER_ASSESSMENT_PLOTTING = lsaBGC_main_directory + '/lsaBGC/plotParameterImpactsOnGCF.R'
RSCRIPT_FOR_TAJIMA = lsaBGC_main_directory + '/lsaBGC/calculateTajimasD.R'
FLANK_SIZE = 500

def calculateMashPairwiseDifferences(fasta_listing_file, outdir, name, sketch_size, cores, logObject):
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


def readInBGCGenbanksPerGCF(gcf_specs_file, logObject, comprehensive_parsing=True):
	sample_index = defaultdict(int)
	bgc_gbk = {}
	bgc_sample = {}
	sample_bgcs = defaultdict(set)
	bgc_genes = {}
	all_genes = set([])
	comp_gene_info = {}
	with open(gcf_specs_file) as obsf:
		for i, line in enumerate(obsf):
			line = line.strip()
			try:
				assert (len(line.split('\t')) == 2)
			except:
				logObject.error(
					"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
							i + 1))
				raise RuntimeError(
					"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
							i + 1))
			sample, gbk = line.split('\t')
			try:
				assert (is_genbank(gbk))
				bgc_id = sample
				if sample_index[sample] > 0:
					bgc_id = sample + '_' + str(sample_index[sample] + 1)
				sample_index[sample] += 1

				bgc_genes_full, bgc_info_full = parseGenbanks(gbk, bgc_id, comprehensive_parsing=comprehensive_parsing)

				comp_gene_info.update(bgc_genes_full)
				bgc_genes[bgc_id] = set(bgc_genes_full.keys())
				all_genes = all_genes.union(bgc_genes[bgc_id])
				bgc_gbk[bgc_id] = gbk
				bgc_sample[bgc_id] = sample
				sample_bgcs[sample].add(bgc_id)

				logObject.info("Incorporating genbank %s for sample %s into analysis." % (gbk, sample))
			except:
				logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis.")
				raise RuntimeWarning("Unable to validate %s as Genbank. Skipping ...")

	return [bgc_gbk, bgc_genes, comp_gene_info, all_genes, bgc_sample, sample_bgcs]

def visualizeGCFViaR(gggenes_track_file, heatmap_track_file, phylogeny_file, result_pdf_file, bgc_genes, gene_to_cog, cog_to_color, comp_gene_info, logObject):

	if not os.path.isfile(gggenes_track_file) or not os.path.isfile(heatmap_track_file):
		gggenes_track_handle = open(gggenes_track_file, 'w')
		heatmap_track_handle = open(heatmap_track_file, 'w')
		logObject.info("Writing gggenes input file to: %s" % gggenes_track_file)
		logObject.info("Writing heatmap input file to: %s" % heatmap_track_file)
		# write header for track files
		gggenes_track_handle.write('label\tgene\tstart\tend\tforward\tog\tog_color\n')
		heatmap_track_handle.write('label\tog\tog_presence\tog_count\n')

		ref_cog_directions = {}

		bgc_gene_counts = defaultdict(int)
		for bgc in bgc_genes:
			bgc_gene_counts[bgc] = len(bgc_genes[bgc])

		tree_obj = Tree(phylogeny_file)
		bgc_weights = defaultdict(int)
		for leaf in tree_obj:
			bgc_weights[str(leaf).strip('\n').lstrip('-')] += 1

		bgc_cog_presence = defaultdict(lambda: defaultdict(lambda: 'Absent'))
		cog_counts = defaultdict(int)
		for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
			bgc = item[0]
			curr_bgc_genes = bgc_genes[bgc]
			last_gene_end = max([comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
			printlist = []
			cog_directions = {}
			cog_lengths = defaultdict(list)
			for lt in curr_bgc_genes:
				ginfo = comp_gene_info[lt]
				cog = 'singleton'
				if lt in gene_to_cog:
					cog = gene_to_cog[lt]

				gstart = ginfo['start']
				gend = ginfo['end']
				forward = "FALSE"
				if ginfo['direction'] == '+': forward = "TRUE"

				cog_color = '"#dbdbdb"'
				if cog in cog_to_color:
					cog_color = '"' + cog_to_color[cog] + '"'

				gene_string = '\t'.join([str(x) for x in [bgc, lt, gstart, gend, forward, cog, cog_color]])
				printlist.append(gene_string)
				if cog != 'singleton':
					bgc_cog_presence[bgc][cog] = cog
					cog_counts[cog] += bgc_weights[bgc]
					cog_directions[cog] = ginfo['direction']
					cog_lengths[cog].append(gend - gstart)
			if i == 0:
				ref_cog_directions = cog_directions
				gggenes_track_handle.write('\n'.join(printlist) + '\n')
			else:
				flip_support = 0
				keep_support = 0
				for c in ref_cog_directions:
					if not c in cog_directions: continue
					cog_weight = statistics.mean(cog_lengths[c])
					if cog_directions[c] == ref_cog_directions[c]:
						keep_support += cog_weight
					else:
						flip_support += cog_weight

				# flip the genbank visual if necessary, first BGC processed is used as reference guide
				if flip_support > keep_support:
					flip_printlist = []
					for gene_string in printlist:
						gene_info = gene_string.split('\t')
						new_forward = 'TRUE'
						if gene_info[4] == 'TRUE': new_forward = 'FALSE'
						new_gstart = int(last_gene_end) - int(gene_info[3])
						new_gend = int(last_gene_end) - int(gene_info[2])
						new_gene_string = '\t'.join([str(x) for x in [gene_info[0], gene_info[1], new_gstart, new_gend, new_forward, gene_info[-2], gene_info[-1]]])
						flip_printlist.append(new_gene_string)
					gggenes_track_handle.write('\n'.join(flip_printlist) + '\n')
				else:
					gggenes_track_handle.write('\n'.join(printlist) + '\n')
		gggenes_track_handle.close()

		for bgc in bgc_cog_presence:
			for cog in cog_counts:
				heatmap_track_handle.write('\t'.join([bgc, cog, bgc_cog_presence[bgc][cog], str(cog_counts[cog])]) + '\n')
		heatmap_track_handle.close()

	rscript_plot_cmd = ["Rscript", RSCRIPT_FOR_BGSEE, phylogeny_file, gggenes_track_file, heatmap_track_file, result_pdf_file]
	logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
	try:
		subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=sys.stderr, stderr=sys.stderr,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
	logObject.info('Plotting completed!')

def createItolBGCSeeTrack(track_file, bgc_genes, gene_to_cog, cog_to_color, comp_gene_info, dataset_label, logObject):
	track_handle = open(track_file, 'w')

	logObject.info("Writing iTol track file to: %s" % track_file)
	logObject.info("Track will have label: %s" % dataset_label)

	# write header for iTol track file
	track_handle.write('DATASET_DOMAINS\n')
	track_handle.write('SEPARATOR TAB\n')
	track_handle.write('DATASET_LABEL\t%s\n' % dataset_label)
	track_handle.write('COLOR\t#000000\n')
	track_handle.write('BORDER_WIDTH\t1\n')
	track_handle.write('BORDER_COLOR\t#000000\n')
	track_handle.write('SHOW_DOMAIN_LABELS\t0\n')
	track_handle.write('DATA\n')

	# write the rest of the iTol track file for illustrating genes across BGC instances
	ref_cog_directions = {}
	bgc_gene_counts = defaultdict(int)
	for bgc in bgc_genes:
		bgc_gene_counts[bgc] = len(bgc_genes[bgc])

	for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
		bgc = item[0]
		curr_bgc_genes = bgc_genes[bgc]
		last_gene_end = max([comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
		printlist = [bgc, str(last_gene_end)]
		cog_directions = {}
		cog_lengths = defaultdict(list)
		for lt in curr_bgc_genes:
			ginfo = comp_gene_info[lt]
			cog = 'singleton'
			if lt in gene_to_cog:
				cog = gene_to_cog[lt]
			shape = 'None'
			if ginfo['direction'] == '+':
				shape = 'TR'
			elif ginfo['direction'] == '-':
				shape = 'TL'
			gstart = ginfo['start']
			gend = ginfo['end']
			cog_color = "#dbdbdb"
			if cog in cog_to_color:
				cog_color = cog_to_color[cog]
			gene_string = '|'.join([str(x) for x in [shape, gstart, gend, cog_color, cog]])
			printlist.append(gene_string)
			if cog != 'singleton':
				cog_directions[cog] = ginfo['direction']
				cog_lengths[cog].append(gend - gstart)
		if i == 0:
			ref_cog_directions = cog_directions
			track_handle.write('\t'.join(printlist) + '\n')
		else:
			flip_support = 0
			keep_support = 0
			for c in ref_cog_directions:
				if not c in cog_directions: continue
				cog_weight = statistics.mean(cog_lengths[c])
				if cog_directions[c] == ref_cog_directions[c]:
					keep_support += cog_weight
				else:
					flip_support += cog_weight

			# flip the genbank visual if necessary, first BGC processed is used as reference guide
			if flip_support > keep_support:
				flip_printlist = printlist[:2]
				for gene_string in printlist[2:]:
					gene_info = gene_string.split('|')
					new_shape = None
					if gene_info[0] == 'TR':
						new_shape = 'TL'
					elif gene_info[0] == 'TL':
						new_shape = 'TR'
					new_gstart = int(last_gene_end) - int(gene_info[2])
					new_gend = int(last_gene_end) - int(gene_info[1])
					new_gene_info = '|'.join([new_shape, str(new_gstart), str(new_gend)] + gene_info[-2:])
					flip_printlist.append(new_gene_info)
				track_handle.write('\t'.join(flip_printlist) + '\n')
			else:
				track_handle.write('\t'.join(printlist) + '\n')
	track_handle.close()


def determineCogOrderIndex(bgc_genes, gene_to_cog, comp_gene_info):
	ref_cog_directions = {}
	bgc_gene_counts = defaultdict(int)
	for bgc in bgc_genes:
		bgc_gene_counts[bgc] = len(bgc_genes[bgc])

	cog_order_scores = defaultdict(int)
	for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
		bgc = item[0]
		curr_bgc_genes = bgc_genes[bgc]
		cog_directions = {}
		cog_lengths = defaultdict(list)
		cog_starts = {}
		for g in curr_bgc_genes:
			ginfo = comp_gene_info[g]
			gstart = ginfo['start']
			gend = ginfo['end']
			if g in gene_to_cog:
				cog = gene_to_cog[g]
				cog_directions[cog] = ginfo['direction']
				cog_lengths[cog].append(gend - gstart)
				cog_starts[cog] = ginfo['start']
		reverse_flag = False
		if i == 0:
			ref_cog_directions = cog_directions
		else:
			flip_support = 0
			keep_support = 0
			for c in ref_cog_directions:
				if not c in cog_directions: continue
				cog_weight = statistics.mean(cog_lengths[c])
				if cog_directions[c] == ref_cog_directions[c]:
					keep_support += cog_weight
				else:
					flip_support += cog_weight

			# reverse ordering
			if flip_support > keep_support:
				reverse_flag = True
		for c in sorted(cog_starts.items(), key=itemgetter(1), reverse=reverse_flag):
			cog_order_scores[c[0]] += c[1]

	return cog_order_scores


def constructBGCPhylogeny(codon_alignments_dir, output_prefix, logObject):
	bgc_sccs = defaultdict(lambda: "")
	fasta_data = []
	fasta_data_tr = []

	for f in os.listdir(codon_alignments_dir):
		cog_align_msa = codon_alignments_dir + f
		# concatenate gene alignments
		with open(cog_align_msa) as opm:
			for rec in SeqIO.parse(opm, 'fasta'):
				bgc_sccs['>' + rec.id] += str(rec.seq).upper()

	for b in bgc_sccs:
		fasta_data.append([b] + list(bgc_sccs[b]))

	for i, ls in enumerate(zip(*fasta_data)):
		if i == 0:
			fasta_data_tr.append(ls)
		else:
			n_count = len([x for x in ls if x == '-'])
			if (float(n_count) / len(ls)) < 0.1:
				fasta_data_tr.append(list(ls))

	scc_fna = output_prefix + '.fna'
	scc_tre = output_prefix + '.tre'
	scc_handle = open(scc_fna, 'w')

	for rec in zip(*fasta_data_tr):
		scc_handle.write(rec[0] + '\n' + ''.join(rec[1:]) + '\n')
	scc_handle.close()

	# use FastTree2 to construct phylogeny
	fasttree_cmd = ['fasttree', '-nt', scc_fna, '>', scc_tre]
	logObject.info('Running FastTree2 with the following command: %s' % ' '.join(fasttree_cmd))
	try:
		subprocess.call(' '.join(fasttree_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(fasttree_cmd))
	except:
		logObject.warning('Had an issue running: %s' % ' '.join(fasttree_cmd))

	return scc_tre


def constructCodonAlignments(bgc_sample, cog_genes, comp_gene_info, outdir, cores, logObject, only_scc=False):
	nucl_seq_dir = os.path.abspath(outdir + 'Nucleotide_Sequences') + '/'
	prot_seq_dir = os.path.abspath(outdir + 'Protein_Sequences') + '/'
	prot_alg_dir = os.path.abspath(outdir + 'Protein_Alignments') + '/'
	codo_alg_dir = os.path.abspath(outdir + 'Codon_Alignments') + '/'
	if not os.path.isdir(nucl_seq_dir): os.system('mkdir %s' % nucl_seq_dir)
	if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
	if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
	if not os.path.isdir(codo_alg_dir): os.system('mkdir %s' % codo_alg_dir)

	all_samples = set(bgc_sample.values())
	try:
		inputs = []
		for cog in cog_genes:
			if len(cog_genes[cog]) < 2: continue
			sample_counts = defaultdict(int)
			gene_sequences = {}
			for gene in cog_genes[cog]:
				gene_info = comp_gene_info[gene]
				bgc_id = gene_info['bgc_name']
				sample_id = bgc_sample[bgc_id]
				nucl_seq = gene_info['nucl_seq']
				prot_seq = gene_info['prot_seq']
				sample_counts[sample_id] += 1
				gid = sample_id + '|' + gene
				if only_scc:
					gid = sample_id
				gene_sequences[gid] = tuple([nucl_seq, prot_seq])
			samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
			# check that cog is single-copy-core
			if only_scc and len(samples_with_single_copy.symmetric_difference(all_samples)) > 0:
				continue
			elif only_scc:
				logObject.info('Homolog group %s detected as SCC across samples (not individual BGCs).' % cog)
			inputs.append([cog, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, logObject])

		p = multiprocessing.Pool(cores)
		p.map(create_codon_msas, inputs)

		return codo_alg_dir
	except:
		logObject.error("Issues with create protein/codon alignments of SCC homologs for BGC.")
		raise RuntimeError("Issues with create protein/codon alignments of SCC homologs for BGC.")


def create_codon_msas(inputs):
	cog, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, logObject = inputs

	cog_nucl_fasta = nucl_seq_dir + '/' + cog + '.fna'
	cog_prot_fasta = prot_seq_dir + '/' + cog + '.faa'
	cog_prot_msa = prot_alg_dir + '/' + cog + '.msa.faa'
	cog_codo_msa = codo_alg_dir + '/' + cog + '.msa.fna'

	cog_nucl_handle = open(cog_nucl_fasta, 'w')
	cog_prot_handle = open(cog_prot_fasta, 'w')
	for s in gene_sequences:
		cog_nucl_handle.write('>' + s + '\n' + str(gene_sequences[s][0]) + '\n')
		cog_prot_handle.write('>' + s + '\n' + str(gene_sequences[s][1]) + '\n')
	cog_nucl_handle.close()
	cog_prot_handle.close()

	mafft_cmd = ['mafft', '--maxiterate', '1000', '--localpair', cog_prot_fasta, '>', cog_prot_msa]
	pal2nal_cmd = ['pal2nal.pl', cog_prot_msa, cog_nucl_fasta, '-output', 'fasta', '>', cog_codo_msa]

	logObject.info('Running mafft with the following command: %s' % ' '.join(mafft_cmd))
	try:
		subprocess.call(' '.join(mafft_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(mafft_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(mafft_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(mafft_cmd))

	logObject.info('Running PAL2NAL with the following command: %s' % ' '.join(pal2nal_cmd))
	try:
		subprocess.call(' '.join(pal2nal_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(pal2nal_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(pal2nal_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(pal2nal_cmd))

	logObject.info('Achieved codon alignment for homolog group %s' % cog)


def mapComprehensiveBGCsList(all_bgcs_file, logObject):
	logObject.info("Reading in file listing comprehensive list of BGCs.")
	comprehensive_bgcs = []
	try:
		with open(all_bgcs_file) as oabf:
			for line in oabf:
				line = line.strip()
				sample, bgc_genbank, bgc_proteome = line.split('\t')
				assert (is_genbank(bgc_genbank))
				assert (is_fasta(bgc_proteome))
				comprehensive_bgcs.append(tuple([sample, bgc_genbank, bgc_proteome]))
		return comprehensive_bgcs
	except:
		logObject.error(
			"Had difficulties in processing file listing comprehensive set of BGCs and paths to genbanks and extracted proteins.")
		raise RuntimeError(
			"Had difficulties in processing file listing comprehensive set of BGCs and paths to genbanks and extracted proteins.")


def runHMMScanAndAssignBGCsToGCF(comprehensive_bgcs, concat_hmm_profiles, scc_homologs, orthofinder_matrix_file, outdir,
								 cores, logObject):
	search_res_dir = os.path.abspath(outdir + 'HMMScan_Results') + '/'
	if not os.path.isdir(search_res_dir): os.system('mkdir %s' % search_res_dir)

	tot_bgc_proteins = defaultdict(int)
	hmmscan_cmds = []
	for i, cb_tuple in enumerate(comprehensive_bgcs):
		sample, bgc_genbank, bgc_proteome = cb_tuple
		with open(bgc_proteome) as obp:
			for rec in SeqIO.parse(obp, 'fasta'):
				tot_bgc_proteins[bgc_genbank] += 1
		result_file = search_res_dir + str(i) + '.txt'
		hmmscan_cmd = ['hmmscan', '--max', '--cpu', '1', '--tblout', result_file, concat_hmm_profiles, bgc_proteome,
					   logObject]
		hmmscan_cmds.append(hmmscan_cmd)

	p = multiprocessing.Pool(cores)
	p.map(multiProcess, hmmscan_cmds)
	p.close()

	protein_hits = defaultdict(list)
	sample_cogs = defaultdict(set)
	bgc_cogs = defaultdict(set)
	sample_bgcs = defaultdict(set)
	sample_cog_proteins = defaultdict(lambda: defaultdict(set))
	for i, cb_tuple in enumerate(comprehensive_bgcs):
		sample, bgc_genbank, bgc_proteome = cb_tuple
		sample_bgcs[sample].add(bgc_genbank)
		result_file = search_res_dir + str(i) + '.txt'
		assert (os.path.isfile(result_file))
		with open(result_file) as orf:
			for line in orf:
				if line.startswith("#"): continue
				line = line.strip()
				ls = line.split()
				cog = ls[0]
				samp, bgc, protein_id = ls[2].split('|')
				eval = float(ls[4])
				if eval <= 1e-5:
					protein_hits[protein_id].append([cog, eval, sample, bgc_genbank])

	for p in protein_hits:
		for i, hits in enumerate(sorted(protein_hits[p], key=itemgetter(1))):
			if i == 0:
				sample_cogs[hits[2]].add(hits[0])
				bgc_cogs[hits[3]].add(hits[0])
				sample_cog_proteins[hits[2]][hits[0]].add(p)

	expanded_orthofinder_matrix_file = outdir + 'Orthogroups.expanded.csv'
	expanded_gcf_list_file = outdir + 'GCF_Expanded.txt'

	expanded_orthofinder_matrix_handle = open(expanded_orthofinder_matrix_file, 'w')
	expanded_gcf_list_handle = open(expanded_gcf_list_file, 'w')

	valid_bgcs = set([])
	for sample in sample_cogs:
		if len(scc_homologs.difference(sample_cogs[sample])) == 0:
			for bgc_gbk in sample_bgcs[sample]:
				if float(len(bgc_cogs[bgc_gbk])) / tot_bgc_proteins[bgc_gbk] >= 0.5 and len(
						bgc_cogs[bgc_gbk]) >= 3 and len(scc_homologs.intersection(bgc_cogs[bgc_gbk])) >= 1:
					valid_bgcs.add(bgc_gbk)

	bgc_cogs = defaultdict(set)
	sample_cog_proteins = defaultdict(lambda: defaultdict(set))
	for p in protein_hits:
		for i, hits in enumerate(sorted(protein_hits[p], key=itemgetter(1))):
			if i != 0 or not hits[3] in valid_bgcs: continue
			bgc_cogs[hits[3]].add(hits[0])
			sample_cog_proteins[hits[2]][hits[0]].add(p)

	all_samples = set([])
	for sample in sample_cogs:
		scc_check = True
		for cog in scc_homologs:
			if len(sample_cog_proteins[sample][cog]) != 1: scc_check = False
		if not scc_check: continue
		for bgc_gbk in sample_bgcs[sample]:
			if float(len(bgc_cogs[bgc_gbk])) / tot_bgc_proteins[bgc_gbk] >= 0.5 and len(bgc_cogs[bgc_gbk]) >= 3 and len(
					scc_homologs.intersection(bgc_cogs[bgc_gbk])) >= 1:
				expanded_gcf_list_handle.write('\t'.join([sample, bgc_gbk]) + '\n')
				all_samples.add(sample)

	original_samples = []
	all_cogs = set([])
	with open(orthofinder_matrix_file) as omf:
		for i, line in enumerate(omf):
			line = line.strip('\n')
			ls = line.split('\t')
			if i == 0:
				original_samples = ls[1:]
				all_samples = all_samples.union(set(original_samples))
			else:
				cog = ls[0]
				all_cogs.add(cog)
				for j, prot in enumerate(ls[1:]):
					sample_cog_proteins[original_samples[j]][cog] = sample_cog_proteins[original_samples[j]][cog].union(
						set(prot.split(', ')))

	header = [''] + [s for s in sorted(all_samples)]
	expanded_orthofinder_matrix_handle.write('\t'.join(header) + '\n')
	for c in sorted(all_cogs):
		printlist = [c]
		for s in sorted(all_samples):
			printlist.append(', '.join(sample_cog_proteins[s][c]))
		expanded_orthofinder_matrix_handle.write('\t'.join(printlist) + '\n')

	expanded_gcf_list_handle.close()
	expanded_orthofinder_matrix_handle.close()

def constructHMMProfiles(bgc_sample, cog_genes, comp_gene_info, outdir, cores, logObject):
	prot_seq_dir = os.path.abspath(outdir + 'Protein_Sequences') + '/'
	prot_alg_dir = os.path.abspath(outdir + 'Protein_Alignments') + '/'
	prot_hmm_dir = os.path.abspath(outdir + 'Profile_HMMs') + '/'
	if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
	if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
	if not os.path.isdir(prot_hmm_dir): os.system('mkdir %s' % prot_hmm_dir)

	all_samples = set(bgc_sample.values())
	scc_homologs = set([])
	try:
		inputs = []
		for cog in cog_genes:
			sample_counts = defaultdict(int)
			sample_sequences = {}
			for gene in cog_genes[cog]:
				gene_info = comp_gene_info[gene]
				bgc_id = gene_info['bgc_name']
				sample_id = bgc_sample[bgc_id]
				prot_seq = gene_info['prot_seq']
				sample_counts[sample_id] += 1
				sample_sequences[sample_id] = prot_seq
			samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
			# check that cog is single-copy-core
			if len(samples_with_single_copy.symmetric_difference(all_samples)) == 0: scc_homologs.add(cog)
			logObject.info('Homolog group %s detected as SCC across samples (not individual BGCs).' % cog)
			inputs.append([cog, sample_sequences, prot_seq_dir, prot_alg_dir, prot_hmm_dir, logObject])

		p = multiprocessing.Pool(cores)
		p.map(create_hmm_profiles, inputs)

		logObject.info(
			"Successfully created profile HMMs for each homolog group. Now beginning concatenation into single file.")
		concatenated_profile_HMM = outdir + 'All_GCF_Homologs.hmm'
		os.system('rm -f %s' % concatenated_profile_HMM)
		for f in os.listdir(prot_hmm_dir):
			os.system('cat %s >> %s' % (prot_hmm_dir + f, concatenated_profile_HMM))

		hmmpress_cmd = ['hmmpress', concatenated_profile_HMM]
		logObject.info(
			'Running hmmpress on concatenated profiles with the following command: %s' % ' '.join(hmmpress_cmd))
		try:
			subprocess.call(' '.join(hmmpress_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran: %s' % ' '.join(hmmpress_cmd))
		except:
			logObject.error('Had an issue running: %s' % ' '.join(hmmpress_cmd))
			raise RuntimeError('Had an issue running: %s' % ' '.join(hmmpress_cmd))

		return scc_homologs, concatenated_profile_HMM
	except:
		logObject.error("Issues with running hmmpress on profile HMMs.")
		raise RuntimeError("Issues with running hmmpress on profile HMMs.")

def create_hmm_profiles(inputs):
	cog, sample_sequences, prot_seq_dir, prot_alg_dir, prot_hmm_dir, logObject = inputs

	cog_prot_fasta = prot_seq_dir + '/' + cog + '.faa'
	cog_prot_msa = prot_alg_dir + '/' + cog + '.msa.faa'
	cog_prot_hmm = prot_hmm_dir + '/' + cog + '.hmm'

	cog_prot_handle = open(cog_prot_fasta, 'w')
	for s in sample_sequences:
		cog_prot_handle.write('>' + s + '\n' + str(sample_sequences[s]) + '\n')
	cog_prot_handle.close()

	mafft_cmd = ['mafft', '--maxiterate', '1000', '--localpair', cog_prot_fasta, '>', cog_prot_msa]
	logObject.info('Running mafft with the following command: %s' % ' '.join(mafft_cmd))
	try:
		subprocess.call(' '.join(mafft_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(mafft_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(mafft_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(mafft_cmd))

	hmmbuild_cmd = ['hmmbuild', '--amino', '-n', cog, cog_prot_hmm, cog_prot_msa]
	logObject.info('Running hmmbuild (from HMMER3) with the following command: %s' % ' '.join(hmmbuild_cmd))
	try:
		subprocess.call(' '.join(hmmbuild_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(hmmbuild_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(hmmbuild_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(hmmbuild_cmd))

	logObject.info('Constructed profile HMM for homolog group %s' % cog)


def modifyPhylogenyForSamplesWithMultipleBGCs(phylogeny, sample_bgcs, output_phylogeny, logObject):
	try:
		number_of_added_leaves = 0
		t = Tree(phylogeny)
		for node in t.traverse('postorder'):
			if node.name in sample_bgcs and len(sample_bgcs[node.name]) > 1:
				og_node_name = node.name
				node.name = node.name + '_INNERNODE'
				for bgc_id in sample_bgcs[og_node_name]:
					# if bgc_id == node.name: continue
					node.add_child(name=bgc_id)
					child_node = t.search_nodes(name=bgc_id)[0]
					child_node.dist = 0
					if bgc_id != og_node_name: number_of_added_leaves += 1
		t.write(format=1, outfile=output_phylogeny)
		logObject.info(
			"New phylogeny with an additional %d leafs to reflect samples with multiple BGCs can be found at: %s." % (
				number_of_added_leaves, output_phylogeny))
	except:
		logObject.error(
			"Had difficulties properly editing phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")
		raise RuntimeError(
			"Had difficulties properly editing phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")


def readInBGCGenbanksComprehensive(bgc_specs_file, logObject, comprehensive_parsing=True):
	bgc_sample = {}
	bgc_product = {}
	bgc_genes = {}
	bgc_core_counts = {}
	all_genes = set([])
	with open(bgc_specs_file) as obsf:
		for i, line in enumerate(obsf):
			line = line.strip()
			try:
				assert (len(line.split('\t')) == 2)
			except:
				logObject.error(
					"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
							i + 1))
				raise RuntimeError(
					"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
							i + 1))
			sample, gbk = line.split('\t')
			try:
				assert (is_genbank(gbk))
				bgc_genes_full, bgc_info_full = parseGenbanks(gbk, gbk, comprehensive_parsing=comprehensive_parsing)
				bgc_product[gbk] = [x['product'] for x in bgc_info_full]
				bgc_core_counts[gbk] = bgc_info_full[0]['count_core_gene_groups']
				bgc_genes[gbk] = set(bgc_genes_full.keys())
				all_genes = all_genes.union(bgc_genes[gbk])
				bgc_sample[gbk] = sample
				logObject.info("Incorporating genbank %s for sample %s into analysis." % (gbk, sample))
			except:
				logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis.")
				raise RuntimeWarning("Unable to validate %s as Genbank. Skipping ...")

	return [bgc_sample, bgc_product, bgc_genes, bgc_core_counts, all_genes]

def getSpeciesRelationshipsFromPhylogeny(species_phylogeny, samples_in_gcf):
	samples_in_phylogeny = set([])
	t = Tree(species_phylogeny)
	for leaf in t:
		samples_in_phylogeny.add(str(leaf).strip('\n').lstrip('-'))

	pairwise_distances = defaultdict(lambda: defaultdict(float))
	for s1 in samples_in_gcf.intersection(samples_in_phylogeny):
		for s2 in samples_in_gcf.intersection(samples_in_phylogeny):
			try:
				s1_s2_dist = t.get_distance(s1, s2)
				pairwise_distances[s1][s2] = s1_s2_dist
			except:
				pass
	return ([pairwise_distances, samples_in_phylogeny])


def readInPopulationsSpecification(pop_specs_file, logObject):
	populations = set([])
	sample_population = defaultdict(lambda: "NA")
	with open(pop_specs_file) as opsf:
		for line in opsf:
			line = line.strip()
			sample, population = line.split('\t')
			sample_population[sample] = population
			populations.add(population)
	logObject.info("Successfully parsed population specifications file. There are %d populations" % len(population))
	return sample_population


"""
def writeLociFileForPegas(codon_msa, outfile, logObject, sample_population=None):
	outfile_handle = open(outfile, 'w')
	number_of_loci = 0
	with open(codon_msa) as ocmf:
		for i, rec in enumerate(SeqIO.parse(ocmf, 'fasta')):
			if i == 0:
				number_of_loci = len(str(rec.seq))
				header = ['Sample'] + ['Pos' + str(p) for p in range(1, len(str(rec.seq))+1)]
				if sample_population:
					header += ['Pop']
				outfile_handle.write('\t'.join(header) + '\n')
			sample_data = [rec.id] + [p.replace('-', 'NA').upper() + '/' + p.replace('-', 'NA').upper() for p in str(rec.seq)]
			if sample_population:
				try:
					sample_data += [sample_population[rec.id]]
				except:
					logObject.warning("Didn't find population for sample %s, filling in population as NA." % rec.id)
					sample_data += ['NA']
			outfile_handle.write('\t'.join(sample_data) + '\n')

	outfile_handle.close()
	return [number_of_loci, outfile]
"""


def parseCodonAlignmentStats(cog, codon_alignment_fasta, comp_gene_info, cog_genes, cog_order_index,
							 cog_median_copy_numbers, plots_dir, popgen_dir, bgc_sample, final_output_handle, logObject,
							 sample_population=None):
	domain_plot_file = plots_dir + cog + '_domain.txt'
	position_plot_file = plots_dir + cog + '_position.txt'
	popgen_plot_file = plots_dir + cog + '_popgen.txt'
	plot_pdf_file = plots_dir + cog + '.pdf'

	domain_plot_handle = open(domain_plot_file, 'w')
	position_plot_handle = open(position_plot_file, 'w')
	popgen_plot_handle = open(popgen_plot_file, 'w')

	seqs = []
	bgc_codons = defaultdict(list)
	num_codons = None
	samples = set([])
	gene_lengths = []
	gene_locs = defaultdict(dict)
	core_counts = defaultdict(int)
	products = set([])
	with open(codon_alignment_fasta) as ocaf:
		for rec in SeqIO.parse(ocaf, 'fasta'):
			sample_id, gene_id = rec.id.split('|')
			if comp_gene_info[gene_id]['core_overlap']:
				core_counts['core'] += 1
			else:
				core_counts['auxiliary'] += 1
			products.add(comp_gene_info[gene_id]['product'])
			real_pos = 1
			seqs.append(list(str(rec.seq)))
			codons = [str(rec.seq)[i:i + 3] for i in range(0, len(str(rec.seq)), 3)]
			num_codons = len(codons)
			bgc_codons[rec.id] = codons
			samples.add(sample_id)
			for msa_pos, bp in enumerate(str(rec.seq)):
				if bp != '-':
					gene_locs[gene_id][real_pos] = msa_pos + 1
					real_pos += 1
			gene_lengths.append(len(str(rec.seq).replace('-', '')))

	is_core = False
	if float(core_counts['core']) / sum(core_counts.values()) >= 0.8: is_core = True

	median_gene_length = statistics.median(gene_lengths)

	variable_sites = set([])
	conserved_sites = set([])
	position_plot_handle.write('\t'.join(['pos', 'num_seqs', 'num_alleles', 'num_gaps', 'maj_allele_freq']) + '\n')
	for i, ls in enumerate(zip(*seqs)):
		al_counts = defaultdict(int)
		for al in ls:
			if al != '-': al_counts[al] += 1
		maj_allele_count = max(al_counts.values())
		tot_count = sum(al_counts.values())
		num_seqs = len(ls)
		num_alleles = len(al_counts.keys())
		num_gaps = num_seqs - tot_count
		maj_allele_freq = float(maj_allele_count) / tot_count
		position_plot_handle.write(
			'\t'.join([str(x) for x in [i + 1, num_seqs, num_alleles, num_gaps, maj_allele_freq]]) + '\n')
		if maj_allele_freq <= 0.90:
			variable_sites.add(i)
		else:
			conserved_sites.add(i)
	position_plot_handle.close()

	differential_domains = set([])
	domain_positions_msa = defaultdict(set)
	domain_min_position_msa = defaultdict(lambda: 1e8)
	all_domains = set([])
	for gene in cog_genes[cog]:
		gene_start = comp_gene_info[gene]['start']
		gene_end = comp_gene_info[gene]['end']
		for domain in comp_gene_info[gene]['gene_domains']:
			domain_start = max(domain['start'], gene_start)
			domain_end = min(domain['end'], gene_end)
			domain_name = domain['aSDomain'] + '_|_' + domain['description']
			relative_start = domain_start - gene_start
			assert (len(gene_locs[gene]) + 3 >= (domain_end - gene_start))
			relative_end = min([len(gene_locs[gene]), domain_end - gene_start])
			domain_range = range(relative_start, relative_end)
			for pos in domain_range:
				msa_pos = gene_locs[gene][pos + 1]
				domain_positions_msa[domain_name].add(msa_pos)
				if domain_min_position_msa[domain_name] > msa_pos:
					domain_min_position_msa[domain_name] = msa_pos
			all_domains.add(domain['type'] + '_|_' + domain['aSDomain'] + '_|_' + domain['description'])

	domain_plot_handle.write('\t'.join(['domain', 'domain_index', 'min_pos', 'max_pos']) + '\n')
	for i, dom in enumerate(sorted(domain_min_position_msa.items(), key=itemgetter(1))):
		tmp = []
		old_pos = None
		for j, pos in enumerate(sorted(domain_positions_msa[dom[0]])):
			if j == 0:
				old_pos = pos - 1
			if pos - 1 != old_pos:
				if len(tmp) > 0:
					min_pos = min(tmp)
					max_pos = max(tmp)
					domain_plot_handle.write('\t'.join([str(x) for x in [dom[0], i, min_pos, max_pos]]) + '\n')
				tmp = []
			tmp.append(pos)
			old_pos = pos
		if len(tmp) > 0:
			min_pos = min(tmp)
			max_pos = max(tmp)
			domain_plot_handle.write('\t'.join([str(x) for x in [dom[0], i, min_pos, max_pos]]) + '\n')
	domain_plot_handle.close()

	popgen_plot_handle.write('\t'.join(['pos', 'type']) + '\n')
	total_core_codons = 0
	total_variable_codons = 0
	nonsynonymous_sites = 0
	synonymous_sites = 0
	for cod_index in range(0, num_codons):
		first_bp = (cod_index + 1) * 3
		aa_count = defaultdict(int)
		aa_codons = defaultdict(set)
		cod_count = defaultdict(int)
		core = True
		for bgc in bgc_codons:
			cod = bgc_codons[bgc][cod_index]
			if '-' in cod or 'N' in cod:
				core = False
			else:
				cod_obj = Seq(cod)
				aa_count[str(cod_obj.translate())] += 1
				cod_count[cod] += 1
				aa_codons[str(cod_obj.translate())].add(cod)

		residues = len([r for r in aa_count if aa_count[r] >= 2])
		residues_with_multicodons = 0
		for r in aa_codons:
			supported_cods = 0
			for cod in aa_codons[r]:
				if cod_count[cod] >= 2: supported_cods += 1
			if supported_cods >= 2: residues_with_multicodons += 1

		maj_allele_count = max(cod_count.values())
		tot_valid_codons = sum(cod_count.values())
		maj_allele_freq = float(maj_allele_count) / tot_valid_codons

		nonsyn_flag = False;
		syn_flag = False
		if maj_allele_freq <= 0.9:
			total_variable_codons += 1
			if len(cod_count.keys()) > 1:
				if residues >= 2: nonsynonymous_sites += 1; nonsyn_flag = True
				if residues_with_multicodons >= 1: synonymous_sites += 1; syn_flag = True

		if core and not nonsyn_flag and not syn_flag: total_core_codons += 1
		if (nonsyn_flag and not syn_flag) or (not nonsyn_flag and syn_flag):
			type = 'S'
			if nonsyn_flag: type = 'NS'
			popgen_plot_handle.write('\t'.join([str(x) for x in [first_bp, type]]) + '\n')
	popgen_plot_handle.close()
	dn_ds = "NA"
	if synonymous_sites > 0: dn_ds = float(nonsynonymous_sites) / synonymous_sites

	rscript_plot_cmd = ["Rscript", RSCRIPT_FOR_PLOTTING, domain_plot_file, position_plot_file, popgen_plot_file,
						plot_pdf_file]
	logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
	try:
		subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
	logObject.info('Plotting completed!')

	tajima_results = popgen_dir + cog + '.tajima.txt'
	rscript_tajimaD_cmd = ["Rscript", RSCRIPT_FOR_TAJIMA, codon_alignment_fasta, tajima_results]
	logObject.info(
		'Running R pegas for calculating Tajima\'s D from codon alignment with the following command: %s' % ' '.join(
			rscript_tajimaD_cmd))
	try:
		subprocess.call(' '.join(rscript_tajimaD_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(rscript_tajimaD_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(rscript_tajimaD_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_tajimaD_cmd))
	logObject.info('Tajima\'s D has been calculated!')

	tajimas_d = "NA"
	if os.path.isfile(tajima_results):
		with open(tajima_results) as otrf:
			for i, line in enumerate(otrf):
				if i == 4:
					try:
						tajimas_d = float(line.strip())
					except:
						pass

	prop_samples_with_cog = len(samples) / float(len(set(bgc_sample.values())))

	cog_info = [cog, '; '.join(products), cog_order_index[cog], cog_median_copy_numbers[cog], median_gene_length,
				is_core, len(seqs), prop_samples_with_cog, tajimas_d, total_core_codons, total_variable_codons,
				nonsynonymous_sites, synonymous_sites, dn_ds, '; '.join(all_domains)]

	if sample_population:

		population_samples = defaultdict(set)
		for sp in sample_population.items():
			population_samples[sp[1]].add(sp[0])

		sample_seqs = {}
		with open(codon_alignment_fasta) as ocaf:
			for rec in SeqIO.parse(ocaf, 'fasta'):
				sample_seqs[rec.id.split('|')[0]] = str(rec.seq).upper()

		pairwise_differences = defaultdict(lambda: defaultdict(int))
		for i, s1 in enumerate(sample_seqs):
			for j, s2 in enumerate(sample_seqs):
				if i < j:
					for p, s1b in enumerate(sample_seqs[s1]):
						s2b = sample_seqs[s2][p]
						if s1b != s2b:
							pairwise_differences[s1][s2] += 1
							pairwise_differences[s2][s1] += 1

		pop_prop_with_cog = defaultdict(float)
		pops_with_cog = 0
		populations_order = []
		within_population_differences = []
		for pop in population_samples:
			pop_prop_with_cog[pop] = len(population_samples[pop].intersection(samples)) / len(population_samples[pop])
			if pop_prop_with_cog[pop] > 0: pops_with_cog += 1
			if len(population_samples[pop]) >= 2:
				within_pop = []
				for i, s1 in enumerate(population_samples[pop]):
					for j, s2 in enumerate(population_samples[pop]):
						if i < j:
							within_pop.append(pairwise_differences[s1][s2])
				within_population_differences.append(within_pop)
				populations_order.append(pop)

		anova_pval = "NA"
		if len(within_population_differences) >= 2:
			F, anova_pval = f_oneway(*within_population_differences)
		cog_population_info = [pops_with_cog,
							   '|'.join([str(x[0]) + '=' + str(x[1]) for x in pop_prop_with_cog.items()]),
							   anova_pval]
		cog_info += cog_population_info

	final_output_handle.write('\t'.join([str(x) for x in cog_info]) + '\n')


def parseGenbanks(gbk, bgc_name, comprehensive_parsing=True):
	"""
	:param gbk: AntiSMASH Genbank file
	:return:
	"""
	bgc_info = []
	domains = []
	core_positions = set([])
	full_sequence = ""
	with open(gbk) as ogbk:
		domain_feature_types = ['PFAM_domain', 'CDS_motif', 'aSDomain']
		for rec in SeqIO.parse(ogbk, 'genbank'):
			full_sequence = str(rec.seq)
			for feature in rec.features:
				if comprehensive_parsing and feature.type in domain_feature_types:
					start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					aSDomain = "NA"
					description = "NA"
					try:
						aSDomain = feature.qualifiers.get('aSDomain')[0]
					except:
						pass
					try:
						description = feature.qualifiers.get('description')[0]
					except:
						pass
					domains.append({'start': start, 'end': end, 'type': feature.type, 'aSDomain': aSDomain,
									'description': description})
				elif feature.type == 'protocluster':
					detection_rule = feature.qualifiers.get('detection_rule')[0]
					try:
						product = feature.qualifiers.get('product')[0]
					except:
						product = "NA"
					contig_edge = feature.qualifiers.get('contig_edge')[0]
					bgc_info.append(
						{'detection_rule': detection_rule, 'product': product, 'contig_edge': contig_edge,
						 'full_sequence': str(rec.seq)})
				elif feature.type == 'proto_core':
					core_start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					core_end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					core_positions = core_positions.union(set(range(core_start, core_end + 1)))

	if len(bgc_info) == 0:
		bgc_info = [{'detection_rule': 'NA', 'product': 'NA', 'contig_edge': 'NA', 'full_sequence': full_sequence}]

	sys.stderr.write('Processing %s\n' % gbk)
	genes = {}
	core_genes = set([])
	gene_order = {}
	with open(gbk) as ogbk:
		for rec in SeqIO.parse(ogbk, 'genbank'):
			for feature in rec.features:
				if feature.type == "CDS":
					lt = feature.qualifiers.get('locus_tag')[0]
					start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					direction = str(feature.location).split('(')[1].split(')')[0]

					try:
						product = feature.qualifiers.get('product')[0]
					except:
						product = "hypothetical protein"

					grange = set(range(start, end + 1))
					core_overlap = False
					if len(grange.intersection(core_positions)) > 0:
						core_overlap = True
						core_genes.add(lt)

					gene_order[lt] = start

					prot_seq, nucl_seq, nucl_seq_with_flanks, relative_start, relative_end, gene_domains = [None]*6
					if comprehensive_parsing:
						prot_seq = feature.qualifiers.get('translation')[0]
						gene_domains = []
						for d in domains:
							drange = set(range(d['start'], d['end'] + 1))
							if len(drange.intersection(grange)) > 0:
								gene_domains.append(d)

						flank_start = start - FLANK_SIZE
						flank_end = end + FLANK_SIZE
						if flank_start < 0: flank_start = 0
						if flank_end >= len(full_sequence): flank_end = None
						if end >= len(full_sequence): end = None
						if end:
							nucl_seq = full_sequence[start:end]
						else:
							nucl_seq = full_sequence[start:]
							end = len(full_sequence)
						if flank_end:
							nucl_seq_with_flanks = full_sequence[flank_start:flank_end]
						else:
							nucl_seq_with_flanks = full_sequence[flank_start:]
							flank_end = len(full_sequence)

						gene_length = end - start

						relative_start = start - flank_start
						relative_end = relative_start + gene_length

						if direction == '-':
							nucl_seq = str(Seq(nucl_seq).reverse_complement())
							nucl_seq_with_flanks = str(Seq(nucl_seq_with_flanks).reverse_complement())
							relative_end = len(nucl_seq_with_flanks) - relative_start
							relative_start = relative_end - gene_length

					genes[lt] = {'bgc_name': bgc_name, 'start': start, 'end': end, 'direction': direction,
								 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
								 'nucl_seq_with_flanks': nucl_seq_with_flanks, 'gene_domains': gene_domains,
								 'core_overlap': core_overlap, 'relative_start': relative_start,
								 'relative_end': relative_end}

	number_of_core_gene_groups = 0
	tmp = []
	for lt in sorted(gene_order.items(), key=itemgetter(1), reverse=True):
		if lt[0] in core_genes:
			tmp.append(lt[0])
		elif len(tmp) > 0:
			number_of_core_gene_groups += 1
			tmp = []
	if len(tmp) > 0:
		number_of_core_gene_groups += 1
		tmp = []

	for i, pc in enumerate(bgc_info):
		bgc_info[i]['count_core_gene_groups'] = number_of_core_gene_groups

	return ([genes, bgc_info])


def extractGeneWithFlanksAndCluster(bgc_genes, comp_gene_info, gene_to_cog, outdir, logObject):
	genes_with_flanks_fasta = outdir + 'GCF_Genes.fasta'
	gwff_handle = open(genes_with_flanks_fasta, 'w')
	for bgc in bgc_genes:
		for gene in bgc_genes[bgc]:
			if gene in gene_to_cog:
				gwff_handle.write('>' + gene + '|' + bgc + '|' + gene_to_cog[gene] + '\n' + comp_gene_info[gene][
					'nucl_seq_with_flanks'] + '\n')
	gwff_handle.close()

	cd_hit_nr_fasta_file = outdir + 'GCF_Genes_NR.fasta'
	cd_hit_clusters_fasta_file = outdir + 'GCF_Genes_Clusters.fasta'

	cd_hit_nr = ['cd-hit-est', '-i', genes_with_flanks_fasta, '-o', cd_hit_nr_fasta_file, '-G', '1', '-g',
				 '1', '-d', '0', '-n', '10', '-M', '2000', '-c', '1.0', '-aL', '0.0', '-aS', '1.0', '-T', '1']
	logObject.info('Running the following command: %s' % ' '.join(cd_hit_nr))
	try:
		subprocess.call(' '.join(cd_hit_nr), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(cd_hit_nr))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(cd_hit_nr))
		raise RuntimeError('Had an issue running: %s' % ' '.join(cd_hit_nr))
	logObject.info('Ran CD-HIT for collapsing redundancy.')

	cd_hit_cluster = ['cd-hit-est', '-i', genes_with_flanks_fasta, '-o', cd_hit_clusters_fasta_file, '-G', '1', '-g',
					  '1', '-d', '0', '-n', '10', '-M', '2000', '-c', '0.98', '-aL', '0.95', '-aS', '0.95', '-T', '1']
	logObject.info('Running the following command: %s' % ' '.join(cd_hit_cluster))
	try:
		subprocess.call(' '.join(cd_hit_cluster), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(cd_hit_cluster))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(cd_hit_cluster))
		raise RuntimeError('Had an issue running: %s' % ' '.join(cd_hit_cluster))
	logObject.info('Ran CD-HIT for clustering genes, with their flanks, into haplotype groups.')

	bowtie2_build = ['bowtie2-build', cd_hit_nr_fasta_file, cd_hit_nr_fasta_file.split('.fasta')[0]]
	logObject.info('Running the following command: %s' % ' '.join(bowtie2_build))
	try:
		subprocess.call(' '.join(bowtie2_build), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(bowtie2_build))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(bowtie2_build))
		raise RuntimeError('Had an issue running: %s' % ' '.join(bowtie2_build))
	logObject.info('Build Bowtie2 database/index for %s' % cd_hit_nr_fasta_file)

	cd_hit_clusters_cltr_file = cd_hit_clusters_fasta_file + '.clstr'
	assert (os.path.isfile(cd_hit_clusters_cltr_file))

	instance_to_haplotype = {}
	cluster = []
	with open(cd_hit_clusters_cltr_file) as off:
		for line in off:
			line = line.strip()
			ls = line.split()
			if line.startswith('>'):
				if len(cluster) > 0:
					for g in cluster:
						instance_to_haplotype[g] = rep
				cluster = []
				rep = None
			else:
				gene_id = ls[2][1:-3]
				cluster.append(gene_id)
				if line.endswith('*'): rep = gene_id
	if len(cluster) > 0:
		if len(cluster) > 0:
			for g in cluster:
				instance_to_haplotype[g] = rep

	return cd_hit_nr_fasta_file, instance_to_haplotype


def createSummaryMatricesForMetaNovelty(paired_end_sequencing_file, results_outdir, main_outdir, logObject):
	samples = set([])
	cog_allele_representatives = set([])
	sample_allele_reads = defaultdict(lambda: defaultdict(int))
	sample_allele_unique_reads = defaultdict(lambda: defaultdict(int))
	sample_allele_novelty_reads = defaultdict(lambda: defaultdict(int))
	sample_allele_unique_novelty_reads = defaultdict(lambda: defaultdict(int))
	with open(paired_end_sequencing_file) as ossf:
		for line in ossf:
			sample = line.strip().split('\t')[0]
			result_file = results_outdir + sample + '.txt'
			samples.add(sample)
			if not os.path.isfile(result_file): continue
			with open(result_file) as orf:
				for i, cog_al in enumerate(orf):
					if i == 0: continue
					cog_al = cog_al.strip()
					cog, allele_representative, reads, reads_with_novelty, reads_uniquely_mapping, reads_uniquely_mapping_with_novelty = cog_al.split(
						'\t')
					car = cog + '|' + allele_representative
					cog_allele_representatives.add(car)
					sample_allele_reads[sample][car] = int(reads)
					sample_allele_unique_reads[sample][car] = int(reads_uniquely_mapping)
					sample_allele_novelty_reads[sample][car] = int(reads_with_novelty)
					sample_allele_unique_novelty_reads[sample][car] = int(reads_uniquely_mapping_with_novelty)

	final_matrix_reads_file = main_outdir + 'Sample_by_OG_Allele_Read_Counts.matrix.txt'
	final_matrix_novelty_reads_file = main_outdir + 'Sample_by_OG_Allele_Novelty_Read_Counts.matrix.txt'
	final_matrix_unique_reads_file = main_outdir + 'Final_OG_Allele_Unique_Read_Counts.matrix.txt'
	final_matrix_unique_and_novelty_reads_file = main_outdir + 'Final_OG_Allele_Unique_and_Novelty_Read_Counts.matrix.txt'

	final_matrix_reads_handle = open(final_matrix_reads_file, 'w')
	final_matrix_novelty_reads_handle = open(final_matrix_novelty_reads_file, 'w')
	final_matrix_unique_reads_handle = open(final_matrix_unique_reads_file, 'w')
	final_matrix_unique_and_novelty_reads_handle = open(final_matrix_unique_and_novelty_reads_file, 'w')

	final_matrix_reads_handle.write('\t'.join(['Sample/OG_Allele'] + list(sorted(cog_allele_representatives))) + '\n')
	final_matrix_novelty_reads_handle.write(
		'\t'.join(['Sample/OG_Allele'] + list(sorted(cog_allele_representatives))) + '\n')
	final_matrix_unique_reads_handle.write(
		'\t'.join(['Sample/OG_Allele'] + list(sorted(cog_allele_representatives))) + '\n')
	final_matrix_unique_and_novelty_reads_handle.write(
		'\t'.join(['Sample/OG_Allele'] + list(sorted(cog_allele_representatives))) + '\n')

	for s in samples:
		printlist_all = [s]
		printlist_uni = [s]
		printlist_nov = [s]
		printlist_uni_nov = [s]
		for c in cog_allele_representatives:
			printlist_all.append(str(sample_allele_reads[s][c]))
			printlist_uni.append(str(sample_allele_unique_reads[s][c]))
			printlist_nov.append(str(sample_allele_novelty_reads[s][c]))
			printlist_uni_nov.append(str(sample_allele_unique_novelty_reads[s][c]))
		final_matrix_reads_handle.write('\t'.join(printlist_all) + '\n')
		final_matrix_novelty_reads_handle.write('\t'.join(printlist_nov) + '\n')
		final_matrix_unique_reads_handle.write('\t'.join(printlist_uni) + '\n')
		final_matrix_unique_and_novelty_reads_handle.write('\t'.join(printlist_uni_nov) + '\n')

	final_matrix_reads_handle.close()
	final_matrix_novelty_reads_handle.close()
	final_matrix_unique_reads_handle.close()


def generateNoveltyReport(results_outdir, codon_alignment_file, cog_prop_multicopy, comp_gene_info, outdir, logObject):
	novelty_report_file = outdir + 'Novelty_Report.txt'
	no_handle = open(novelty_report_file, 'w')

	gene_pos_to_msa_pos = defaultdict(lambda: defaultdict(dict))
	msa_pos_alleles = defaultdict(lambda: defaultdict(set))
	with open(codon_alignment_file) as ocaf:
		for line in ocaf:
			line = line.strip()
			cog, cod_alignment = line.split('\t')
			with open(cod_alignment) as oca:
				for rec in SeqIO.parse(oca, 'fasta'):
					sample_id, gene_id = rec.id.split('|')
					real_pos = 1
					for msa_pos, bp in enumerate(str(rec.seq)):
						if bp != '-':
							gene_pos_to_msa_pos[cog][gene_id][real_pos] = msa_pos + 1
							real_pos += 1
							msa_pos_alleles[cog][msa_pos + 1].add(bp.upper())

	for f in os.listdir(results_outdir):
		if not f.endswith('.snvs'): continue
		metsample = f.split('.snvs')[0]
		mges = set(['transp', 'integrase'])

		with open(results_outdir + f) as of:
			for i, line in enumerate(of):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				snv_count = ls[1]
				gsc, ref_pos, ref_al, alt_al = ls[0].split('_|_')
				gene, sample, cog = gsc.split('|')
				if cog_prop_multicopy[cog] >= 0.05: continue
				if any(word in comp_gene_info[gene]['product'].lower() for word in mges): continue
				# if any(word in ' '.join(comp_gene_info[gene]['gene_domains']).lower() for word in mges): continue

				# offset = comp_gene_info[gene]['relative_start']
				ref_pos = int(ref_pos)  # - offset + 2
				print(ref_pos)
				# print(offset)
				print(line)
				print(comp_gene_info[gene])
				if not int(ref_pos) in gene_pos_to_msa_pos[cog][gene]: continue
				msa_pos = gene_pos_to_msa_pos[cog][gene][int(ref_pos)]
				msa_pos_als = msa_pos_alleles[cog][msa_pos]
				print(msa_pos)
				print(msa_pos_als)
				print(msa_pos_alleles[cog][msa_pos - 1])
				print(msa_pos_alleles[cog][msa_pos + 1])
				print('\t'.join([metsample, cog, str(msa_pos), sample, gene, str(ref_pos), ref_al, alt_al, snv_count]))
				assert (ref_al in msa_pos_alleles[cog][msa_pos])
				if not alt_al in msa_pos_als:
					no_handle.write('\t'.join(
						[metsample, cog, str(msa_pos), alt_al, snv_count, sample, gene, str(ref_pos), ref_al]) + '\n')

	no_handle.close()


def runSNVMining(cog_genes, comp_gene_info, bowtie2_ref_fasta, paired_end_sequencing_file, instance_to_haplotype,
				 bowtie2_outdir, results_outdir, cores, logObject):
	process_args = []
	with open(paired_end_sequencing_file) as opesf:
		for line in opesf:
			line = line.strip()
			sample, frw_read, rev_read = line.split('\t')
			process_args.append([sample, bowtie2_outdir + sample + '.filtered.sorted.bam',
								 bowtie2_ref_fasta, instance_to_haplotype, results_outdir, cog_genes, comp_gene_info,
								 logObject])

	p = multiprocessing.Pool(cores)
	p.map(snv_miner, process_args)
	p.close()


def snv_miner(input_args):
	# try:
	sample, bam_alignment, ref_fasta, cog_gene_to_rep, res_dir, bgc_cog_genes, comp_gene_info, logObject = input_args
	cog_rep_genes = defaultdict(set)
	for g, r in cog_gene_to_rep.items():
		cog_rep_genes[r].add(g)

	if not os.path.isfile(bam_alignment): return
	snvs_file = res_dir + sample + '.snvs'
	snv_outf = open(snvs_file, 'w')
	result_file = res_dir + sample + '.txt'
	outf = open(result_file, 'w')
	outf.write('\t'.join(['# cog', 'allele_representative', 'reads', 'reads_with_novelty', 'reads_uniquely_mapping',
						  'reads_uniquely_mapping_with_novelty']) + '\n')

	bam_handle = pysam.AlignmentFile(bam_alignment, 'rb')

	topaligns_file = res_dir + sample + '_topaligns.bam'
	topaligns_file_sorted = res_dir + sample + '_topaligns.sorted.bam'
	topaligns_handle = pysam.AlignmentFile(topaligns_file, "wb", template=bam_handle)

	unialigns_file = res_dir + sample + '_unialigns.bam'
	unialigns_file_sorted = res_dir + sample + '_unialigns.sorted.bam'
	unialigns_handle = pysam.AlignmentFile(unialigns_file, "wb", template=bam_handle)

	for cog, cog_genes in bgc_cog_genes.items():
		read_ascores_per_allele = defaultdict(list)
		read_genes_mapped = defaultdict(set)
		snv_counts = defaultdict(set)
		cog_genes_covered = 0
		rep_alignments = defaultdict(lambda: defaultdict(set))
		with open(ref_fasta) as opff:
			for rec in SeqIO.parse(opff, 'fasta'):
				if rec.id.split('|')[-1] != cog: continue
				g, sample, _ = rec.id.split('|')
				ginfo = comp_gene_info[g]
				gstart = ginfo['relative_start']
				gend = ginfo['relative_end']
				offset = gstart
				gene_length = gend - gstart + 1
				gene_covered_1 = 0
				gene_covered_3 = 0
				for pileupcolumn in bam_handle.pileup(contig=rec.id, start=gstart, stop=gend + 1, stepper="nofilter",
													  truncate=True):
					pos_depth = 0
					for pileupread in pileupcolumn.pileups:
						read = pileupread.alignment
						if pileupread.is_del or pileupread.is_refskip or not read.is_proper_pair: continue
						if read.query_qualities[pileupread.query_position] < 20: continue
						pos_depth += 1
						if pos_depth >= 1:
							gene_covered_1 += 1
							if pos_depth >= 3:
								gene_covered_3 += 1

				gene_coverage_1 = gene_covered_1 / float(gene_length)
				gene_coverage_3 = gene_covered_3 / float(gene_length)
				if gene_coverage_1 < 0.95: continue
				cog_genes_covered += 1

				for read1_alignment, read2_alignment in read_pair_generator(bam_handle, region_string=rec.id,
																			start=gstart, stop=gend):
					read_name = read1_alignment.query_name
					read1_ascore = read1_alignment.tags[0][1]
					read2_ascore = read2_alignment.tags[0][1]
					combined_ascore = read1_ascore + read2_ascore

					snvs = set([])
					g_rep = cog_gene_to_rep[rec.id]

					read1_ref_positions = set(read1_alignment.get_reference_positions())
					read2_ref_positions = set(read2_alignment.get_reference_positions())

					read_intersect = len(read1_ref_positions.intersection(read2_ref_positions))
					min_read_length = min(len(read1_ref_positions), len(read2_ref_positions))
					read_overlap_prop = read_intersect / min_read_length

					min_read1_ref_pos = min(read1_ref_positions)
					read1_referseq = read1_alignment.get_reference_sequence().upper()
					read1_queryseq = read1_alignment.query_sequence
					read1_queryqua = read1_alignment.query_qualities

					min_read2_ref_pos = min(read2_ref_positions)
					read2_referseq = read2_alignment.get_reference_sequence().upper()
					read2_queryseq = read2_alignment.query_sequence
					read2_queryqua = read2_alignment.query_qualities

					alignment_has_indel = False
					mismatch_count = 0
					for b in read1_alignment.get_aligned_pairs(with_seq=True):
						if b[0] == None or b[1] == None:
							alignment_has_indel = True
						elif b[2].islower():
							que_qual = read1_queryqua[b[0]]
							if que_qual >= 30: mismatch_count += 1
					for b in read2_alignment.get_aligned_pairs(with_seq=True):
						if b[0] == None or b[1] == None:
							alignment_has_indel = True
						elif b[2].islower():
							que_qual = read2_queryqua[b[0]]
							if que_qual >= 30: mismatch_count += 1

					for b in read1_alignment.get_aligned_pairs(with_seq=True):
						if b[0] == None or b[1] == None: continue
						if not b[2].islower(): continue
						ref_pos = b[1]
						que_qual = read1_queryqua[b[0]]
						alt_al = read1_queryseq[b[0]].upper()
						ref_al = read1_referseq[b[1] - min_read1_ref_pos].upper()
						# print(alt_al)
						# print(ref_al)
						# print(read1_queryqual[b[0]])
						# print(str(rec.seq).upper()[b[1]])
						# print(b[2])
						assert (ref_al == str(rec.seq).upper()[b[1]])
						assert (alt_al != ref_al)
						if que_qual >= 30 and not alignment_has_indel and mismatch_count <= 5 and read_overlap_prop <= 0.25 and min_read_length >= 75:
							snvs.add(str(rec.id) + '_|_' + str(ref_pos - offset + 1) + '_|_' + ref_al + '_|_' + alt_al)
							snv_counts[
								str(rec.id) + '_|_' + str(ref_pos - offset + 1) + '_|_' + ref_al + '_|_' + alt_al].add(
								read_name)

					for b in read2_alignment.get_aligned_pairs(with_seq=True):
						if b[0] == None or b[1] == None: continue
						if not b[2].islower(): continue
						ref_pos = b[1]
						que_qual = read2_queryqua[b[0]]
						alt_al = read2_queryseq[b[0]].upper()
						ref_al = read2_referseq[b[1] - min_read2_ref_pos].upper()
						# print(alt_al)
						# print(ref_al)
						# print(read2_queryqual[b[0]])
						# print(str(rec.seq).upper()[b[1]])
						# print(b[2])
						assert (ref_al == str(rec.seq).upper()[b[1]])
						assert (alt_al != ref_al)
						if que_qual >= 30 and not alignment_has_indel and mismatch_count <= 5 and read_overlap_prop <= 0.25 and min_read_length >= 75:
							snvs.add(str(rec.id) + '_|_' + str(ref_pos - offset + 1) + '_|_' + ref_al + '_|_' + alt_al)
							snv_counts[
								str(rec.id) + '_|_' + str(ref_pos - offset + 1) + '_|_' + ref_al + '_|_' + alt_al].add(
								read_name)

					read_genes_mapped[read1_alignment.query_name].add(rec.id)
					rep_alignments[rec.id][read1_alignment.query_name].add(tuple([read1_alignment, read2_alignment]))
					read_ascores_per_allele[read1_alignment.query_name].append(
						[g_rep.split('|')[0], g_rep.split('|')[1], combined_ascore, snvs, g])

		if cog_genes_covered / float(len(cog_genes)) < 0.80: continue

		supported_snvs = set([])
		allele_reads = defaultdict(set)
		allele_reads_with_mismatch = defaultdict(set)
		multi_partitioned_reads = set([])
		for read in read_ascores_per_allele:
			top_score = -1000000
			top_score_grep = None
			score_sorted_alignments = sorted(read_ascores_per_allele[read], key=itemgetter(2), reverse=True)
			for i, align in enumerate(score_sorted_alignments):
				g_rep = align[0] + '|' + align[1] + '|' + cog
				if i == 0: top_score = align[2]; top_score_grep = g_rep
				if (i == 0 and align[2] == top_score) and (
						len(score_sorted_alignments) == 1 or align[2] > score_sorted_alignments[i + 1][2]):
					for snv in align[3]:
						if len(snv_counts[snv]) >= 5:
							supported_snvs.add(snv)
				if align[2] == top_score:
					g_map = read_genes_mapped[read].intersection(cog_rep_genes[g_rep])
					g_map_prop = len(g_map) / float(len(cog_rep_genes[g_rep]))
					if g_map_prop < 0.80: continue

					allele_reads[g_rep].add(read)
					for snv in align[3]:
						if len(snv_counts[snv]) >= 5:
							allele_reads_with_mismatch[g_rep].add(read)

				if g_rep != top_score_grep and align[2] == top_score and i > 0:
					multi_partitioned_reads.add(read)

		for snv in supported_snvs:
			snv_outf.write(snv + '\t' + str(len(snv_counts[snv])) + '\n')

		for al in allele_reads:
			for r in allele_reads[al]:
				for pa in rep_alignments[al][r]:
					topaligns_handle.write(pa[0])
					topaligns_handle.write(pa[1])
					if not r in multi_partitioned_reads:
						unialigns_handle.write(pa[0])
						unialigns_handle.write(pa[1])
			outf.write(
				'\t'.join([str(x) for x in [cog, al, len(allele_reads[al]), len(allele_reads_with_mismatch[al]),
											len(allele_reads[al].difference(multi_partitioned_reads)),
											len(allele_reads_with_mismatch[al].difference(
												multi_partitioned_reads))]]) + '\n')

	snv_outf.close()
	outf.close()
	topaligns_handle.close()
	unialigns_handle.close()
	bam_handle.close()

	os.system("samtools sort -@ %d %s -o %s" % (1, unialigns_file, unialigns_file_sorted))
	os.system("samtools index %s" % unialigns_file_sorted)

	os.system("samtools sort -@ %d %s -o %s" % (1, topaligns_file, topaligns_file_sorted))
	os.system("samtools index %s" % topaligns_file_sorted)


# except: pass

def runBowtie2Alignments(bowtie2_reference, paired_end_sequencing_file, bowtie2_outdir, cores, logObject):
	bowtie2_cores = cores
	bowtie2_pool_size = 1
	if cores >= 4:
		bowtie2_cores = 4
		bowtie2_pool_size = int(cores / 4)

	bowtie2_inputs = []
	with open(paired_end_sequencing_file) as opesf:
		for line in opesf:
			line = line.strip()
			sample, frw_read, rev_read = line.split('\t')
			bowtie2_inputs.append(
				[sample, frw_read, rev_read, bowtie2_reference.split('.fasta')[0], bowtie2_outdir, bowtie2_cores,
				 logObject])

	p = multiprocessing.Pool(bowtie2_pool_size)
	p.map(bowtie2_alignment, bowtie2_inputs)
	p.close()


def run_cmd(cmd, logObject, stdout=subprocess.DEVNULL):
	logObject.info('Running the following command: %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))


def bowtie2_alignment(input_args):
	sample, frw_read, rev_read, bgc_cdhit, bowtie2_dir, cores, logObject = input_args

	sam_file = bowtie2_dir + sample + '.sam'
	bam_file = bowtie2_dir + sample + '.bam'
	bam_file_sorted = bowtie2_dir + sample + '.sorted.bam'
	bam_file_filtered = bowtie2_dir + sample + '.filtered.bam'
	bam_file_filtered_sorted = bowtie2_dir + sample + '.filtered.sorted.bam'

	bowtie2_cmd = ['bowtie2', '--very-sensitive', '--no-mixed', '--no-discordant', '--no-unal', '-a', '-x',
				   bgc_cdhit, '-1', frw_read, '-2', rev_read, '-S', sam_file, '-p', str(cores)]
	samtools_view_cmd = ['samtools', 'view', '-h', '-Sb', sam_file, '>', bam_file]
	samtools_sort_cmd = ['samtools', 'sort', '-@', str(cores), bam_file, '-o', bam_file_sorted]
	samtools_index_cmd = ['samtools', 'index', bam_file_sorted]
	samtools_sort_cmd_2 = ['samtools', 'sort', '-@', str(cores), bam_file_filtered, '-o', bam_file_filtered_sorted]
	samtools_index_cmd_2 = ['samtools', 'index', bam_file_filtered_sorted]

	try:
		run_cmd(bowtie2_cmd, logObject)
		run_cmd(samtools_view_cmd, logObject)
		run_cmd(samtools_sort_cmd, logObject)
		run_cmd(samtools_index_cmd, logObject)

		bam_handle = pysam.AlignmentFile(bam_file_sorted, 'rb')
		filt_bam_handle = pysam.AlignmentFile(bam_file_filtered, "wb", template=bam_handle)

		for read in bam_handle.fetch():
			if not read.is_proper_pair or read.is_supplementary: continue
			filt_bam_handle.write(read)

		filt_bam_handle.close()
		bam_handle.close()

		run_cmd(samtools_sort_cmd_2, logObject)
		run_cmd(samtools_index_cmd_2, logObject)

		os.system(
			"rm -f %s %s %s %s %s" % (sam_file, bam_file, bam_file_sorted, bam_file_filtered, bam_file_sorted + '.bai'))
	except:
		os.system('rm -f %s/%s*' % (bowtie2_dir, sample))
		raise RuntimeError()


def read_pair_generator(bam, region_string=None, start=None, stop=None):
	"""
    Function taken from: https://www.biostars.org/p/306041/
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
	read_dict = defaultdict(lambda: [None, None])
	for read in bam.fetch(region_string, start=start, stop=stop):
		if not read.is_proper_pair or read.is_supplementary:
			continue
		qname = read.query_name
		if qname not in read_dict:
			if read.is_read1:
				read_dict[qname][0] = read
			else:
				read_dict[qname][1] = read
		else:
			if read.is_read1:
				yield read, read_dict[qname][1]
			else:
				yield read_dict[qname][0], read
			del read_dict[qname]