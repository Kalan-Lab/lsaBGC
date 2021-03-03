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

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_PLOTTING = lsaBGC_main_directory + 'plotCogConservation.R'
RSCRIPT_FOR_POPGEN = lsaBGC_main_directory + 'popGenStatisticsUsingPegas.R'

def multiProcess(input):
	input_cmd = input[:-1]
	logObject = input[-1]
	logObject.info('Running the following command: %s' % ' '.join(input_cmd))
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(input_cmd))
	except:
		logObject.warning('Had an issue running: %s' % ' '.join(input_cmd))

def readInBGCGenbanksPerGCF(gcf_specs_file, logObject):
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
				logObject.error("More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
				raise RuntimeError("More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
			sample, gbk = line.split('\t')
			try:
				assert (is_genbank(gbk))
				bgc_id = sample
				if sample_index[sample] > 0:
					bgc_id = sample + '_' + str(sample_index[sample] + 1)
				sample_index[sample] += 1

				bgc_genes_full, bgc_info_full = parseGenbanks(gbk, bgc_id)

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

def createItolBGCSeeTrack(bgc_genes, gene_to_cog, cog_to_color, comp_gene_info, dataset_label, outdir, logObject):
	track_file = outdir + 'BGCs_Visualization.iTol.txt'
	track_handle = open(track_file, 'w')

	logObject.info("Writing track file to: %s" % track_file)
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
				cog_lengths[cog].append(gend-gstart)
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
				bgc_end = int(printlist[1])
				for gene_string in printlist[2:]:
					gene_info = gene_string.split('|')
					new_shape = None
					if gene_info[0] == 'TR':
						new_shape = 'TL'
					elif gene_info[0] == 'TL':
						new_shape = 'TR'
					new_gstart = int(bgc_end) - int(gene_info[2])
					new_gend = int(bgc_end) - int(gene_info[1])
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

	cog_order_scores = {}
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
				cog_lengths[cog].append(gend-gstart)
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

def constructCodonAlignments(bgc_sample, cog_genes, comp_gene_info, outdir, cores, logObject):

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
			sample_counts = defaultdict(int)
			sample_sequences = {}
			for gene in cog_genes[cog]:
				gene_info = comp_gene_info[gene]
				bgc_id = gene_info['bgc_name']
				sample_id = bgc_sample[bgc_id]
				nucl_seq = gene_info['nucl_seq']
				prot_seq = gene_info['prot_seq']
				sample_counts[sample_id] += 1
				sample_sequences[sample_id] = tuple([nucl_seq, prot_seq])
			samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
			# check that cog is single-copy-core
			if len(samples_with_single_copy.symmetric_difference(all_samples)) > 0: continue
			logObject.info('Homolog group %s detected as SCC across samples (not individual BGCs).' % cog)
			inputs.append([cog, sample_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, logObject])


		p = multiprocessing.Pool(cores)
		p.map(create_codon_msas, inputs)

		return codo_alg_dir
	except:
		logObject.error("Issues with create protein/codon alignments of SCC homologs for BGC.")
		raise RuntimeError("Issues with create protein/codon alignments of SCC homologs for BGC.")


def create_codon_msas(inputs):
	cog, sample_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, logObject = inputs

	cog_nucl_fasta = nucl_seq_dir + '/' + cog + '.fna'
	cog_prot_fasta = prot_seq_dir + '/' + cog + '.faa'
	cog_prot_msa = prot_alg_dir + '/' + cog + '.msa.faa'
	cog_codo_msa = codo_alg_dir + '/' + cog + '.msa.fna'

	cog_nucl_handle = open(cog_nucl_fasta, 'w')
	cog_prot_handle = open(cog_prot_fasta, 'w')
	for s in sample_sequences:
		cog_nucl_handle.write('>' + s + '\n' + str(sample_sequences[s][0]) + '\n')
		cog_prot_handle.write('>' + s + '\n' + str(sample_sequences[s][1]) + '\n')
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
				assert(is_genbank(bgc_genbank))
				assert(is_fasta(bgc_proteome))
				comprehensive_bgcs.append(tuple([sample, bgc_genbank, bgc_proteome]))
		return comprehensive_bgcs
	except:
		logObject.error("Had difficulties in processing file listing comprehensive set of BGCs and paths to genbanks and extracted proteins.")
		raise RuntimeError("Had difficulties in processing file listing comprehensive set of BGCs and paths to genbanks and extracted proteins.")

def runHMMScanAndAssignBGCsToGCF(comprehensive_bgcs, concat_hmm_profiles, scc_homologs, outdir, cores, logObject):
	search_res_dir = os.path.abspath(outdir + 'HMMScan_Results') + '/'
	if not os.path.isdir(search_res_dir): os.system('mkdir %s' % search_res_dir)

	hmmscan_cmds = []
	for i, cb_tuple in enumerate(comprehensive_bgcs):
		sample, bgc_genbank, bgc_proteome = cb_tuple
		result_file = search_res_dir + str(i) + '.txt'
		hmmscan_cmd = ['hmmscan', '--max', '--cpu', '1', '--tblout', result_file, concat_hmm_profiles, bgc_proteome, logObject]
		hmmscan_cmds.append(hmmscan_cmd)

	p = multiprocessing.Pool(cores)
	p.map(multiProcess, hmmscan_cmds)
	p.close()

	sample_cogs = defaultdict(set)
	bgc_cogs = defaultdict(set)
	sample_bgcs = defaultdict(set)
	sample_cog_proteins = defaultdict(lambda: defaultdict(set))
	for i, cb_tuple in enumerate(comprehensive_bgcs):
		sample, bgc_genbank, bgc_proteome = cb_tuple
		sample_bgcs[sample].add(bgc_genbank)
		result_file = search_res_dir + str(i) + '.txt'
		assert(os.path.isfile(result_file))
		with open(result_file) as orf:
			for line in orf:
				if line.startswith("#"): continue
				line = line.strip()
				ls = line.split()
				cog = ls[0]
				samp, bgc, protein_id = ls[2]
				eval = float(ls[4])
				if eval <= 1e-5:
					sample_cogs[sample].add(cog)
					bgc_cogs[bgc_genbank].add(cog)
					sample_cog_proteins[sample][cog].add(protein_id)

	expanded_orthofinder_matrix_file = outdir + 'Orthogroups.expanded.csv'
	expanded_bgc_list_file = outdir + 'GCF_Expanded.txt'

	expanded_orthofinder_matrix_handle = open(expanded_orthofinder_matrix_file, 'w')
	expanded_bgc_list_handle = open(expanded_bgc_list_file, 'w')

	for s in sample_cogs:
		if len(scc_homologs.difference(sample_cogs[s])) == 0:
			for bgc_gbk in sample_bgcs:
				if len(bgc_cogs[bgc_gbk]) >= 3 and len(scc_homologs.intersection(bgc_cogs[bgc_gbk])) >= 1:
					expanded_bgc_list_handle.write('\t'.join([sample, bgc_gbk]) + '\n')




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

		logObject.info("Successfully created profile HMMs for each homolog group. Now beginning concatenation into single file.")
		concatenated_profile_HMM = outdir + 'All_GCF_Homologs.hmm'
		os.system('rm -f %s' % concatenated_profile_HMM)
		for f in os.listdir(prot_hmm_dir):
			os.system('cat %s >> %s' % (prot_hmm_dir + f, concatenated_profile_HMM))

		hmmpress_cmd = ['hmmpress', concatenated_profile_HMM]
		logObject.info('Running hmmpress on concatenated profiles with the following command: %s' % ' '.join(hmmpress_cmd))
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
		logObject.info("New phylogeny with an additional %d leafs to reflect samples with multiple BGCs can be found at: %s." % (number_of_added_leaves, output_phylogeny))
	except:
		logObject.error("Had difficulties properly editing phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")
		raise RuntimeError("Had difficulties properly editing phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")

def readInBGCGenbanksComprehensive(bgc_specs_file, logObject):
	bgc_sample = {}
	bgc_product = {}
	bgc_genes = {}
	all_genes = set([])
	with open(bgc_specs_file) as obsf:
		for i, line in enumerate(obsf):
			line = line.strip()
			try:
				assert (len(line.split('\t')) == 2)
			except:
				logObject.error(
					"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
				raise RuntimeError(
					"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
			sample, gbk = line.split('\t')
			try:
				assert (is_genbank(gbk))
				bgc_genes_full, bgc_info_full = parseGenbanks(gbk, gbk)
				bgc_product[gbk] = [x['product'] for x in bgc_info_full]
				bgc_genes[gbk] = set(bgc_genes_full.keys())
				all_genes = all_genes.union(bgc_genes[gbk])
				bgc_sample[gbk] = sample
				logObject.info("Incorporating genbank %s for sample %s into analysis." % (gbk, sample))
			except:
				logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis.")
				raise RuntimeWarning("Unable to validate %s as Genbank. Skipping ...")

	return [bgc_sample, bgc_product, bgc_genes, all_genes]


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
	try:
		populations = set([])
		sample_population = {}
		with open(pop_specs_file) as opsf:
			for line in opsf:
				line = line.strip()
				sample, population = line.split('\t')
				sample_population[sample] = population
				populations.add(population)
		logObject.info("Successfully parsed population specifications file. There are %d populations" % len(population))
		return sample_population
	except Exception as e:
		logObject.error('Failed to upload to ftp: ' + str(e))
		logObject.error("Had issues parsing parsed population specifications file. Exiting now.")
		raise RuntimeError("Had issues parsing parsed population specifications file. Exiting now.")

def writeLociFileForPegas(codon_msa, outfile, logObject, sample_populations=None):
	outfile_handle = open(outfile, 'w')
	number_of_loci = 0
	with open(codon_msa) as ocmf:
		for i, rec in enumerate(SeqIO.parse(ocmf, 'fasta')):
			if i == 0:
				number_of_loci = len(str(rec.seq))
				header = ['Sample'] + ['Pos' + str(p) for p in range(1, len(str(rec.seq))+1)]
				if sample_populations:
					header += ['population']
				outfile_handle.write('\t'.join(header) + '\n')
			sample_data = [rec.id] + '\t' + [p.replace('-', 'NA').upper() for p in str(rec.seq)]
			if sample_populations:
				try:
					sample_data += [sample_populations[rec.id]]
				except:
					logObject.warning("Didn't find population for sample %s, filling in population as NA." % rec.id)
					sample_data += ['NA']
			outfile_handle.write('\t'.join(sample_data) + '\n')

	outfile_handle.close()
	return [number_of_loci, outfile]

def parseCodonAlignmentStats(cog, codon_alignment_fasta, comp_gene_info, cog_genes, cog_order_index, cog_median_copy_numbers, plots_dir, popgen_dir, logObject, sample_population=None):
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
			if comp_gene_info[gene_id]['core_overlap']: core_counts['core'] += 1
			else: core_counts['auxiliary'] += 1
			products.add(comp_gene_info[gene_id]['product'])
			real_pos = 1
			seqs.append(list(str(rec.seq)))
			codons = [str(rec.seq)[i:i+3] for i in range(0, len(str(rec.seq)), 3)]
			num_codons = len(codons)
			bgc_codons[rec.id] = codons
			samples.add(sample_id)
			for msa_pos, bp in enumerate(str(rec.seq)):
				if bp != '-':
					gene_locs[gene_id][real_pos] = msa_pos+1
					real_pos += 1
			gene_lengths.append(len(str(rec.seq).replace('-', '')))

	is_core = False
	if float(core_counts['core'])/sum(core_counts.values()) >= 0.8: is_core = True

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
		maj_allele_freq = float(maj_allele_count)/tot_count
		position_plot_handle.write('\t'.join([str(x) for x in [i+1, num_seqs, num_alleles, num_gaps, maj_allele_freq]]) + '\n')
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
			relative_start = domain_start-gene_start
			assert(len(gene_locs[gene])+3 >= (domain_end-gene_start))
			relative_end = min([len(gene_locs[gene]), domain_end-gene_start])
			domain_range = range(relative_start, relative_end)
			for pos in domain_range:
				msa_pos = gene_locs[gene][pos+1]
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
				old_pos = pos-1
			if pos-1 != old_pos:
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
	core_codons = 0
	total_variable_codons = 0
	nonsynonymous_sites = 0
	synonymous_sites = 0
	for cod_index in range(0, num_codons):
		first_bp = (cod_index+1)*3
		aa_count = defaultdict(int)
		cod_count = defaultdict(int)
		core = True
		for bgc in bgc_codons:
			cod = bgc_codons[bgc][cod_index]
			if '-' in cod or 'N' in cod:
				core = False
			else:
				cod_obj = Seq(cod)
				aa_count[str(cod_obj.translate().seq)] += 1
				cod_count[cod] += 1
		maj_allele_count = max(cod_count.values())
		tot_valid_codons = sum(cod_count.values())
		residues = len(aa_count.keys())
		maj_allele_freq = float(maj_allele_count)/tot_valid_codons

		nonsyn_flag = False; syn_flag = False
		if maj_allele_freq <= 0.9:
			total_variable_codons += 1
			if len(cod_count.keys()) > 1:
				if residues > 1: nonsynonymous_sites += 1; nonsyn_flag = True
				else: synonymous_sites += 1; syn_flag = True

		if core: core_codons += 1
		if nonsyn_flag or syn_flag:
			type = 'S'
			if nonsyn_flag: type = 'NS'
			popgen_plot_handle.write('\t'.join([str(x) for x in [first_bp, type]]) + '\n')
	popgen_plot_handle.close()
	dn_ds = "NA"
	if synonymous_sites > 0: dn_ds =  float(nonsynonymous_sites) / synonymous_sites

	os.system('Rscript %s %s %s %s %s' % (RSCRIPT_FOR_PLOTTING, domain_plot_file, position_plot_file, popgen_plot_file, plot_pdf_file))

	popgen_loci_file = popgen_dir + cog + '_loci.txt'
	number_of_loci = writeLociFileForPegas(codon_alignment_fasta, popgen_loci_file, logObject, sample_population=sample_population)

	header = ['cog', 'annotation', 'cog_order_index', 'cog_median_copy_count', 'is_core', 'median_gene_length',
			  'bgcs_with_cog',
			  'samples_with_cog', 'Tajimas_D', 'core_codons', 'total_variable_codons', 'nonsynonymous_codons',
			  'synonymous_codons',
			  'dn_ds', 'all_domains']
	if sample_population:
		header += ['populations_with_cog', 'max_population_specificity', 'Fst']



	return (
		['; '.join(products), is_core, median_gene_length, len(seqs), len(samples), cs_stats['D'], core_codons, total_variable_codons,
		 nonsynonymous_sites, synonymous_sites,
		 dn_ds, cd.num_codons_eff, cd.num_pol_NS, cd.num_pol_S,
		 effective_dn_ds, tot_phylogenetic_breadth, '; '.join(differential_domains), '; '.join(all_domains)])



def parseGenbanks(gbk, bgc_name):
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
			for feature in rec.features:
				if feature.type in domain_feature_types:
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
					product = feature.qualifiers.get('product')[0]
					contig_edge = feature.qualifiers.get('contig_edge')[0]
					full_sequence = str(rec.seq)
					bgc_info.append({'detection_rule': detection_rule, 'product': product, 'contig_edge': contig_edge,
									 'full_sequence': str(rec.seq)})
				elif feature.type == 'proto_core':
					core_start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					core_end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					core_positions = core_positions.union(set(range(core_start, core_end + 1)))

	if len(bgc_info) == 0:
		bgc_info = [{'detection_rule': 'NA', 'product': 'NA', 'contig_edge': 'NA', 'full_sequence': full_sequence}]

	genes = {}
	with open(gbk) as ogbk:
		for rec in SeqIO.parse(ogbk, 'genbank'):
			for feature in rec.features:
				if feature.type == "CDS":
					lt = feature.qualifiers.get('locus_tag')[0]
					prot_seq = feature.qualifiers.get('translation')[0]
					start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					direction = str(feature.location).split('(')[1].split(')')[0]
					product = feature.qualifiers.get('product')[0]
					grange = set(range(start, end + 1))
					gene_domains = []
					for d in domains:
						drange = set(range(d['start'], d['end'] + 1))
						if len(drange.intersection(grange)) > 0:
							gene_domains.append(d)
					nucl_seq = full_sequence[start:end]
					core_overlap = False
					if len(grange.intersection(core_positions)) > 0: core_overlap = True
					if direction == '-':
						nucl_seq = str(Seq(full_sequence[start:end]).reverse_complement())
					genes[lt] = {'bgc_name': bgc_name, 'start': start, 'end': end, 'direction': direction,
								 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
								 'gene_domains': gene_domains, 'core_overlap': core_overlap}
	return ([genes, bgc_info])

def calculateBGCPairwiseRelations(bgc_genes, gene_to_cog, avg_nz_cog_counts, outdir, logObject):
	pair_relations_txt_file = outdir + 'bgc_pair_relationships.txt'
	pair_relations_mci_file = outdir + 'bgc_pair_relationships.mci'
	pair_relations_tab_file = outdir + 'bgc_pair_relationships.tab'

	bgc_cogs = defaultdict(set)
	for bgc in bgc_genes:
		for gene in bgc_genes[bgc]:
			if gene in gene_to_cog:	bgc_cogs[bgc].add(gene_to_cog[gene])

	try:
		prf_handle = open(pair_relations_txt_file, 'w')
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
		logObject.info("Calculated pairwise relations and wrote to: %s" % pair_relations_txt_file)
	except:
		logObject.error("Problem with creating relations file between pairs of BGCs.")
		raise RuntimeError("Problem with creating relations file between pairs of BGCs.")

	mcxload_cmd = ['mcxload', '-abc', pair_relations_txt_file, '--stream-mirror', '-write-tab', pair_relations_tab_file,
				   '-o', pair_relations_mci_file]

	logObject.info('Running the following command: %s' % ' '.join(mcxload_cmd))
	try:
		subprocess.call(' '.join(mcxload_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(mcxload_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(mcxload_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(mcxload_cmd))
	logObject.info('Converted format of pair relationship file via mxcload.')

	return bgc_cogs, pairwise_relations, pair_relations_txt_file, pair_relations_mci_file, pair_relations_tab_file

def runMCLAndReportGCFs(mip, outdir, sf_handle, pairwise_relations, pair_relations_mci_file, pair_relations_tab_file, bgc_cogs, bgc_product, bgc_sample, inflation_testing, cores, logObject):
	relations_mcl_file = outdir + 'mcl.' + str(mip).replace('.', '_') + '.out'
	mcxdump_out_file = outdir + 'final_mcl.' + str(mip).replace('.', '_') + '.out'
	mcl_cmd = ['mcl', pair_relations_mci_file, '-I', str(mip), '-o', relations_mcl_file, '-te', str(cores)]
	mcxdump_cmd = ['mcxdump', '-icl', relations_mcl_file, '-tabr', pair_relations_tab_file, '-o', mcxdump_out_file]

	logObject.info('Running MCL and MCXDUMP with inflation parameter set to %f' % mip)
	logObject.info('Running the following command: %s' % ' '.join(mcl_cmd))
	try:
		subprocess.call(' '.join(mcl_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
		assert(os.path.isfile(relations_mcl_file) and os.path.getsize(relations_mcl_file) > 100)
		logObject.info('Successfully ran: %s' % ' '.join(mcl_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(mcl_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(mcl_cmd))

	sys.stderr.write('Dumping results in human-readable format ...')
	logObject.info('Running the following command: %s' % ' '.join(mcxdump_cmd))
	try:
		subprocess.call(' '.join(mcxdump_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
		assert(os.path.isfile(mcxdump_out_file) and os.path.getsize(mcxdump_out_file) > 100)
		logObject.info('Successfully ran: %s' % ' '.join(mcl_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(mcxdump_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(mcxdump_cmd))

	logObject.info('Successfully ran MCL and MCXDUMP with inflatimulti_same_sampleon parameter set to %f!' % mip)

	clustered_bgcs = set([])
	with open(mcxdump_out_file) as omo:
		for j, gcf in enumerate(omo):
			gcf = gcf.strip()
			gcf_mems = gcf.split()
			if len(gcf_mems) < 2: continue
			diffs = set([])
			samp_counts = defaultdict(int)
			samp_ogs = defaultdict(set)
			products = set([])
			for a, bgc1 in enumerate(gcf_mems):
				samp1 = bgc_sample[bgc1]
				samp_counts[samp1] += 1
				for prod in bgc_product[bgc1]: products.add(prod)
				samp_ogs[samp1] = samp_ogs[samp1].union(bgc_cogs[bgc1])
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
			stdev = "NA"
			mean = "NA"
			try:
				mean = statistics.mean(num_ogs)
				stdev = statistics.stdev(num_ogs)
			except:
				pass
			gcf_stats = ['GCF_' + str(j+1), len(gcf_mems), multi_same_sample, len(scc), mean, stdev, min(diffs),
						 max(diffs), '; '.join(products)]
			if inflation_testing:
				gcf_stats = [mip] + gcf_stats
			sf_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
	singleton_bgcs = set([])
	for bgc in bgc_cogs:
		if not bgc in clustered_bgcs: singleton_bgcs.add(bgc)
	singleton_stats = ['singletons', len(singleton_bgcs)] + (['NA'] * 7)
	if inflation_testing:
		singleton_stats = [mip] + singleton_stats
	sf_handle.write('\t'.join([str(x) for x in singleton_stats]) + '\n')

	if not inflation_testing:
		logObject.info("Writing list of BGCs for each GCF, which will be used as input for downstream programs in the suite!")
		gcf_listing_dir = outdir + 'GCF_Listings/'
		if not os.path.isdir(gcf_listing_dir): os.system('mkdir %s' % gcf_listing_dir)
		with open(mcxdump_out_file) as omo:
			for j, gcf in enumerate(omo):
				gcf = gcf.strip()
				gcf_mems = gcf.split()
				if len(gcf_mems) < 2: continue
				outf_list = open(gcf_listing_dir + 'GCF_' + str(j + 1) + '.txt', 'w')
				for bgc in gcf_mems:
					sname = bgc_sample[bgc]
					gbkpath = bgc
					outf_list.write('%s\t%s\n' % (sname, gbkpath))
				outf_list.close()
		logObject.info("Successfully wrote lists of BGCs for each GCF.")
	return

def parseOrthoFinderMatrix(orthofinder_matrix_file, all_gene_lts):
	"""
	:param orthofinder_matrix: OrthoFinderV2 matrix Orthogroups.csv file, should also include singleton orthologroups
	:param all_gene_lts: Set of all the relevant gene locus tag identifiers found in BGC Genbanks
	:return: dictionary mapping gene locus tags to CoGs
	"""
	gene_to_cog = {}
	cog_genes = defaultdict(set)
	cog_median_gene_counts = defaultdict(lambda: 'NA')
	with open(orthofinder_matrix_file) as ofm:
		for i, line in enumerate(ofm):
			if i == 0: continue
			line = line.strip()
			ls = line.split('\t')
			cog = ls[0]
			flag_in_bgc = False
			for sgs in ls[1:]:
				for g in sgs.split(', '):
					if g in all_gene_lts:
						flag_in_bgc = True
						gene_to_cog[g] = cog
						cog_genes[cog].add(g)
			if flag_in_bgc:
				gene_counts = []
				for sgs in ls[1:]:
					gene_counts.append(len(sgs.split(', ')))
				cog_median_gene_counts[cog] = statistics.median(gene_counts)
	return ([gene_to_cog, cog_genes, cog_median_gene_counts])

def readInAssemblyListing(assembly_listing_file, logObject):
	sample_assembly_paths = {}
	try:
		with open(assembly_listing_file) as oalf:
			for line in oalf:
				line = line.strip()
				sample, assembly = line.split('\t')
				sample = sample.replace(' ', '_').replace('|', '_').replace('"', '_').replace("'", '_')
				try:
					assert (is_fasta(assembly))
					sample_assembly_paths[sample] = assembly
				except:
					logObject.info(
						'Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.' % sample)
					sys.stderr.write(
						'Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.\n' % sample)
		assert (len(sample_assembly_paths) >= 2)
		return (sample_assembly_paths)
	except:
		raise RuntimeError('Input file listing the location of assemblies for samples leads to incorrect path was provided. Exiting now ...')

def runProkka(sample_assemblies, prokka_outdir, prokka_proteomes, prokka_genbanks, prokka_load_code, dry_run_flag,
			  lineage, cores, logObject):
	task_file = prokka_outdir + 'prokka.cmds'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing Prokka commands to following: %s" % task_file)
		if os.path.isfile(task_file): os.system('rm -f %s' % task_file)
	alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
	possible_locustags = list(itertools.product(alphabet, repeat=3))
	prokka_cmds = []
	for i, sample in enumerate(sample_assemblies):
		sample_assembly = sample_assemblies[sample]
		sample_outdir = prokka_outdir + sample + '/'
		prokka_cmd = [prokka_load_code, 'prokka', '--cpus', '1', '--outdir', sample_outdir, '--prefix', sample,
					  '--genus', lineage, '--locustag', ''.join(list(possible_locustags[i])),
					  sample_assembly, ';', 'mv', sample_outdir + sample + '.gbf', prokka_genbanks + sample + '.gbk',
					  ';', 'mv',
					  sample_outdir + sample + '.faa', prokka_proteomes]
		prokka_cmds.append(prokka_cmd + [logObject])
		if dry_run_flag:
			task_handle = open(task_file, 'a+')
			task_handle.write(' '.join(prokka_cmd) + '\n')
			task_handle.close()

	if not dry_run_flag:
		p = multiprocessing.Pool(cores)
		p.map(multiProcess, prokka_cmds)
		p.close()

def extractBGCProteomes(s, bgc_genbank, bgc_proteomes_outdir, logObject):
	b = bgc_genbank.split('/')[-1].split('.gbk')[0]
	bgc_genbank_proteome = bgc_proteomes_outdir + s + '_BGC_' + b + '.faa'
	try:
		bgc_genbank_proteome_handle = open(bgc_genbank_proteome, 'w')
		with open(bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == "CDS":
						lt = feature.qualifiers.get('locus_tag')[0]
						prot_seq = feature.qualifiers.get('translation')[0]
						prot_id = s + '|' + b + '|' + lt
						bgc_genbank_proteome_handle.write('>' + prot_id + '\n' + str(prot_seq) + '\n')
		bgc_genbank_proteome_handle.close()
		return bgc_genbank_proteome
	except:
		logObject.error("Had problems extracting protein sequences of CDS features in genbank %s" % bgc_genbank)
		raise RuntimeError("Had problems extracting protein sequences of CDS features in genbank %s" % bgc_genbank)

def runAntiSMASHFromAssemblies(sample_assemblies, antismash_outdir, antismash_load_code, dry_run_flag, cores, logObject, barebone=True):
	task_file = antismash_outdir + 'antismash.cmds'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing AntiSMASH commands to following: %s" % task_file)
		if os.path.isfile(task_file): os.system('rm -f %s' % task_file)

	asm_cores = 1
	pool_size = 1
	if cores < 4:
		asm_cores = cores
	else:
		asm_cores = 4
		pool_size = int(cores / 4)
	antismash_cmds = []
	for sample in sample_assemblies:
		sample_resdir = antismash_outdir + sample + '/'
		antismash_cmd = []
		if barebone:
			antismash_cmd = [antismash_load_code, 'antismash', '--taxon', 'bacteria', '--genefinding-tool',
							'prodigal-m', '--output-dir', sample_resdir, '-c', str(asm_cores), sample_assemblies[sample]]
		else:
			antismash_cmd = [antismash_load_code, 'antismash', '--taxon', 'bacteria', '--genefinding-tool', 'prodigal',
							 '--output-dir', sample_resdir, '--fullhmmer', '--asf', '--cb-general', '--cb-subclusters',
							 '--cb-knownclusters', '--cf-create-clusters', '-c', str(asm_cores), sample_assemblies[sample]]
		antismash_cmds.append(antismash_cmd + [logObject])
		if dry_run_flag:
			task_handle = open(task_file, 'a+')
			task_handle.write(' '.join(antismash_cmd) + '\n')
			task_handle.close()

	if not dry_run_flag:
		p = multiprocessing.Pool(pool_size)
		p.map(multiProcess, antismash_cmds)
		p.close()

def runAntiSMASH(prokka_genbanks_dir, antismash_outdir, antismash_load_code, dry_run_flag, cores, logObject):
	task_file = antismash_outdir + 'antismash.cmds'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing AntiSMASH commands to following: %s" % task_file)
		if os.path.isfile(task_file): os.system('rm -f %s' % task_file)

	asm_cores = 1
	pool_size = 1
	if cores < 4:
		asm_cores = cores
	else:
		asm_cores = 4
		pool_size = int(cores / 4)
	antismash_cmds = []
	for sample_gbk in os.listdir(prokka_genbanks_dir):
		sample = sample_gbk.split('.gbk')[0]
		sample_resdir = antismash_outdir + sample + '/'
		antismash_cmd = [antismash_load_code, 'antismash', '--taxon', 'bacteria', '--genefinding-tool', 'prodigal',
						'--output-dir', sample_resdir, '--fullhmmer', '--asf', '--cb-general', '--cb-subclusters',
						'--cb-knownclusters', '--cf-create-clusters', '-c', str(asm_cores),
							 prokka_genbanks_dir + sample_gbk]
		antismash_cmds.append(antismash_cmd + [logObject])
		if dry_run_flag:
			task_handle = open(task_file, 'a+')
			task_handle.write(' '.join(antismash_cmd) + '\n')
			task_handle.close()

	if not dry_run_flag:
		p = multiprocessing.Pool(pool_size)
		p.map(multiProcess, antismash_cmds)
		p.close()

def runOrthoFinder(prokka_proteomes_dir, orthofinder_outdir, orthofinder_load_code, dry_run_flag, cores, logObject):
	task_file = orthofinder_outdir + 'orthofinder.cmd'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing OrthoFinder command to following: %s" % task_file)

	orthofinder_cmd = [orthofinder_load_code, 'orthofinder', '-f', prokka_proteomes_dir, '-t', str(cores)]
	if dry_run_flag:
		task_handle = open(task_file, 'w')
		task_handle.write(' '.join(orthofinder_cmd) + '\n')
		task_handle.close()
	else:
		try:
			logObject.info('Running the following command: %s' % ' '.join(orthofinder_cmd))
			subprocess.call(' '.join(orthofinder_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran OrthoFinder!')
			tmp_orthofinder_dir = os.path.abspath(
				[prokka_proteomes_dir + f for f in os.listdir(prokka_proteomes_dir) if f.startswith('Results')][
					0]) + '/'
			os.system('mv %s %s' % (tmp_orthofinder_dir, orthofinder_outdir))
		except:
			logObject.warning('Had an issue running: %s' % ' '.join(orthofinder_cmd))
			raise RuntimeError('Had an issue running: %s' % ' '.join(orthofinder_cmd))

def is_fasta(fasta):
	try:
		with open(fasta) as of:
			SeqIO.parse(of, 'fasta')
		return True
	except:
		return False


def is_genbank(fasta):
	try:
		with open(fasta) as of:
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
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		sys.stderr.write(pn + ': ' + str(pv) + '\n')


def logParametersToFile(parameter_file, parameter_names, parameter_values):
	parameter_handle = open(parameter_file, 'w')
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		parameter_handle.write(pn + ': ' + str(pv) + '\n')
	parameter_handle.close()

def logParametersToObject(logObject, parameter_names, parameter_values):
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		logObject.info(pn + ': ' + str(pv))


def assignColorsToCOGs(cogs):
	"""
	:param cogs: set of CoGs.
	:return: dictionary mapping each CoG to a hex color value.
	"""

	# read in list of colors
	dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
	colors_file = dir_path + 'colors_200.txt'
	colors = []
	with open(colors_file) as ocf:
		colors = [x.strip() for x in ocf.readlines()]
	random.shuffle(colors)

	cog_to_color = {}
	for i, c in enumerate(set(cogs)):
		cog_to_color[c] = colors[i]
	return(cog_to_color)

