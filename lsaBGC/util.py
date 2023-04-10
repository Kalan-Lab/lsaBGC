import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import logging
import subprocess
import statistics
from operator import itemgetter
from collections import defaultdict
import traceback
import multiprocessing
import copy
from scipy import stats
from ete3 import Tree
import itertools
import math
import numpy as np
import gzip
import pathlib
import operator

valid_alleles = set(['A', 'C', 'G', 'T'])
curr_dir = os.path.abspath(pathlib.Path(__file__).parent.resolve()) + '/'
main_dir = '/'.join(curr_dir.split('/')[:-2]) + '/'


def writeRefinedProteomes(s, sample_bgcs, refined_proteomes_outdir, logObject):
	try:
		refined_proteome_handle = open(refined_proteomes_outdir + s + '.faa', 'w')
		for bgc in sample_bgcs:
			with open(bgc) as obgc:
				for rec in SeqIO.parse(obgc, 'genbank'):
					for feature in rec.features:
						if feature.type == "CDS":
							lt = feature.qualifiers.get('locus_tag')[0]
							prot_seq = feature.qualifiers.get('translation')[0]
							refined_proteome_handle.write('>' + lt + '\n' + str(prot_seq) + '\n')
		refined_proteome_handle.close()
	except:
		logObject.warning(
			"Had issues writing sample %s's BGC-specific proteome file. Will be ignored in OrthoFinder analysis.")


def getSampleRetentionSet(sample_retention_file):
	sample_retention_set = None
	if sample_retention_file:
		sample_retention_set = set([])
		with open(sample_retention_file) as osrf:
			for line in osrf:
				line = line.strip()
				sample_retention_set.add(line)

	return sample_retention_set


def determineOutliersByGeneLength(gene_sequences, logObject):
	filtered_gene_sequences = {}
	try:
		og_gene_nucl_seq_lens = []
		for g in gene_sequences:
			sample, gene = g.split('|')
			if len(gene.split('_')[0]) == 3:
				gene_nucl_seq = gene_sequences[g][0]
				gene_nucl_seq_len = len(gene_nucl_seq)
				og_gene_nucl_seq_lens.append(gene_nucl_seq_len)

		median_gene_nucl_seq_lens = statistics.median(og_gene_nucl_seq_lens)
		mad_gene_nucl_seq_lens = max(stats.median_abs_deviation(og_gene_nucl_seq_lens, scale="normal"), 25)

		for g in gene_sequences:
			gene_nucl_seq = gene_sequences[g][0]
			gene_nucl_seq_len = len(gene_nucl_seq)
			if abs(gene_nucl_seq_len - median_gene_nucl_seq_lens) <= mad_gene_nucl_seq_lens:
				filtered_gene_sequences[g] = gene_sequences[g]
	except:
		logObject.warning(
			"Unable to filter gene sequences to remove outliers, possibly because there are too few sequences.")
		filtered_gene_sequences = gene_sequences
	return filtered_gene_sequences


def determineNonUniqueRegionsAlongCodonAlignment(outdir, initial_sample_prokka_data, codon_alignments_file, cpus=1,
												 logObject=None):
	"""
	Wrapper function to determine regions along
	"""
	outdir = os.path.abspath(outdir) + '/'
	prot_seq_dir = outdir + 'Protein_Sequences/'

	if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)

	try:
		gcf_protein_ids = set([])
		gcf_protein_to_hg = {}
		gene_pos_to_msa_pos = defaultdict(lambda: defaultdict(dict))
		with open(codon_alignments_file) as ocaf:
			for line in ocaf:
				line = line.strip()
				hg, cmsa_fasta = line.split('\t')
				with open(cmsa_fasta) as ocf:
					for rec in SeqIO.parse(ocf, 'fasta'):
						samp, gene_id = rec.id.split('|')
						if len(gene_id.split('_')[0]) != 3: continue
						gcf_protein_to_hg[gene_id] = hg
						gcf_protein_ids.add(gene_id)
						real_pos = 1
						for msa_pos, bp in enumerate(str(rec.seq)):
							if bp != '-':
								gene_pos_to_msa_pos[hg][gene_id][real_pos] = msa_pos + 1
								real_pos += 1

		all_gcf_proteins_fasta_file = outdir + 'All_GCF_Proteins.faa'
		all_gcf_proteins_fasta_db = outdir + 'All_GCF_Proteins'
		all_comp_gcf_proteins_fasta_file = outdir + 'Complement_All_GCF_Proteins.faa'
		diamond_outfmt6_result_file = outdir + 'Diamond_CompProts_vs_Prots.txt'

		all_gcf_proteins_fasta_handle = open(all_gcf_proteins_fasta_file, 'w')
		all_comp_gcf_proteins_fasta_handle = open(all_comp_gcf_proteins_fasta_file, 'w')

		original_samples = set([])
		for sample in initial_sample_prokka_data:
			sample_proteome = initial_sample_prokka_data[sample]['predicted_proteome']
			with open(sample_proteome) as osp:
				for rec in SeqIO.parse(osp, 'fasta'):
					if len(rec.id.split('_')[0]) != 3: continue
					original_samples.add(sample)
					if rec.id in gcf_protein_ids:
						all_gcf_proteins_fasta_handle.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
					else:
						all_comp_gcf_proteins_fasta_handle.write('>' + rec.id + '\n' + str(rec.seq) + '\n')

		all_gcf_proteins_fasta_handle.close()
		all_comp_gcf_proteins_fasta_handle.close()

		makedb_cmd = ['diamond', 'makedb', '--in', all_gcf_proteins_fasta_file, '-d', all_gcf_proteins_fasta_db]
		if logObject:
			logObject.info(
				'Running Diamond makedb on proteins from GCF with the following command: %s' % ' '.join(makedb_cmd))
		try:
			subprocess.call(' '.join(makedb_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if logObject:
				logObject.info('Successfully ran: %s' % ' '.join(makedb_cmd))
		except:
			if logObject:
				logObject.error('Had an issue running: %s' % ' '.join(makedb_cmd))
				logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(makedb_cmd))

		diamond_cmd = ['diamond', 'blastp', '--threads', str(cpus), '--ultra-sensitive', '--db',
					   all_gcf_proteins_fasta_db,
					   '--query', all_comp_gcf_proteins_fasta_file, '--outfmt', '6', '--out',
					   diamond_outfmt6_result_file,
					   '--max-target-seqs', '0']
		if logObject:
			logObject.info(
				'Running Diamond blastp between proteins not found in GCF against proteins found in GCF: %s' % ' '.join(
					diamond_cmd))
		try:
			subprocess.call(' '.join(diamond_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if logObject:
				logObject.info('Successfully ran: %s' % ' '.join(diamond_cmd))
		except:
			if logObject:
				logObject.error('Had an issue running: %s' % ' '.join(diamond_cmd))
				logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(diamond_cmd))

		hg_msa_pos_aligned = defaultdict(lambda: defaultdict(set))
		with open(diamond_outfmt6_result_file) as orf:
			for line in orf:
				line = line.strip()
				ls = line.split()
				qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = ls
				pident = float(pident)
				length = int(length)
				if (length >= 20 and pident >= 99.0) or (length >= 33 and pident >= 90.0):
					sstart = int(sstart)
					send = int(send)
					hg = gcf_protein_to_hg[sseqid]
					for pos in range(sstart, send + 1):
						msa_pos = gene_pos_to_msa_pos[hg][sseqid][pos]
						hg_msa_pos_aligned[hg][msa_pos].add(qseqid.split('_')[0])

		hg_differentiation_file = outdir + 'Non_Unique_Codon_Alignment_Positions.txt'
		hg_differentiation_handle = open(hg_differentiation_file, 'w')
		hg_nonunique_positions = defaultdict(set)
		for hg in hg_msa_pos_aligned:
			nonunique_positions = set([])
			for msa_pos in hg_msa_pos_aligned[hg]:
				if (len(hg_msa_pos_aligned[hg][msa_pos]) / float(len(original_samples))) >= 0.05:
					nonunique_positions.add(msa_pos)
			hg_nonunique_positions[hg] = nonunique_positions
			hg_differentiation_handle.write(
				'\t'.join([hg, ','.join([str(x) for x in sorted(nonunique_positions)])]) + '\n')
		hg_differentiation_handle.close()
		return hg_nonunique_positions
	except:
		if logObject:
			logObject.error("Issues with determining non-unique positions on profile HMMs.")
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())


def determineSeqSimProteinAlignment(protein_alignment_file, use_only_core=True):
	protein_sequences = {}
	with open(protein_alignment_file) as ocaf:
		for rec in SeqIO.parse(ocaf, 'fasta'):
			protein_sequences[rec.id] = str(rec.seq).upper()

	pair_seq_matching = defaultdict(lambda: defaultdict(lambda: 0.0))
	for i, g1 in enumerate(sorted(protein_sequences)):
		s1 = g1.split('|')[0]
		g1s = protein_sequences[g1]
		for j, g2 in enumerate(sorted(protein_sequences)):
			if i >= j: continue
			s2 = g2.split('|')[0]
			if s1 == s2: continue
			g2s = protein_sequences[g2]
			tot_comp_pos = 0
			match_pos = 0
			for pos, g1a in enumerate(g1s):
				g2a = g2s[pos]
				if g1a != '-' or g2a != '-':
					if not use_only_core or (use_only_core and g1a != '-' and g2a != '-'):
						tot_comp_pos += 1
						if g1a == g2a:
							match_pos += 1
			general_matching_percentage = 0.0
			if tot_comp_pos > 0:
				general_matching_percentage = float(match_pos) / float(tot_comp_pos)
			if pair_seq_matching[s1][s2] < general_matching_percentage and pair_seq_matching[s2][
				s1] < general_matching_percentage:
				pair_seq_matching[s1][s2] = general_matching_percentage
				pair_seq_matching[s2][s1] = general_matching_percentage

	return pair_seq_matching


def determineSeqSimCodonAlignment(codon_alignment_file, use_translation=False, use_only_core=True):
	gene_sequences = {}
	with open(codon_alignment_file) as ocaf:
		for i, rec in enumerate(SeqIO.parse(ocaf, 'fasta')):
			if use_translation:
				gene_sequences[rec.id] = str(rec.seq.upper().translate().upper())
			else:
				gene_sequences[rec.id] = str(rec.seq).upper()

	pair_seq_matching = defaultdict(lambda: defaultdict(lambda: 0.0))
	for i, g1 in enumerate(sorted(gene_sequences)):
		s1 = g1.split('|')[0]
		g1s = gene_sequences[g1]
		for j, g2 in enumerate(sorted(gene_sequences)):
			if i >= j: continue
			s2 = g2.split('|')[0]
			if s1 == s2: continue
			g2s = gene_sequences[g2]
			tot_comp_pos = 0
			match_pos = 0
			for pos, g1a in enumerate(g1s):
				g2a = g2s[pos]
				if g1a != '-' or g2a != '-':
					if not use_only_core or (use_only_core and g1a != '-' and g2a != '-'):
						tot_comp_pos += 1
						if g1a == g2a:
							match_pos += 1

			general_matching_percentage = 0.0
			if tot_comp_pos > 0:
				general_matching_percentage = float(match_pos) / float(tot_comp_pos)
			if pair_seq_matching[s1][s2] < general_matching_percentage and pair_seq_matching[s2][
				s1] < general_matching_percentage:
				pair_seq_matching[s1][s2] = general_matching_percentage
				pair_seq_matching[s2][s1] = general_matching_percentage

	return pair_seq_matching


def determineBGCSequenceSimilarity(input):
	hg, s1, g1s, s2, g2s, i, j, use_only_core, comparisons_managed = input

	comparison_id = '_|_'.join([hg, s1, s2, str(i), str(j)])

	tot_comp_pos = 0
	match_pos = 0
	for pos, g1a in enumerate(g1s):
		g2a = g2s[pos]
		if g1a != '-' or g2a != '-':
			if not use_only_core or (use_only_core and g1a != '-' and g2a != '-'):
				tot_comp_pos += 1
				if g1a == g2a:
					match_pos += 1
	general_matching_percentage = 0.0
	if tot_comp_pos > 0:
		general_matching_percentage = float(match_pos) / float(tot_comp_pos)
	comparisons_managed[comparison_id] = general_matching_percentage


def determineBGCSequenceSimilarityFromCodonAlignments(codon_alignments_file, cpus=1, use_translation=False,
													  use_only_core=True):
	sample_hgs = defaultdict(set)
	comparisons_managed = None
	pair_seq_matching = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0.0)))
	with multiprocessing.Manager() as manager:
		comparisons_managed = manager.dict()
		multiprocess_inputs = []
		with open(codon_alignments_file) as ocaf:
			for line in ocaf:
				line = line.strip()
				hg, codon_alignment = line.split('\t')
				gene_sequences = {}
				allele_identifiers = {}
				with open(codon_alignment) as oca:
					for i, rec in enumerate(SeqIO.parse(oca, 'fasta')):
						if use_translation:
							gene_sequences[rec.id] = str(rec.seq.upper().translate().upper())
						else:
							gene_sequences[rec.id] = str(rec.seq).upper()
						allele_identifiers[rec.id] = i
						sample = rec.id.split('|')[0]
						sample_hgs[sample].add(hg)

				for i, g1 in enumerate(gene_sequences):
					s1 = g1.split('|')[0]
					g1s = gene_sequences[g1]
					for j, g2 in enumerate(gene_sequences):
						if i >= j: continue
						s2 = g2.split('|')[0]
						if s1 == s2: continue
						g2s = gene_sequences[g2]
						multiprocess_inputs.append([hg, s1, g1s, s2, g2s, i, j, use_only_core, comparisons_managed])
		with manager.Pool(cpus) as pool:
			pool.map(determineBGCSequenceSimilarity, multiprocess_inputs)

		for comp in comparisons_managed:
			hg, s1, s2, i, j = comp.split('_|_')
			general_matching_percentage = comparisons_managed[comp]
			if pair_seq_matching[s1][s2][hg] < general_matching_percentage and pair_seq_matching[s2][s1][
				hg] < general_matching_percentage:
				pair_seq_matching[s1][s2][hg] = general_matching_percentage
				pair_seq_matching[s2][s1][hg] = general_matching_percentage

	bgc_pairwise_similarities = defaultdict(lambda: defaultdict(lambda: ["NA", "NA"]))
	for i, s1 in enumerate(sorted(sample_hgs)):
		for j, s2 in enumerate(sorted(sample_hgs)):
			if i >= j: continue
			common_hgs = sample_hgs[s1].intersection(sample_hgs[s2])
			total_hgs = sample_hgs[s1].union(sample_hgs[s2])
			sum_pair_seq_matching = 0.0
			for hg in common_hgs:
				sum_pair_seq_matching += pair_seq_matching[s1][s2][hg]
			if len(common_hgs) > 0:
				bgc_pairwise_similarities[s1][s2] = [sum_pair_seq_matching / float(len(common_hgs)),
													 float(len(common_hgs)) / float(len(total_hgs))]
				bgc_pairwise_similarities[s2][s1] = [sum_pair_seq_matching / float(len(common_hgs)),
													 float(len(common_hgs)) / float(len(total_hgs))]
			else:
				bgc_pairwise_similarities[s1][s2] = ["NA", 0.0]
				bgc_pairwise_similarities[s2][s1] = ["NA", 0.0]

	return bgc_pairwise_similarities


def determineAllelesFromCodonAlignment(codon_alignment, max_mismatch=10, matching_percentage_cutoff=0.99,
									   filter_by_genelength=True):
	gene_sequences = {}
	gene_sequences_lengths = []
	allele_identifiers = {}
	seqs_comprehensive = set([])
	with open(codon_alignment) as oca:
		for i, rec in enumerate(SeqIO.parse(oca, 'fasta')):
			gene_sequences_lengths.append(len(str(rec.seq).upper().replace('N', '').replace('-', '')))
	median_length = statistics.median(gene_sequences_lengths)
	mad_length = max(stats.median_abs_deviation(gene_sequences_lengths, scale="normal"), 5)
	with open(codon_alignment) as oca:
		for i, rec in enumerate(SeqIO.parse(oca, 'fasta')):
			gene_seq_len = len(str(rec.seq).upper().replace('N', '').replace('-', ''))
			if filter_by_genelength and (
					gene_seq_len < (median_length - mad_length) or gene_seq_len > (median_length + mad_length)):
				continue
			gene_sequences[rec.id] = str(rec.seq).upper()
			allele_identifiers[rec.id] = i
			seqs_comprehensive.add(rec.id)

	pairs = []
	seqs_paired = set([])
	pair_matching = defaultdict(lambda: defaultdict(float))
	for i, g1 in enumerate(gene_sequences):
		g1s = gene_sequences[g1]
		for j, g2 in enumerate(gene_sequences):
			if i >= j: continue
			g2s = gene_sequences[g2]
			tot_comp_pos = 0
			g1_comp_pos = 0
			g2_comp_pos = 0
			match_pos = 0
			mismatch_pos = 0
			for pos, g1a in enumerate(g1s):
				g2a = g2s[pos]
				if g1a in valid_alleles and g2a in valid_alleles:
					if g1a != g2a:
						mismatch_pos += 1
				if g1a in valid_alleles or g2a in valid_alleles:
					tot_comp_pos += 1
					if g1a == g2a:
						match_pos += 1
				if g1a in valid_alleles:
					g1_comp_pos += 1
				if g2a in valid_alleles:
					g2_comp_pos += 1
			general_matching_percentage = float(match_pos) / float(tot_comp_pos)
			g1_matching_percentage = float(match_pos) / float(g1_comp_pos)
			g2_matching_percentage = float(match_pos) / float(g2_comp_pos)
			pair_matching[g1][g2] = general_matching_percentage
			if general_matching_percentage >= matching_percentage_cutoff or g1_matching_percentage >= matching_percentage_cutoff or g2_matching_percentage >= matching_percentage_cutoff:
				if mismatch_pos <= max_mismatch:
					seqs_paired.add(g1)
					seqs_paired.add(g2)
					pairs.append(sorted([g1, g2]))

	"""	
	Solution for single-linkage clustering taken from mimomu's repsonse in the stackoverflow page:
	https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements?lq=1
	"""
	L = pairs
	LL = set(itertools.chain.from_iterable(L))
	for each in LL:
		components = [x for x in L if each in x]
		for i in components:
			L.remove(i)
		L += [list(set(itertools.chain.from_iterable(components)))]

	for seq in seqs_comprehensive:
		if not seq in seqs_paired:
			L.append([seq])

	allele_cluster_min_id = {}
	for allele_cluster in L:
		gene_identifiers = set([])
		for gene in allele_cluster:
			gene_identifiers.add(allele_identifiers[gene])
		min_gi = min(gene_identifiers)
		allele_cluster_min_id[min_gi] = allele_cluster

	allele_clusters = defaultdict(set)
	for i, aci in enumerate(sorted(allele_cluster_min_id.keys())):
		for gene in allele_cluster_min_id[aci]:
			allele_clusters['Allele_Cluster_' + str(i + 1)].add(gene)

	return [allele_clusters, pair_matching]


def cleanUpSampleName(original_name):
	return original_name.replace('#', '').replace('*', '_').replace(':', '_').replace(';', '_').replace(' ',
																										'_').replace(
		':', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_").replace('-', '_').replace('(',
																													'').replace(
		')', '').replace('/', '').replace('\\', '').replace('[', '').replace(']', '').replace(',', '')


def read_pair_generator_defunct(bam, region_string=None, start=None, stop=None):
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


def read_pair_generator(bam, region_string, reference_length, max_insert_size=500):
	"""
    Function adapted from: https://www.biostars.org/p/306041/
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """

	read_reverse_counts = defaultdict(int)
	read_forward_counts = defaultdict(int)
	for read in bam.fetch(region_string):
		if read.is_supplementary: continue
		qname = read.query_name
		if read.is_reverse:
			read_reverse_counts[qname] += 1
		else:
			read_forward_counts[qname] += 1

	visited = set([])
	read_dict = defaultdict(lambda: [None, None])
	for read in bam.fetch(region_string):
		if read.is_supplementary: continue
		qname = read.query_name
		if read_reverse_counts[qname] > 1 or read_forward_counts[qname] > 1: continue
		if qname not in read_dict:
			if read.is_read1:
				read_dict[qname][0] = read
			else:
				read_dict[qname][1] = read
		else:
			# report proper pair paired-end reads, discard cases where both reads map to reference contig
			# but in an improper fashion.
			if read.is_proper_pair:
				if read.is_read1:
					yield read, read_dict[qname][1]
				else:
					yield read_dict[qname][0], read
			visited.add(qname)
			del read_dict[qname]

	for read in bam.fetch(region_string):
		qname = read.query_name
		if qname in visited: continue
		if read.is_supplementary or (not read.mate_is_unmapped): continue

		first_real_alignment_pos = None
		last_real_alignment_pos = None
		indel_positions = set([])
		match = 0
		aligned = 0
		for b in read.get_aligned_pairs(with_seq=True):
			if b[0] != None and b[1] != None and b[2] != None:
				if first_real_alignment_pos == None:
					first_real_alignment_pos = b[0]
				last_real_alignment_pos = b[0]
				if b[2].isupper():
					match += 1
				aligned += 1
			else:
				indel_positions.add(b[0])
		main_alignment_positions = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
		has_indel = len(main_alignment_positions.intersection(indel_positions)) > 0

		matching_percentage = match / aligned
		positionally_checks_out = False
		if read.is_reverse:
			min_position = min(main_alignment_positions)
			positionally_checks_out = (min_position - max_insert_size) < 0
		else:
			max_position = max(main_alignment_positions)
			positionally_checks_out = (max_position + max_insert_size) > reference_length

		if (positionally_checks_out):  # and (matching_percentage >= 0.95):
			yield read_dict[qname][0], read_dict[qname][1]


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


def runBowtie2Alignments(bowtie2_reference, paired_end_sequencing_file, bowtie2_outdir, logObject, cpus=1):
	"""
	Wrapper function for running Bowtie2 alignments to reference database/index

	:param bowtie2_reference: path to Bowtie2 reference/index
	:param paired_end_sequencing_file: tab delimited file with three columns: (1) sample name (2) path to forward
									   reads and (3) path to reverse reads
	:param bowtie2_outdir: Path to directory where Bowtie2 results should be written
	:param logObject: logging object for documentation
	:param cpus: number of cpus (total) to use. If more than 4 cpus provided, then parallel Bowtie2 jobs with 4 cpus
				  each will be started.
	"""
	bowtie2_cpus = cpus
	bowtie2_pool_size = 1
	if cpus >= 4:
		bowtie2_cpus = 4
		bowtie2_pool_size = int(cpus / 4)

	try:
		bowtie2_inputs = []
		with open(paired_end_sequencing_file) as opesf:
			for line in opesf:
				line = line.strip()
				sample = line.split('\t')[0]
				reads = line.split('\t')[1:]
				bowtie2_inputs.append([sample, reads, bowtie2_reference, bowtie2_outdir, bowtie2_cpus, logObject])
		p = multiprocessing.Pool(bowtie2_pool_size)
		p.map(bowtie2_alignment, bowtie2_inputs)
		p.close()
	except Exception as e:
		logObject.error("Issues in setting up and running Bowtie2 alignments.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())


def bowtie2_alignment(input_args):
	"""
	Function to perform Bowtie2 alignment of paired-end reads to a database/reference and post-processing of alignment
	file with samtools (e.g. convert to BAM format, sort BAM file, and index it).
	"""
	sample, reads, bowtie2_reference, bowtie2_outdir, bowtie2_cpus, logObject = input_args

	sam_file = bowtie2_outdir + sample + '.sam'
	bam_file = bowtie2_outdir + sample + '.bam'
	bam_file_sorted = bowtie2_outdir + sample + '.sorted.bam'
	# bam_file_filtered = bowtie2_outdir + sample + '.filtered.bam'
	# am_file_filtered_sorted = bowtie2_outdir + sample + '.filtered.sorted.bam'

	bowtie2_cmd = ['bowtie2', '--very-sensitive-local', '--no-unal', '-a', '-x', bowtie2_reference, '-U',
				   ','.join(reads), '-S', sam_file, '-p', str(bowtie2_cpus)]

	samtools_view_cmd = ['samtools', 'view', '-h', '-Sb', sam_file, '>', bam_file]
	samtools_sort_cmd = ['samtools', 'sort', '-@', str(bowtie2_cpus), bam_file, '-o', bam_file_sorted]
	samtools_index_cmd = ['samtools', 'index', bam_file_sorted]

	try:
		run_cmd(bowtie2_cmd, logObject)
		run_cmd(samtools_view_cmd, logObject)
		run_cmd(samtools_sort_cmd, logObject)
		run_cmd(samtools_index_cmd, logObject)

		"""
		bam_handle = pysam.AlignmentFile(bam_file_sorted, 'rb')
		filt_bam_handle = pysam.AlignmentFile(bam_file_filtered, "wb", template=bam_handle)

		for read in bam_handle.fetch():
			if not read.is_proper_pair or read.is_supplementary: continue
			filt_bam_handle.write(read)

		filt_bam_handle.close()
		bam_handle.close()

		run_cmd(samtools_sort_cmd_2, logObject)
		run_cmd(samtools_index_cmd_2, logObject)
		"""

		os.system(
			"rm -f %s %s" % (sam_file, bam_file))  # bam_file_sorted, bam_file_filtered, bam_file_sorted + '.bai'))
	except Exception as e:
		if bowtie2_outdir != "" and sample != "":
			os.system('rm -f %s/%s*' % (bowtie2_outdir, sample))
		raise RuntimeError(traceback.format_exc())


def createBGCGenbank(full_genbank_file, new_genbank_file, scaffold, start_coord, end_coord):
	"""
	Function to prune full genome-sized Genbank for only features in BGC of interest.

	:param full_genbank_file: Prokka generated Genbank file for full genome.
	:param new_genbank_file: Path to BGC specific Genbank to be created
	:param scaffold: Scaffold identifier.
	:param start_coord: Start coordinate.
	:param end_coord: End coordinate.
	"""
	try:
		ngf_handle = open(new_genbank_file, 'w')
		pruned_coords = set(range(start_coord, end_coord + 1))
		with open(full_genbank_file) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
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
					start = None
					end = None
					direction = None
					all_coords = []

					if not 'join' in str(feature.location) and not 'order' in str(feature.location):
						start = min([int(x.strip('>').strip('<')) for x in
									 str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max(
							[int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
						direction = str(feature.location).split('(')[1].split(')')[0]
						all_coords.append([start, end, direction])
					elif 'order' in str(feature.location):
						all_starts = []
						all_ends = []
						all_directions = []
						for exon_coord in str(feature.location)[6:-1].split(', '):
							start = min(
								[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
							end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
							direction = exon_coord.split('(')[1].split(')')[0]
							all_starts.append(start)
							all_ends.append(end)
							all_directions.append(direction)
							all_coords.append([start, end, direction])
						start = min(all_starts)
						end = max(all_ends)
					else:
						all_starts = []
						all_ends = []
						all_directions = []
						for exon_coord in str(feature.location)[5:-1].split(', '):
							start = min(
								[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
							end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
							direction = exon_coord.split('(')[1].split(')')[0]
							all_starts.append(start)
							all_ends.append(end)
							all_directions.append(direction)
							all_coords.append([start, end, direction])
						start = min(all_starts)
						end = max(all_ends)

					feature_coords = set(range(start, end + 1))
					if len(feature_coords.intersection(pruned_coords)) > 0:
						fls = []
						for sc, ec, dc in all_coords:
							updated_start = sc - start_coord + 1
							updated_end = ec - start_coord + 1
							if ec > end_coord:
								# note overlapping genes in prokaryotes are possible so avoid proteins that overlap
								# with boundary proteins found by the HMM.
								if feature.type == 'CDS':
									continue
								else:
									updated_end = end_coord - start_coord + 1  # ; flag1 = True
							if sc < start_coord:
								if feature.type == 'CDS':
									continue
								else:
									updated_start = 1  # ; flag2 = True

							strand = 1
							if dc == '-':
								strand = -1
							fls.append(FeatureLocation(updated_start - 1, updated_end, strand=strand))
						if len(fls) > 0:
							updated_location = fls[0]
							if len(fls) > 1:
								updated_location = sum(fls)
							feature.location = updated_location
							updated_features.append(feature)
				updated_rec.features = updated_features
				SeqIO.write(updated_rec, ngf_handle, 'genbank')
		ngf_handle.close()
	except Exception as e:
		raise RuntimeError(traceback.format_exc())


def parseGenbankAndFindBoundaryGenes(inputs):
	"""
	Function to parse Genbanks from Prokka and return a dictionary of genes per scaffold, gene to scaffold, and a
	set of genes which lie on the boundary of scaffolds.

	:param sample_genbank: Prokka generated Genbank file.
	:param distance_to_scaffold_boundary: Distance to scaffold edge considered as boundary.

	:return gene_to_scaffold: Dictionary mapping each gene's locus tag to the scaffold it is found on.
	:return scaffold_genes: Dictionary with keys as scaffolds and values as a set of genes found on that scaffold.
	:return boundary_genes: Set of gene locus tag ids which are found within proximity to scaffold edges.
	"""

	distance_to_scaffold_boundary = 2000
	gene_location = {}
	scaffold_genes = defaultdict(set)
	boundary_genes = set([])
	gene_id_to_order = defaultdict(dict)
	gene_order_to_id = defaultdict(dict)

	sample, sample_genbank, sample_gbk_info = inputs
	osg = None
	if sample_genbank.endswith('.gz'):
		osg = gzip.open(sample_genbank, 'rt')
	else:
		osg = open(sample_genbank)
	for rec in SeqIO.parse(osg, 'genbank'):
		scaffold = rec.id
		scaffold_length = len(str(rec.seq))
		boundary_ranges = set(range(1, distance_to_scaffold_boundary + 1)).union(
			set(range(scaffold_length - distance_to_scaffold_boundary, scaffold_length + 1)))
		gene_starts = []
		for feature in rec.features:
			if not feature.type == 'CDS': continue
			locus_tag = feature.qualifiers.get('locus_tag')[0]

			start = None
			end = None
			direction = None
			if not 'join' in str(feature.location):
				start = min(
					[int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
				end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
				direction = str(feature.location).split('(')[1].split(')')[0]
			else:
				all_starts = []
				all_ends = []
				all_directions = []
				for exon_coord in str(feature.location)[5:-1].split(', '):
					start = min([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
					end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
					direction = exon_coord.split('(')[1].split(')')[0]
					all_starts.append(start)
					all_ends.append(end)
					all_directions.append(direction)
				start = min(all_starts)
				end = max(all_ends)
				direction = all_directions[0]

			gene_location[locus_tag] = {'scaffold': scaffold, 'start': start, 'end': end, 'direction': direction}
			scaffold_genes[scaffold].add(locus_tag)

			gene_range = set(range(start, end + 1))
			if len(gene_range.intersection(boundary_ranges)) > 0:
				boundary_genes.add(locus_tag)

			gene_starts.append([locus_tag, start])

		for i, g in enumerate(sorted(gene_starts, key=itemgetter(1))):
			gene_id_to_order[scaffold][g[0]] = i
			gene_order_to_id[scaffold][i] = g[0]
	osg.close()
	sample_gbk_info[sample] = [gene_location, dict(scaffold_genes), boundary_genes, dict(gene_id_to_order),
							   dict(gene_order_to_id)]


def chunks(lst, n):
	"""
    Yield successive n-sized chunks from lst.
    Solution taken from: https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """
	for i in range(0, len(lst), n):
		yield lst[i:i + n]


def calculateMashPairwiseDifferences(fasta_listing_file, outdir, name, sketch_size, cpus, logObject, prune_set=None):
	"""
	Calculate MASH pairwise distances (estimated ANI) between FASTA files.

	:param fasta_listing_file: A tab-delimited listing file with two columns: (1) sample name (2) path to FASTA file
	:param outdir: The output directory where to write results
	:param name: Name of analysis scope
	:param sketch_size: Sketch size (a parameter of MASH)
	:param cpus: Number of cpus/threads to use
	:param logObject: The logging object.
	"""
	mash_db = outdir + name
	fastas = []
	fasta_to_name = {}
	try:
		with open(fasta_listing_file) as oflf:
			for line in oflf:
				line = line.strip()
				ls = line.split('\t')
				if prune_set != None and not ls[0] in prune_set: continue
				fastas.append(ls[1])
				fasta_to_name[ls[1]] = ls[0]
	except:
		error_message = "Had issues reading the FASTA listing file %s" % fasta_listing_file
		logObject.error(error_message)
		raise RuntimeError(error_message)

	mash_input_file = outdir + 'MASH_Input.txt'
	mash_input_handle = open(mash_input_file, 'w')
	mash_input_handle.write('\n'.join(fastas))
	mash_input_handle.close()

	# create mash database (using mash sketch)
	mash_sketch_cmd = ['mash', 'sketch', '-p', str(cpus), '-s', str(sketch_size), '-o', mash_db, '-l', mash_input_file]
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
	mash_dist_cmd = ['mash', 'dist', '-s', str(sketch_size), '-p', str(cpus), mash_db, mash_db, '>',
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

	pairwise_similarities = defaultdict(lambda: defaultdict(float))
	try:
		with open(outdir + name + '.out') as of:
			for line in of:
				line = line.strip()
				ls = line.split('\t')
				f1, f2, dist = ls[:3]
				dist = float(dist)
				n1 = fasta_to_name[f1]
				n2 = fasta_to_name[f2]
				pairwise_similarities[n1][n2] = (1.0 - dist)
	except:
		error_message = 'Had issues reading the output of MASH dist anlaysis in: %s.out' % outdir + name
		logObject.error(error_message)
		raise RuntimeError(error_message)
	return pairwise_similarities


def runFastANI(fasta_listing_file, outdir, fastani_output_file, cpus, logObject, prune_set=None):
	"""
	Calculate ANI estimate between pairs of samples using FastANI.

	:param fasta_listing_file: A tab-delimited listing file with two columns: (1) sample name (2) path to FASTA file
	:param outdir: The output directory where to write results
	:param fastani_output_file: Output file name
	:param cpus: Number of cpus/threads to use
	:param logObject: The logging object.
	"""
	fastas = []
	fasta_to_name = {}
	try:
		with open(fasta_listing_file) as oflf:
			for line in oflf:
				line = line.strip()
				ls = line.split('\t')
				if prune_set != None and not ls[0] in prune_set: continue
				fastas.append(ls[1])
				fasta_to_name[ls[1]] = ls[0]
	except:
		error_message = "Had issues reading the FASTA listing file %s" % fasta_listing_file
		logObject.error(error_message)
		raise RuntimeError(error_message)

	fastani_input_file = outdir + 'FastANI_Input.txt'
	fastani_input_handle = open(fastani_input_file, 'w')
	fastani_input_handle.write('\n'.join(fastas))
	fastani_input_handle.close()

	if not os.path.isfile(fastani_output_file):
		fastani_cmd = ['fastANI', '-t', str(cpus), '--ql', fastani_input_file, '--rl', fastani_input_file, '-o',
					   fastani_output_file]
		logObject.info('Running fastANI  with the following command: %s' % ' '.join(fastani_cmd))
		try:
			subprocess.call(' '.join(fastani_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran: %s' % ' '.join(fastani_cmd))
		except:
			error_message = 'Had an issue running: %s' % ' '.join(fastani_cmd)
			logObject.error(error_message)
			raise RuntimeError(error_message)
		try:
			assert (os.path.isfile(fastani_output_file))
		except:
			error_message = "Had issue validating that FastANI ran properly, couldn't find: %s" % fastani_output_file
			logObject.error(error_message)
			raise RuntimeError(error_message)
	else:
		logObject.info('fastANI result already exists, assuming valid and avoiding rerunning. ')

	pairwise_similarities = defaultdict(lambda: defaultdict(float))
	pairwise_comparisons = defaultdict(lambda: defaultdict(float))
	try:
		with open(fastani_output_file) as of:
			for line in of:
				line = line.strip()
				ls = line.split('\t')
				f1, f2, sim, comp_seg, tota_seg = ls
				sim = float(sim)
				comp_prop = float(comp_seg) / float(tota_seg)
				n1 = fasta_to_name[f1]
				n2 = fasta_to_name[f2]
				pairwise_similarities[n1][n2] = sim / 100.0
				pairwise_comparisons[n1][n2] = comp_prop
	except:
		error_message = 'Had issues reading the output of FastANI analysis at: %s' % fastani_output_file
		logObject.error(error_message)
		raise RuntimeError(error_message)
	return [pairwise_similarities, pairwise_comparisons]


def runCompareM(fasta_listing_file, comparem_results_dir, cpus, logObject, prune_set=None):
	"""
	Calculate AAI estimate between pairs of samples using CompareM.

	:param fasta_listing_file: A tab-delimited listing file with two columns: (1) sample name (2) path to FASTA file
	:param comparem_results_dir: The output directory where to write results
	:param cpus: Number of cpus/threads to use
	:param logObject: The logging object.
	"""
	fastas = []
	fasta_to_name = {}
	try:
		with open(fasta_listing_file) as oflf:
			for line in oflf:
				line = line.strip()
				ls = line.split('\t')
				if prune_set != None and not ls[0] in prune_set: continue
				fastas.append(ls[1])
				fasta_to_name[ls[1]] = ls[0]
	except:
		error_message = "Had issues reading the FASTA listing file %s" % fasta_listing_file
		logObject.error(error_message)
		raise RuntimeError(error_message)

	comparem_input_file = comparem_results_dir + 'CompareM_Input.txt'
	comparem_input_handle = open(comparem_input_file, 'w')
	comparem_input_handle.write('\n'.join(fastas))
	comparem_input_handle.close()

	tmp_dir = comparem_results_dir + 'tmp/'
	if not os.path.isdir(tmp_dir): os.system('mkdir %s' % tmp_dir)

	comparem_result_file = comparem_results_dir + 'aai/aai_summary.tsv'
	if not os.path.isfile(comparem_result_file):
		comparem_cmd = ['comparem', 'aai_wf', '--tmp_dir', tmp_dir, '--cpus', str(cpus), comparem_input_file,
						comparem_results_dir]
		logObject.info('Running CompareM  with the following command: %s' % ' '.join(comparem_cmd))
		try:
			subprocess.call(' '.join(comparem_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran: %s' % ' '.join(comparem_cmd))
		except:
			error_message = 'Had an issue running: %s' % ' '.join(comparem_cmd)
			logObject.error(error_message)
			raise RuntimeError(error_message)
		try:
			assert (os.path.isfile(comparem_result_file))
		except:
			error_message = "Had issue validating that CompareM ran properly, couldn't find: %s" % comparem_result_file
			logObject.error(error_message)
			raise RuntimeError(error_message)
	else:
		logObject.info('CompareM result already exists, assuming valid and avoiding rerunning.')

	pairwise_similarities = defaultdict(lambda: defaultdict(float))
	pairwise_comparisons = defaultdict(lambda: defaultdict(float))
	try:
		with open(comparem_result_file) as of:
			for i, line in enumerate(of):
				if i == 0: continue
				line = line.strip()
				s1, s1g, s2, s2g, com_genes, sim, sim_std, ortho_frac = line.split('\t')
				sim = float(sim)
				pairwise_similarities[s1][s2] = sim / 100.0
				pairwise_similarities[s2][s1] = sim / 100.0
				pairwise_comparisons[s1][s2] = float(com_genes) / float(s1g)
				pairwise_comparisons[s2][s1] = float(com_genes) / float(s2g)
	except:
		error_message = 'Had issues reading the output of CompareM analysis at: %s' % comparem_result_file
		logObject.error(error_message)
		raise RuntimeError(error_message)
	return [pairwise_similarities, pairwise_comparisons]


def parseOrthoFinderMatrix(orthofinder_matrix_file, relevant_gene_lts, all_primary=False):
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
					# critical for calculating homolog group stats, like median gene counts, multicopy proportion
					# use only genes from the original set of genomes used to conduct full orthofinder analysis.
					gene_counts.append(len([x for x in sgs.split(', ') if (len(x.split('_')[0]) == 3 or all_primary)]))

				hg_multicopy_proportion[hg] = float(sum([1 for x in gene_counts if x > 1])) / sum(
					[1 for x in gene_counts if x > 0])
				hg_median_gene_counts[hg] = statistics.median(gene_counts)

	return ([gene_to_hg, hg_genes, hg_median_gene_counts, hg_multicopy_proportion])


def run_cmd(cmd, logObject, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL):
	"""
	Simple function to run a single command through subprocess with logging.
	"""
	logObject.info('Running the following command: %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=stderr, executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except Exception as e:
		logObject.error('Had an issue running: %s' % ' '.join(cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))


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
	except Exception as e:
		logObject.warning('Had an issue running: %s' % ' '.join(input_cmd))
		logObject.warning(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())


def processGenomes(sample_genomes, prodigal_outdir, prodigal_proteomes, prodigal_genbanks, logObject, cpus=1,
				   locus_tag_length=3, use_pyrodigal=False, avoid_locus_tags=set([])):
	"""
	Void function to run Prodigal based gene-calling and annotations.

	:param sample_genomes: dictionary with keys as sample names and values as genomic assembly paths.
	:param prodigal_outdir: full path to directory where Prokka results will be written.
	:param prodigal_proteomes: full path to directory where Prokka generated predicted-proteome FASTA files will be moved after prodigal has run.
	:param prodigal_genbanks: full path to directory where Prokka generated Genbank (featuring predicted CDS) files will be moved after prodigal has run.
	:param taxa: name of the taxonomic clade of interest.
	:param logObject: python logging object handler.
	:param cpus: number of cpus to use in multiprocessing Prokka cmds.
	:param locus_tag_length: length of locus tags to generate using unique character combinations.

	Note length of locus tag must be 3 beause this is substituting for base lsaBGC analysis!!
	"""

	prodigal_cmds = []
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list(
			set([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))]).difference(
				avoid_locus_tags)))
		for i, sample in enumerate(sample_genomes):
			sample_assembly = sample_genomes[sample]
			sample_locus_tag = ''.join(list(possible_locustags[i]))

			prodigal_cmd = ['runProdigalAndMakeProperGenbank.py', '-i', sample_assembly, '-s', sample,
							'-l', sample_locus_tag, '-o', prodigal_outdir]
			if use_pyrodigal:
				prodigal_cmd += ['-py']
			prodigal_cmds.append(prodigal_cmd + [logObject])

		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, prodigal_cmds)
		p.close()

		for sample in sample_genomes:
			try:
				assert (os.path.isfile(prodigal_outdir + sample + '.faa') and os.path.isfile(
					prodigal_outdir + sample + '.gbk'))
				os.system('mv %s %s' % (prodigal_outdir + sample + '.gbk', prodigal_genbanks))
				os.system('mv %s %s' % (prodigal_outdir + sample + '.faa', prodigal_proteomes))
			except:
				raise RuntimeError(
					"Unable to validate successful genbank/predicted-proteome creation for sample %s" % sample)
	except Exception as e:
		logObject.error(
			"Problem with creating commands for running prodigal via script runProdigalAndMakeProperGenbank.py. Exiting now ...")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())


def parseSampleGenomes(genome_listing_file, logObject):
	try:
		sample_genomes = {}
		all_genbanks = True
		all_fastas = True
		at_least_one_genbank = False
		at_least_one_fasta = False
		with open(genome_listing_file) as oglf:
			for line in oglf:
				line = line.strip()
				ls = line.split('\t')
				sample, genome_file = ls
				try:
					assert (os.path.isfile(genome_file))
				except:
					logObject.warning(
						"Problem with finding genome file %s for sample %s, skipping" % (genome_file, sample))
					continue
				if sample in sample_genomes:
					logObject.warning(
						'Skipping genome %s for sample %s because a genome file was already provided for this sample' % (
						genome_file, sample))
					continue

				sample_genomes[sample] = genome_file
				if not is_fasta(genome_file):
					all_fastas = False
				else:
					at_least_one_fasta = True
				if not is_genbank(genome_file):
					all_genbanks = False
				else:
					at_least_one_genbank = True

		format_prediction = 'mixed'
		if all_genbanks and at_least_one_genbank:
			format_prediction = 'genbank'
		elif all_fastas and at_least_one_fasta:
			format_prediction = 'fasta'

		return ([sample_genomes, format_prediction])

	except Exception as e:
		logObject.error("Problem with creating commands for running Prodigal. Exiting now ...")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())


def processGenomesAsGenbanks(sample_genomes, proteomes_directory, genbanks_directory, gene_name_mapping_outdir,
							 logObject, cpus=1, locus_tag_length=3, avoid_locus_tags=set([])):
	"""
	Extracts CDS/proteins from existing Genbank files and recreates
	"""

	sample_genomes_updated = {}
	process_cmds = []
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list(
			set([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))]).difference(
				avoid_locus_tags)))
		lacking_cds_gbks = set([])

		for i, sample in enumerate(sample_genomes):
			sample_genbank = sample_genomes[sample]
			ogh = None
			if sample_genbank.endswith('.gz'):
				ogh = gzip.open(sample_genbank, 'rt')
			else:
				ogh = open(sample_genbank)
			cds_flag = False
			for rec in SeqIO.parse(ogh, 'genbank'):
				for feature in rec.features:
					if feature.type == 'CDS':
						cds_flag = True
						break
			ogh.close()
			sample_locus_tag = ''.join(list(possible_locustags[i]))
			if not cds_flag:
				lacking_cds_gbks.add(sample)
				logObject.warning('NCBI genbank file %s for sample %s lacks CDS features' % (sample_genbank, sample))
				continue
			process_cmd = ['processAndReformatNCBIGenbanks.py', '-i', sample_genbank, '-s', sample,
						   '-l', sample_locus_tag, '-g', genbanks_directory, '-p', proteomes_directory, '-n',
						   gene_name_mapping_outdir]
			process_cmds.append(process_cmd + [logObject])

		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, process_cmds)
		p.close()

		for sample in sample_genomes:
			if sample in lacking_cds_gbks:
				continue
			try:
				assert (os.path.isfile(proteomes_directory + sample + '.faa') and
						os.path.isfile(genbanks_directory + sample + '.gbk') and
						os.path.isfile(gene_name_mapping_outdir + sample + '.txt'))
				sample_genomes_updated[sample] = genbanks_directory + sample + '.gbk'
			except:
				raise RuntimeError(
					"Unable to validate successful genbank/predicted-proteome creation for sample %s" % sample)
	except Exception as e:
		logObject.error("Problem with processing existing Genbanks to (re)create genbanks/proteomes. Exiting now ...")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return sample_genomes


def splitDeepBGCGenbank(bgc_genbank_listing_file, deepbgc_split_directory, outdir, logObject):
	updated_bgc_genbank_listing_file = outdir + 'DeepBGC_Genbanks_Split_Listing.txt'
	updated_bgc_genbank_listing_handle = open(updated_bgc_genbank_listing_file, 'w')
	try:
		with open(bgc_genbank_listing_file) as obglf:
			for line in obglf:
				sample, bgc = line.strip().split('\t')
				logObject.info('Splitting DeepBGC Genbank %s for sample %s.' % (bgc, sample))
				with open(bgc) as obgbk:
					for rec in SeqIO.parse(obgbk, 'genbank'):
						for feature in rec.features:
							if feature.type == "CDS":
								start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
								end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
								direction = str(feature.location).split('(')[1].split(')')[0]
								if end:
									nucl_seq = str(rec.seq)[start - 1:end]
								else:
									nucl_seq = str(rec.seq)[start - 1:]
								if direction == '-':
									nucl_seq = str(Seq(nucl_seq).reverse_complement())
								feature.qualifiers['translation'] = Seq(nucl_seq).translate()
						out_gbk = deepbgc_split_directory + rec.id.replace('.', '_') + '.gbk'
						updated_bgc_genbank_listing_handle.write(sample + '\t' + out_gbk + '\n')
						out_gbk_handle = open(out_gbk, 'w')
						SeqIO.write(rec, out_gbk_handle, 'genbank')
						out_gbk_handle.close()
		updated_bgc_genbank_listing_handle.close()
	except Exception as e:
		logObject.error("Problem with parsing DeepBGC Genbanks to split.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return updated_bgc_genbank_listing_file


def findBGCInFullGenbank(inputs):
	try:
		bgc_gbk, full_gbk, outf = inputs

		bgc_seq = None
		with open(bgc_gbk) as obgf:
			for rec in SeqIO.parse(obgf, 'genbank'):
				bgc_seq = str(rec.seq.upper())

		scaff_id, scaff_start = [None] * 2
		off = None
		if full_gbk.endswith('.gz'):
			off = gzip.open(full_gbk, 'rt')
		else:
			off = open(full_gbk)
		outf_handle = open(outf, 'w')
		for rec in SeqIO.parse(off, 'genbank'):
			if bgc_seq in str(rec.seq).upper():
				scaff_id = rec.id
				scaff_start = str(rec.seq).find(bgc_seq)
				outf_handle.write(scaff_id + '\t' + str(scaff_start) + '\n')
		outf_handle.close()
		if off:
			off.close()
	except Exception as e:
		raise RuntimeError(traceback.format_exc())


def mapBGCtoGenomeBySequence(bgc_genbank_listing_file, sample_genomes, outdir, logObject, cpus=1):
	bgc_mappings = {}
	locate_bgc_directory = outdir + 'Map_BGC_to_Full_Genbanks/'
	bgcs_without_mappings_handle = open(outdir + 'BGCs_Unable_to_be_Mapped_to_Genome.txt', 'w')
	bgcs_with_multiple_mappings_handle = open(outdir + 'BGCs_with_Multiple_Perfect_Mappings_to_Genome.txt', 'w')
	setupReadyDirectory([locate_bgc_directory])
	try:
		outf_to_info = {}
		find_inputs = []
		with open(bgc_genbank_listing_file) as obglf:
			for line in obglf:
				line = line.strip()
				sample, bgc_gbk = line.split('\t')
				full_gbk = sample_genomes[sample]
				outf = locate_bgc_directory + sample + '_BGC-ID_' + bgc_gbk.split('/')[-1] + '.txt'
				outf_to_info[outf] = [bgc_gbk, sample]
				find_inputs.append([bgc_gbk, full_gbk, outf])

		p = multiprocessing.Pool(cpus)
		p.map(findBGCInFullGenbank, find_inputs)
		p.close()

		loc_tuples = set([])
		for f in os.listdir(locate_bgc_directory):
			outf = locate_bgc_directory + f
			bgc_gbk, sample = outf_to_info[outf]
			count = 0
			with open(outf) as of:
				for i, line in enumerate(of):
					scaff, start = line.strip().split('\t')
					loc_tuple = tuple([sample, scaff, start])
					count += 1
					if loc_tuple in loc_tuples: continue
					bgc_mappings[bgc_gbk] = [scaff, int(start)]
					loc_tuples.add(loc_tuple)
			if count == 0:
				logObject.warning('Dropping BGC %s for sample %s because the BGC sequence did not match any scaffold directly which should not be the case!' % (bgc_gbk, sample))
				bgcs_without_mappings_handle.write(sample + '\t' + bgc_gbk + '\n')
			elif count > 1:
				bgcs_with_multiple_mappings_handle.write(
					sample + '\t' + bgc_gbk + '\t' + str(count) + '\t' + str(bgc_mappings[bgc_gbk]) + '\n')

	except Exception as e:
		logObject.error("Issues with parsing out protein sequences from BGC Genbanks.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	bgcs_with_multiple_mappings_handle.close()
	bgcs_without_mappings_handle.close()
	return bgc_mappings


def processBGCGenbanks(bgc_listing_file, bgc_mappings, bgc_prediction_software, sample_genomes,
					   antismash_bgcs_directory, proteomes_directory, logObject):
	"""
	Creates local versions of BGC prediction Genbanks where proteins have locus tags updated.
	"""

	sample_bgcs = defaultdict(set)
	bgc_to_sample = {}
	sample_protein_iter = defaultdict(lambda: 1000000)
	try:
		with open(bgc_listing_file) as oalf:
			for line in oalf:
				line = line.strip()
				sample, og_bgc_genbank = line.split('\t')
				prot_fasta = proteomes_directory + sample + '.faa'

				if not os.path.isfile(prot_fasta):
					logObject.warning(
						"Issues with finding proteome for sample %s, thus skipping inclusion of BGCs belonging to this sample." % sample)
					continue

				loc_to_tag = {}
				seq_to_tag = {}
				sample_lt = None
				with open(prot_fasta) as opf:
					for rec in SeqIO.parse(opf, 'fasta'):
						location_tuple = tuple(list(rec.description.split()[1:-1]))
						loc_to_tag[location_tuple] = rec.id
						seq_to_tag[str(rec.seq).replace('*', '')] = rec.id
						sample_lt = rec.id.split('_')[0]

				cp_bgc_sample_dir = antismash_bgcs_directory + sample + '/'
				if not os.path.isdir(cp_bgc_sample_dir):
					os.system('mkdir %s' % cp_bgc_sample_dir)

				cp_bgc_genbank = cp_bgc_sample_dir + og_bgc_genbank.split('/')[-1]
				if not og_bgc_genbank in bgc_mappings: continue
				sample_bgcs[sample].add(cp_bgc_genbank)
				bgc_to_sample[cp_bgc_genbank] = sample
				cp_bgc_genbank_handle = open(cp_bgc_genbank, 'w')

				scaff_id, scaff_start = bgc_mappings[og_bgc_genbank]
				with open(og_bgc_genbank) as obgf:
					for rec in SeqIO.parse(obgf, 'genbank'):
						rec.name = scaff_id
						rec.description = 'BGC prediction on scaffold %s by %s, starts at: %d' % (
						scaff_id, bgc_prediction_software, scaff_start)
						rec.id = scaff_id
						updated_rec = copy.deepcopy(rec)
						updated_features = []
						for feature in rec.features:
							if feature.type == 'CDS':
								try:
									start = None
									end = None
									all_starts = []
									all_ends = []
									all_directions = []
									all_coords = []
									if not 'join' in str(feature.location):
										start = min(
											[int(x.strip('>').strip('<')) for x in
											 str(feature.location)[1:].split(']')[0].split(':')]) + 1
										end = max([int(x.strip('>').strip('<')) for x in
												   str(feature.location)[1:].split(']')[0].split(':')])
										direction = str(feature.location).split('(')[1].split(')')[0]
										all_starts.append(start)
										all_ends.append(end)
									else:
										all_starts = []
										all_ends = []
										all_directions = []
										for exon_coord in str(feature.location)[5:-1].split(', '):
											start = min([int(x.strip('>').strip('<')) for x in
														 exon_coord[1:].split(']')[0].split(':')]) + 1
											end = max([int(x.strip('>').strip('<')) for x in
													   exon_coord[1:].split(']')[0].split(':')])
											direction = exon_coord.split('(')[1].split(')')[0]
											all_starts.append(start)
											all_ends.append(end)
									start = min(all_starts)
									end = max(all_ends)
									start = scaff_start + start
									end = scaff_start + end
									loc_tuple = tuple([scaff_id, str(start), str(end)])
									prot_seq = str(feature.qualifiers.get('translation')[0]).replace('*', '')
									prot_id = loc_to_tag[loc_tuple]
									prot_id_by_seq = seq_to_tag[prot_seq]
									assert (prot_id_by_seq == prot_id)
									feature.qualifiers.get('locus_tag')[0] = prot_id
								except:
									feature.qualifiers.get('locus_tag')[0] = sample_lt + '_' + str(
										sample_protein_iter[sample_lt])
									sample_protein_iter[sample_lt] += 1
								updated_features.append(feature)
							else:
								updated_features.append(feature)
						updated_rec.features = updated_features
						SeqIO.write(updated_rec, cp_bgc_genbank_handle, 'genbank')
				cp_bgc_genbank_handle.close()
	except Exception as e:
		logObject.error("Problem with parsing BGC Genbank listings.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return ([sample_bgcs, bgc_to_sample])


def extractProteinsFromBGCs(sample_bgcs, bgc_prot_directory, logObject):
	sample_bgc_prots = defaultdict(lambda: defaultdict(set))
	try:
		for sample in sample_bgcs:
			samp_bgc_prot_file = bgc_prot_directory + sample + '.faa'
			samp_bgc_prot_handle = open(samp_bgc_prot_file, 'w')
			for bgc in sample_bgcs[sample]:
				with open(bgc) as obf:
					for rec in SeqIO.parse(obf, 'genbank'):
						scaff_id = rec.id
						scaff_start = int(rec.description.split()[-1].strip())
						for feature in rec.features:
							if feature.type == 'CDS':
								prot_lt = feature.qualifiers.get('locus_tag')[0]
								start = None
								end = None
								direction = None
								all_starts = []
								all_ends = []
								all_directions = []
								if not 'join' in str(feature.location):
									start = min(
										[int(x.strip('>').strip('<')) for x in
										 str(feature.location)[1:].split(']')[0].split(':')]) + 1
									end = max([int(x.strip('>').strip('<')) for x in
											   str(feature.location)[1:].split(']')[0].split(':')])
									direction = str(feature.location).split('(')[1].split(')')[0]
									all_starts.append(start);
									all_ends.append(end);
									all_directions.append(direction)
								else:
									all_starts = []
									all_ends = []
									all_directions = []
									for exon_coord in str(feature.location)[5:-1].split(', '):
										start = min([int(x.strip('>').strip('<')) for x in
													 exon_coord[1:].split(']')[0].split(':')]) + 1
										end = max([int(x.strip('>').strip('<')) for x in
												   exon_coord[1:].split(']')[0].split(':')])
										direction = exon_coord.split('(')[1].split(')')[0]
										all_starts.append(start);
										all_ends.append(end);
										all_directions.append(direction)
								assert (len(set(all_directions)) == 1)
								start = min(all_starts)
								end = max(all_ends)
								direction = all_directions[0]
								prot_seq = str(feature.qualifiers.get('translation')[0]).replace('*', '')
								sample_bgc_prots[sample][bgc].add(prot_lt)
								samp_bgc_prot_handle.write('>' + ' '.join([str(x) for x in
																		   [prot_lt, scaff_id, start, end,
																			direction]]) + '\n' + prot_seq + '\n')
			samp_bgc_prot_handle.close()
	except Exception as e:
		logObject.error("Issues with parsing out protein sequences from BGC Genbanks.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return sample_bgc_prots


def performKOFamAndPGAPAnnotation(sample_bgc_proteins, bgc_prot_directory, annot_directory, kofam_hmm_file,
								  pgap_hmm_file,
								  kofam_pro_list, pgap_inf_list, logObject, cpus=1):
	sample_protein_annotations = defaultdict(lambda: defaultdict(lambda: "hypothetical protein"))
	try:
		ko_score_cutoffs = {}
		ko_score_types = {}
		ko_descriptions = defaultdict(lambda: "NA")
		with open(kofam_pro_list) as okpl:
			for i, line in enumerate(okpl):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				ko = ls[0]
				if ls[1] == '-': continue
				score_cutoff = float(ls[1])
				profile_type = ls[2]
				description = ls[-1]
				ko_score_cutoffs[ko] = score_cutoff
				ko_score_types[ko] = profile_type
				ko_descriptions[ko] = description

		pgap_descriptions = defaultdict(lambda: "NA")
		with open(pgap_inf_list) as opil:
			for i, line in enumerate(opil):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				label = ls[2]
				description = ls[10]
				pgap_descriptions[label] = description

		hmmsearch_cmds = []
		for f in os.listdir(bgc_prot_directory):
			sample = f.split('.faa')[0]
			prot_file = bgc_prot_directory + f
			ko_annot_result = annot_directory + sample + '.ko.txt'
			pgap_annot_result = annot_directory + sample + '.pgap.txt'
			hmmsearch_ko_cmd = ['hmmsearch', '--cpu', '2', '--tblout', ko_annot_result, kofam_hmm_file,
								prot_file, logObject]
			hmmsearch_pgap_cmd = ['hmmsearch', '--cpu', '2', '--tblout', pgap_annot_result, pgap_hmm_file,
								  prot_file, logObject]
			hmmsearch_cmds.append(hmmsearch_ko_cmd)
			hmmsearch_cmds.append(hmmsearch_pgap_cmd)

		p = multiprocessing.Pool(math.floor(cpus / 2))
		p.map(multiProcess, hmmsearch_cmds)
		p.close()

		for sample in sample_bgc_proteins:
			ko_annot_result = annot_directory + sample + '.ko.txt'
			pgap_annot_result = annot_directory + sample + '.pgap.txt'
			hits_per_prot = defaultdict(list)

			with open(ko_annot_result) as orf:
				for line in orf:
					if line.startswith("#"): continue
					line = line.strip()
					ls = line.split()
					ko = ls[2]
					prot_id = ls[0]
					full_eval = float(ls[4])
					full_score = float(ls[5])
					dom_score = float(ls[8])
					if ko in ko_score_types:
						if ko_score_types[ko] == 'full' and full_score >= ko_score_cutoffs[ko]:
							hits_per_prot[prot_id].append(['ko', ko, -full_score, full_eval])
						elif ko_score_types[ko] == 'domain' and dom_score >= ko_score_cutoffs[ko]:
							hits_per_prot[prot_id].append(['ko', ko, -dom_score, full_eval])
					if full_eval < 1e-10:
						hits_per_prot[prot_id].append(['ko', ko, 1E10, full_eval])

			with open(pgap_annot_result) as orf:
				for line in orf:
					if line.startswith("#"): continue
					line = line.strip()
					ls = line.split()
					pgap = ls[2]
					prot_id = ls[0]
					full_eval = float(ls[4])
					if full_eval < 1e-10:
						hits_per_prot[prot_id].append(['pgap', pgap, 1E10, full_eval])

			best_hits = defaultdict(lambda: 'hypothetical protein')
			for pi in hits_per_prot:
				for i, hit in enumerate(sorted(hits_per_prot[pi], key=lambda h: (h[2], h[3]), reverse=False)):
					if i == 0:
						conf = ' (High Confidence)'
						if hit[2] == 1E10 and hit[0] == 'ko':
							conf = ' (KO; Low Confidence : e-val < 1e-10)'
						elif hit[2] == 1E10 and hit[0] == 'pgap':
							conf = ' (PGAP; Low Confidence : e-val < 1e-10)'
						if hit[0] == 'ko':
							best_hits[pi] = hit[1] + ': ' + ko_descriptions[hit[1]] + conf
						else:
							best_hits[pi] = hit[1] + ': ' + pgap_descriptions[hit[1]] + conf

			for bgc in sample_bgc_proteins[sample]:
				for pi in sample_bgc_proteins[sample][bgc]:
					sample_protein_annotations[sample][pi] = best_hits[pi]

	except Exception as e:
		logObject.error("Issues with performing KOfam / PGAP based annotations.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return dict(sample_protein_annotations)


def runOrthoFinder2Full(bgc_prot_directory, orthofinder_outdir, logObject, cpus=1):
	result_file = orthofinder_outdir + 'Orthogroups_BGC_Comprehensive.tsv'
	try:
		orthofinder_cmd = ['orthofinder', '-f', bgc_prot_directory, '-t', str(cpus)]

		logObject.info('Running the following command: %s' % ' '.join(orthofinder_cmd))
		subprocess.call(' '.join(orthofinder_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran OrthoFinder!')
		tmp_orthofinder_dir = os.path.abspath(
			[bgc_prot_directory + 'OrthoFinder/' + f for f in os.listdir(bgc_prot_directory + 'OrthoFinder/') if
			 f.startswith('Results')][0]) + '/'

		os.system('mv %s %s' % (tmp_orthofinder_dir, orthofinder_outdir))
		main_file = orthofinder_outdir + 'Orthogroups/Orthogroups.tsv'
		singletons_file = orthofinder_outdir + 'Orthogroups/Orthogroups_UnassignedGenes.tsv'
		n0_file = orthofinder_outdir + 'Phylogenetic_Hierarchical_Orthogroups/N0.tsv'

		genomes = []
		og_genes_in_hog = defaultdict(set)
		result_handle = open(result_file, 'w')
		with open(n0_file) as n0f:
			for i, line in enumerate(n0f):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					genomes = ls[3:]
					result_handle.write('Orthogroup\t' + '\t'.join(genomes) + '\n')
				else:
					hog = ls[0].split('N0.')[1]
					og = ls[1]
					for j, gs in enumerate(ls[3:]):
						genome = genomes[j]
						for gene in gs.split(', '):
							gene = gene.strip()
							gene_genome_pair = tuple([gene, genome])
							og_genes_in_hog[og].add(gene_genome_pair)
					result_handle.write(hog + '\t' + '\t'.join(ls[3:]) + '\n')

		genomes_mf = []
		with open(main_file) as omf:
			for i, line in enumerate(omf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					genomes_mf = ls[1:]
				else:
					og = ls[0]
					printlist = [og]
					value_count = 0
					for j, gs in enumerate(ls[1:]):
						genome = genomes_mf[j]
						updated_gs = []
						for gene in gs.split(', '):
							gene_genome_pair = tuple([gene, genome])
							if not gene_genome_pair in og_genes_in_hog[og]:
								updated_gs.append(gene)
								value_count += 1
						printlist.append(', '.join(updated_gs))
					result_handle.write('\t'.join(printlist) + '\n')

		genomes_sf = []
		with open(singletons_file) as osf:
			for i, line in enumerate(osf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					genomes_sf = ls[1:]
				else:
					result_handle.write(line + '\n')
		result_handle.close()

		assert (genomes == genomes_mf and genomes == genomes_sf)
		assert (os.path.isfile(result_file))
	except Exception as e:
		logObject.error("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return result_file


def runOrthoFinder2(bgc_prot_directory, orthofinder_outdir, logObject, cpus=1):
	result_file = orthofinder_outdir + 'Orthogroups_BGC_Comprehensive.tsv'
	try:
		orthofinder_cmd = ['orthofinder', '-f', bgc_prot_directory, '-t', str(cpus), '-og']
		logObject.info('Running the following command: %s' % ' '.join(orthofinder_cmd))
		subprocess.call(' '.join(orthofinder_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran OrthoFinder!')
		tmp_orthofinder_dir = os.path.abspath(
			[bgc_prot_directory + 'OrthoFinder/' + f for f in os.listdir(bgc_prot_directory + 'OrthoFinder/') if
			 f.startswith('Results')][0]) + '/'
		os.system('mv %s %s' % (tmp_orthofinder_dir, orthofinder_outdir))
		main_file = orthofinder_outdir + 'Orthogroups/Orthogroups.tsv'
		singletons_file = orthofinder_outdir + 'Orthogroups/Orthogroups_UnassignedGenes.tsv'
		result_handle = open(result_file, 'w')
		with open(main_file) as omf:
			for line in omf:
				result_handle.write(line)
		with open(singletons_file) as osf:
			for i, line in enumerate(osf):
				if i == 0: continue
				result_handle.write(line)
		result_handle.close()
		assert (os.path.isfile(result_file))
	except Exception as e:
		logObject.error("Problem with running OrthoFinder2 cmd: %s." % ' '.join(orthofinder_cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return result_file


def determineParalogyThresholds(orthofinder_bgc_matrix_file, bgc_prot_directory, blast_directory, logObject, cpus=1):
	paralogy_thresholds = defaultdict(
		lambda: [90.0, 90.0])  # item 1 : percent identiy threshold; item 2 : query coverage threshold
	samp_hg_lts = defaultdict(lambda: defaultdict(set))
	lt_to_hg = {}
	# create homolog group fasta files
	try:
		tmp_hg_seq_dir = blast_directory + 'HG_FASTAs_tmp/'
		if not os.path.isdir(tmp_hg_seq_dir):
			os.system('mkdir %s' % tmp_hg_seq_dir)

		prot_lt_to_seq = {}
		for f in os.listdir(bgc_prot_directory):
			if not f.endswith('.faa'): continue
			with open(bgc_prot_directory + f) as obpdf:
				for rec in SeqIO.parse(obpdf, 'fasta'):
					prot_lt_to_seq[rec.id] = str(rec.seq)

		samples = []
		with open(orthofinder_bgc_matrix_file) as oobmf:
			for i, line in enumerate(oobmf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					samples = ls[1:]
				else:
					hg = ls[0]
					hg_fasta_file = tmp_hg_seq_dir + hg + '.faa'
					hg_fasta_handle = open(hg_fasta_file, 'w')
					for j, lts in enumerate(ls[1:]):
						for lt in lts.split(','):
							lt = lt.strip()
							if lt == '': continue
							samp = samples[j]
							samp_hg_lts[samp][hg].add(lt)
							lt_to_hg[lt] = hg
							hg_fasta_handle.write('>' + lt + '\n' + prot_lt_to_seq[lt] + '\n')
					hg_fasta_handle.close()
	except Exception as e:
		logObject.error("Problem with processing OrthoFinder matrix: %s." % orthofinder_bgc_matrix_file)
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	# run diamond self alignment of HG FASTAs
	try:
		tmp_hg_dia_dir = blast_directory + 'HG_Self_Alignment/'
		if not os.path.isdir(tmp_hg_dia_dir):
			os.system('mkdir %s' % tmp_hg_dia_dir)

		diamond_db_cmds = []
		diamond_bp_cmds = []
		for f in os.listdir(tmp_hg_seq_dir):
			hg = f.split('.faa')[0]
			hg_fasta_file = tmp_hg_seq_dir + f
			hg_diamond_db = tmp_hg_seq_dir + hg
			hg_diamond_tsv = tmp_hg_dia_dir + hg + '.tsv'

			diamond_makedb_cmd = ['diamond', 'makedb', '--in', hg_fasta_file, '-d', hg_diamond_db, logObject]
			diamond_blastp_cmd = ['diamond', 'blastp', '-p', '1', '-d', hg_diamond_db, '-q', hg_fasta_file, '-o',
								  hg_diamond_tsv,
								  '-f',
								  '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp',
								  logObject]

			diamond_db_cmds.append(diamond_makedb_cmd)
			diamond_bp_cmds.append(diamond_blastp_cmd)

		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, diamond_db_cmds)
		p.close()

		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, diamond_bp_cmds)
		p.close()

	except Exception as e:
		logObject.error("Problem with running Diamond self-alignment of homolog group protein instance FASTAs.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	try:
		for f in os.listdir(tmp_hg_dia_dir):
			hg = f.split('.tsv')[0]
			hg_diamond_tsv = tmp_hg_dia_dir + f
			assert (os.path.isfile(hg_diamond_tsv))
			best_alignments_between_proteins = defaultdict(lambda: [0.0, None, None])
			with open(hg_diamond_tsv) as ohdt:
				for line in ohdt:
					line = line.strip()
					ls = line.split('\t')
					query = ls[0]
					subject = ls[1]
					if query == subject: continue
					pair = '|'.join(sorted([query, subject]))
					pident = float(ls[2])
					qcovhsp = float(ls[-1])
					bitscore = float(ls[-2])
					if bitscore > best_alignments_between_proteins[pair][0]:
						best_alignments_between_proteins[pair] = [bitscore, pident, qcovhsp]

			min_bitscore = 100000.0
			for pair in best_alignments_between_proteins:
				bitscore, pident, qcovhsp = best_alignments_between_proteins[pair]
				if bitscore < min_bitscore:
					min_bitscore = bitscore
					paralogy_thresholds[hg] = [pident, qcovhsp]
	except Exception as e:
		logObject.error("Problem with parsing self alignment results." % orthofinder_bgc_matrix_file)
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	os.system('rm -rf %s %s' % (tmp_hg_seq_dir, tmp_hg_dia_dir))
	return [samp_hg_lts, lt_to_hg, paralogy_thresholds]


def updateSampleGenomesWithGenbanks(genbanks_directory):
	sample_genomes = {}
	for f in os.listdir(genbanks_directory):
		if not f.endswith('.gbk'): continue
		gbk = genbanks_directory + f
		sample = f.split('.gbk')[0]
		sample_genomes[sample] = gbk
	return sample_genomes


def incorporateBGCProteinsIntoProteomesAndGenbanks(sample_bgc_proteins, sample_genomes, protein_annotations,
												   bgc_prot_directory,
												   proteomes_directory, final_proteomes_directory,
												   final_genbanks_directory,
												   sample_listing_file, bgc_listing_file, logObject):
	sample_listing_handle = open(sample_listing_file, 'w')
	bgc_listing_handle = open(bgc_listing_file, 'w')
	try:
		for sample in sample_bgc_proteins:

			gw_prot_file = proteomes_directory + sample + '.faa'
			gw_prot_to_location = {}
			scaff_glts = defaultdict(set)
			with open(gw_prot_file) as ogpf:
				for rec in SeqIO.parse(ogpf, 'fasta'):
					gw_prot_to_location[rec.id] = [rec.description.split()[1], int(rec.description.split()[2]),
												   int(rec.description.split()[3])]
					scaff_glts[rec.description.split()[1]].add(rec.id)

			sample_lts_to_prune = set([])
			sample_lts_to_add_protein_sequences = {}
			sample_lts_to_add_genbank_features = defaultdict(list)
			for bgc in sample_bgc_proteins[sample]:

				bgc_listing_handle.write(sample + '\t' + bgc + '\n')

				bgc_prots = sample_bgc_proteins[sample][bgc]
				bgc_prots_1x = set([x for x in bgc_prots if x.split('_')[1][0] == '1'])

				scaff_start = None
				with open(bgc) as obf:
					for rec in SeqIO.parse(obf, 'genbank'):
						scaff_id = rec.id
						scaff_start = int(rec.description.split()[-1].strip())
						for feature in rec.features:
							if feature.type == 'CDS':
								prot_lt = feature.qualifiers.get('locus_tag')[0]
								if prot_lt in bgc_prots_1x:
									all_coords = []
									if not 'join' in str(feature.location):
										start = scaff_start + min(
											[int(x.strip('>').strip('<')) for x in
											 str(feature.location)[1:].split(']')[0].split(':')]) + 1
										end = scaff_start + max([int(x.strip('>').strip('<')) for x in
																 str(feature.location)[1:].split(']')[0].split(':')])
										direction = str(feature.location).split('(')[1].split(')')[0]
										all_coords.append([start, end, direction])
									else:
										for exon_coord in str(feature.location)[5:-1].split(', '):
											start = scaff_start + min([int(x.strip('>').strip('<')) for x in
																	   exon_coord[1:].split(']')[0].split(':')]) + 1
											end = scaff_start + max([int(x.strip('>').strip('<')) for x in
																	 exon_coord[1:].split(']')[0].split(':')])
											direction = exon_coord.split('(')[1].split(')')[0]
											all_coords.append([start, end, direction])
									fls = []
									for sc, ec, dc in all_coords:
										dir = 1
										if dc == '-': dir = -1
										fls.append(FeatureLocation(sc - 1, ec, strand=dir))
									feat_loc = fls[0]
									if len(fls) > 1:
										feat_loc = sum(fls)
									feature.location = feat_loc
									sample_lts_to_add_genbank_features[scaff_id].append(feature)

				bgc_prot_file = bgc_prot_directory + sample + '.faa'
				bgc_prot_to_location = {}
				with open(bgc_prot_file) as obpf:
					for rec in SeqIO.parse(obpf, 'fasta'):
						bgc_prot_to_location[rec.id] = [rec.description.split()[1],
														int(rec.description.split()[2]),
														int(rec.description.split()[3])]
						if rec.id in bgc_prots_1x:
							sample_lts_to_add_protein_sequences[rec.id] = [
								rec.id + ' ' + rec.description.split()[1] + ' ' + str(
									int(rec.description.split()[2]) + scaff_start) + ' ' + str(
									int(rec.description.split()[3]) + scaff_start), str(rec.seq)]

				for blt in bgc_prots_1x:
					blt_scaff, blt_start, blt_end = bgc_prot_to_location[blt]
					blt_range = set(range(blt_start + scaff_start, blt_end + scaff_start + 1))
					for glt in scaff_glts[blt_scaff]:
						glt_scaff, glt_start, glt_end = gw_prot_to_location[glt]
						glt_range = set(range(glt_start, glt_end + 1))
						if glt_scaff == blt_scaff and float(len(blt_range.intersection(glt_range))) / float(
								len(glt_range)) >= 0.25:
							if not glt in bgc_prots:
								sample_lts_to_prune.add(glt)

			final_gw_sample_faa = final_proteomes_directory + sample + '.faa'
			final_gw_sample_gbk = final_genbanks_directory + sample + '.gbk'

			sample_listing_handle.write(sample + '\t' + final_gw_sample_gbk + '\t' + final_gw_sample_faa + '\n')

			faa_handle = open(final_gw_sample_faa, 'w')
			with open(proteomes_directory + sample + '.faa') as of:
				for rec in SeqIO.parse(of, 'fasta'):
					if not rec.id in sample_lts_to_prune:
						faa_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
			for rec in sample_lts_to_add_protein_sequences.items():
				faa_handle.write('>' + rec[1][0] + '\n' + str(rec[1][1]) + '\n')
			faa_handle.close()

			gbk_handle = open(final_gw_sample_gbk, 'w')
			og = None
			if sample_genomes[sample].endswith('.gz'):
				og = gzip.open(sample_genomes[sample], 'rt')
			else:
				og = open(sample_genomes[sample])
			for rec in SeqIO.parse(og, 'genbank'):
				updated_features = []
				cds_iter = 0
				starts = []
				for feature in rec.features:
					if feature.type == 'CDS':
						prot_lt = feature.qualifiers.get('locus_tag')[0]
						if protein_annotations == None:
							feature.qualifiers['product'] = 'hypothetical protein'
						else:
							feature.qualifiers['product'] = [protein_annotations[sample][prot_lt]]
						feature.qualifiers.move_to_end('translation')
						if prot_lt in sample_lts_to_prune: continue
						updated_features.append(feature)
						all_starts = []
						if not 'join' in str(feature.location):
							start = min([int(x.strip('>').strip('<')) for x in
										 str(feature.location)[1:].split(']')[0].split(':')]) + 1
							all_starts.append(start)
						else:
							for exon_coord in str(feature.location)[5:-1].split(', '):
								start = min(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
								all_starts.append(start)
						start = min(all_starts)
						starts.append([cds_iter, start])
						cds_iter += 1
				for feature in sample_lts_to_add_genbank_features[rec.id]:
					if feature.type == 'CDS':
						prot_lt = feature.qualifiers.get('locus_tag')[0]
						if protein_annotations == None:
							feature.qualifiers['product'] = 'hypothetical protein'
						else:
							feature.qualifiers['product'] = [protein_annotations[sample][prot_lt]]
						feature.qualifiers.move_to_end('translation')
						updated_features.append(feature)
						all_starts = []
						if not 'join' in str(feature.location):
							start = min([int(x.strip('>').strip('<')) for x in
										 str(feature.location)[1:].split(']')[0].split(':')]) + 1
							all_starts.append(start)
						else:
							for exon_coord in str(feature.location)[5:-1].split(', '):
								start = min(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
								all_starts.append(start)
						start = min(all_starts)
						starts.append([cds_iter, start])
						cds_iter += 1
				sorted_updated_features = []
				for sort_i in sorted(starts, key=itemgetter(1)):
					sorted_updated_features.append(updated_features[sort_i[0]])
				rec.features = sorted_updated_features
				SeqIO.write(rec, gbk_handle, 'genbank')
		if og != None:
			og.close()
		gbk_handle.close()

	except Exception as e:
		logObject.error("Problem with creating updated.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	sample_listing_handle.close()
	bgc_listing_handle.close()


def identifyParalogsAndCreateResultFiles(samp_hg_lts, lt_to_hg, sample_bgc_proteins, paralogy_thresholds,
										 protein_annotations,
										 bgc_prot_directory, blast_directory, proteomes_directory, sample_genomes,
										 final_proteomes_directory,
										 final_genbanks_directory, sample_listing_file, bgc_listing_file,
										 orthofinder_matrix_file, logObject, cpus=1):
	sample_listing_handle = open(sample_listing_file, 'w')
	bgc_listing_handle = open(bgc_listing_file, 'w')
	try:
		tmp_diamond_dir = blast_directory + 'Diamond_BGC_to_GenomeWide_Proteomes_tmp/'
		if not os.path.isdir(tmp_diamond_dir): os.system('mkdir %s' % tmp_diamond_dir)
		diamond_db_cmds = []
		diamond_bp_cmds = []
		for f in os.listdir(proteomes_directory):
			if not f.endswith('.faa'): continue
			sample = f.split('.faa')[0]
			gw_prot_file = proteomes_directory + f
			bgc_prot_file = bgc_prot_directory + sample + '.faa'
			diamond_db = tmp_diamond_dir + sample
			diamond_tsv = tmp_diamond_dir + sample + '.tsv'

			diamond_makedb_cmd = ['diamond', 'makedb', '--in', gw_prot_file, '-d', diamond_db, logObject]
			diamond_blastp_cmd = ['diamond', 'blastp', '-p', '1', '-d', diamond_db, '-q', bgc_prot_file, '-o',
								  diamond_tsv, '-f',
								  '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp',
								  logObject]

			diamond_db_cmds.append(diamond_makedb_cmd)
			diamond_bp_cmds.append(diamond_blastp_cmd)

		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, diamond_db_cmds)
		p.close()

		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, diamond_bp_cmds)
		p.close()

	except Exception as e:
		logObject.error(
			"Problem with running diamond alignments between BGC proteins and genome-wide proteomes for samples.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	try:
		for f in os.listdir(tmp_diamond_dir):
			if not f.endswith('.tsv'): continue
			sample = f.split('.tsv')[0]
			diamond_tsv_file = tmp_diamond_dir + f
			hg_queries = defaultdict(set)
			hg_subjects = defaultdict(set)
			with open(diamond_tsv_file) as odtf:
				for line in odtf:
					line = line.strip()
					ls = line.split('\t')
					query = ls[0]
					query_hg = lt_to_hg[query]
					subject = ls[1]
					pident = float(ls[2])
					qcovhsp = float(ls[-1])
					hg_queries[query_hg].add(query)
					if pident >= paralogy_thresholds[query_hg][0] and qcovhsp >= paralogy_thresholds[query_hg][1]:
						hg_subjects[query_hg].add(subject)

			gw_prot_file = proteomes_directory + sample + '.faa'
			gw_prot_to_location = {}
			with open(gw_prot_file) as ogpf:
				for rec in SeqIO.parse(ogpf, 'fasta'):
					gw_prot_to_location[rec.id] = [rec.description.split()[1], int(rec.description.split()[2]),
												   int(rec.description.split()[3])]

			sample_lts_to_prune = set([])
			sample_lts_to_add_protein_sequences = {}
			sample_lts_to_add_genbank_features = defaultdict(list)
			for bgc in sample_bgc_proteins[sample]:

				bgc_listing_handle.write(sample + '\t' + bgc + '\n')

				bgc_prots = sample_bgc_proteins[sample][bgc]
				bgc_prots_1x = set([x for x in bgc_prots if x.split('_')[1][0] == '1'])

				bgc_prot_file = bgc_prot_directory + sample + '.faa'
				bgc_prot_to_location = {}
				with open(bgc_prot_file) as obpf:
					for rec in SeqIO.parse(obpf, 'fasta'):
						bgc_prot_to_location[rec.id] = [rec.description.split()[1],
														int(rec.description.split()[2]),
														int(rec.description.split()[3])]
						if rec.id in bgc_prots_1x:
							sample_lts_to_add_protein_sequences[rec.id] = [rec.description, str(rec.seq)]

				with open(bgc) as obf:
					for rec in SeqIO.parse(obf, 'genbank'):
						scaff_id = rec.id
						scaff_start = int(rec.description.split()[-1].strip())
						for feature in rec.features:
							if feature.type == 'CDS':
								prot_lt = feature.qualifiers.get('locus_tag')[0]
								if not prot_lt in bgc_prots_1x: continue
								all_coords = []
								if not 'join' in str(feature.location):
									start = scaff_start + min(
										[int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
									end = scaff_start + max(
										[int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
									direction = str(feature.location).split('(')[1].split(')')[0]
									all_coords.append([start, end, direction])
								else:
									for exon_coord in str(feature.location)[5:-1].split(', '):
										start = scaff_start + min([int(x.strip('>').strip('<')) for x in
																   exon_coord[1:].split(']')[0].split(':')]) + 1
										end = scaff_start + max([int(x.strip('>').strip('<')) for x in
																 exon_coord[1:].split(']')[0].split(':')])
										direction = exon_coord.split('(')[1].split(')')[0]
										all_coords.append([start, end, direction])
								fls = []
								for sc, ec, dc in all_coords:
									dir = 1
									if dc == '-': dir = -1
									fls.append(FeatureLocation(sc - 1, ec, strand=dir))
								feat_loc = fls[0]
								if len(fls) > 1:
									feat_loc = sum(fls)
								feature.location = feat_loc
								sample_lts_to_add_genbank_features[scaff_id].append(feature)

				for hg in hg_queries:
					bgc_hg_lts = hg_queries[hg].intersection(bgc_prots)
					gw_hg_lts = hg_subjects[hg].difference(hg_queries[hg])
					for glt in gw_hg_lts:
						glt_scaff, glt_start, glt_end = gw_prot_to_location[glt]
						glt_range = set(range(glt_start, glt_end + 1))
						for blt in bgc_hg_lts:
							if blt.split('_')[1][0] == '0': continue
							blt_scaff, blt_start, blt_end = bgc_prot_to_location[blt]
							blt_range = set(range(blt_start, blt_end + 1))
							if glt_scaff == blt_scaff and float(len(blt_range.intersection(glt_range))) / float(
									len(glt_range)) >= 0.25:
								sample_lts_to_prune.add(glt)

			final_gw_sample_faa = final_proteomes_directory + sample + '.faa'
			final_gw_sample_gbk = final_genbanks_directory + sample + '.gbk'

			sample_listing_handle.write(sample + '\t' + final_gw_sample_gbk + '\t' + final_gw_sample_faa + '\n')

			faa_handle = open(final_gw_sample_faa, 'w')
			with open(proteomes_directory + sample + '.faa') as of:
				for rec in SeqIO.parse(of, 'fasta'):
					if not rec.id in sample_lts_to_prune:
						faa_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
			for rec in sample_lts_to_add_protein_sequences.items():
				faa_handle.write('>' + rec[1][0] + '\n' + str(rec[1][1]) + '\n')
			faa_handle.close()

			gbk_handle = open(final_gw_sample_gbk, 'w')
			og = None
			if sample_genomes[sample].endswith('.gz'):
				og = gzip.open(sample_genomes[sample], 'rt')
			else:
				og = open(sample_genomes[sample])
			for rec in SeqIO.parse(og, 'genbank'):
				updated_features = []
				cds_iter = 0
				starts = []
				for feature in rec.features:
					if feature.type == 'CDS':
						prot_lt = feature.qualifiers.get('locus_tag')[0]
						if protein_annotations == None:
							feature.qualifiers['product'] = 'hypothetical protein'
						else:
							feature.qualifiers['product'] = [protein_annotations[sample][prot_lt]]
						feature.qualifiers.move_to_end('translation')
						if prot_lt in sample_lts_to_prune: continue
						updated_features.append(feature)
						all_starts = []
						if not 'join' in str(feature.location):
							start = min([int(x.strip('>').strip('<')) for x in
										 str(feature.location)[1:].split(']')[0].split(':')]) + 1
							all_starts.append(start)
						else:
							for exon_coord in str(feature.location)[5:-1].split(', '):
								start = min(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
								all_starts.append(start)
						start = min(all_starts)
						starts.append([cds_iter, start])
						cds_iter += 1
				for feature in sample_lts_to_add_genbank_features[rec.id]:
					if feature.type == 'CDS':
						prot_lt = feature.qualifiers.get('locus_tag')[0]
						if protein_annotations == None:
							feature.qualifiers['product'] = 'hypothetical protein'
						else:
							feature.qualifiers['product'] = [protein_annotations[sample][prot_lt]]
						feature.qualifiers.move_to_end('translation')
						updated_features.append(feature)
						all_starts = []
						if not 'join' in str(feature.location):
							start = min([int(x.strip('>').strip('<')) for x in
										 str(feature.location)[1:].split(']')[0].split(':')]) + 1
							all_starts.append(start)
						else:
							for exon_coord in str(feature.location)[5:-1].split(', '):
								start = min(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
								all_starts.append(start)
						start = min(all_starts)
						starts.append([cds_iter, start])
						cds_iter += 1
				sorted_updated_features = []
				for sort_i in sorted(starts, key=itemgetter(1)):
					sorted_updated_features.append(updated_features[sort_i[0]])
				rec.features = sorted_updated_features
				SeqIO.write(rec, gbk_handle, 'genbank')
			og.close()
			gbk_handle.close()

			for hg in samp_hg_lts[sample]:
				updated_lts = set(samp_hg_lts[sample][hg])
				for lt in hg_subjects[hg].difference(hg_queries[hg]):
					if not lt in sample_lts_to_prune:
						updated_lts.add(lt)
				samp_hg_lts[sample][hg] = updated_lts

	except Exception as e:
		logObject.error("Problem with creating updated.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	try:
		all_hgs = set([])
		for s in samp_hg_lts:
			for h in samp_hg_lts[s]:
				all_hgs.add(h)

		orthofinder_matrix_handle = open(orthofinder_matrix_file, 'w')
		orthofinder_matrix_handle.write('Orthogroup\t' + '\t'.join([sample for sample in sorted(samp_hg_lts)]) + '\n')
		for hg in sorted(all_hgs):
			hg_row = [hg]
			for sample in sorted(samp_hg_lts):
				hg_row.append(', '.join(samp_hg_lts[sample][hg]))
			orthofinder_matrix_handle.write('\t'.join(hg_row) + '\n')
		orthofinder_matrix_handle.close()

	except Exception as e:
		logObject.error("Problem with updating OrthoFinder2 sample vs. homolog group matrix.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	sample_listing_handle.close()
	bgc_listing_handle.close()


def createGCFListingsDirectory(sample_bgcs, bgc_to_sample, bigscape_results_dir, gcf_listings_directory, logObject):
	try:
		final_stats_file = '/'.join(gcf_listings_directory.split('/')[:-2]) + '/GCF_Details.txt'
		sf_handle = open(final_stats_file, 'w')
		sf_handle.write('\t'.join(['GCF id', 'number of BGCs', 'number of samples', 'samples with multiple BGCs in GCF',
								   'size of the SCC', 'mean number of OGs', 'stdev for number of OGs',
								   'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity',
								   'number of core gene aggregates', 'annotations']) + '\n')
		bgc_paths = {}
		for sample in sample_bgcs:
			for bgc in sample_bgcs[sample]:
				bgc_paths[bgc.split('/')[-1].split('.gbk')[0]] = [sample, bgc]

		nf_bigscape_results_dir = bigscape_results_dir + 'network_files/'
		assert (os.path.isdir(nf_bigscape_results_dir))
		newest_date = [0, 0, 0]
		newest_time = [0, 0, 0]
		most_recent_results_dir = None
		for sd in os.listdir(nf_bigscape_results_dir):
			full_sd_dir = nf_bigscape_results_dir + sd + '/'
			if not os.path.isdir(full_sd_dir): continue
			date = [int(x) for x in sd.split('_')[0].split('-')]
			time = [int(x) for x in sd.split('_')[1].split('-')]
			if date[0] > newest_date[0] or (date[0] == newest_date[0] and date[1] > newest_date[1]) or (
					date[0] == newest_date[0] and date[1] == newest_date[1] and date[2] > newest_date[2]):
				newest_date = date
				newest_time = time
				most_recent_results_dir = full_sd_dir
			elif date[0] == newest_date[0] and date[1] == newest_date[1] and date[2] == newest_date[2]:
				if time[0] > newest_time[0] or (time[0] == newest_time[0] and time[1] > newest_time[1]) or (
						time[0] == newest_time[0] and time[1] == newest_time[1] and time[2] > newest_time[2]):
					newest_date = date
					newest_time = time
					most_recent_results_dir = full_sd_dir
		assert (os.path.isdir(most_recent_results_dir))

		gcfs = defaultdict(set)
		gcf_classes = defaultdict(set)
		for sd in os.listdir(most_recent_results_dir):
			class_sd = most_recent_results_dir + sd + '/'
			if not os.path.isdir(class_sd): continue
			for f in os.listdir(class_sd):
				if not (f.endswith('.tsv') and '_clustering_' in f): continue
				bgc_class = f.split('_clustering_')[0]
				class_gcf_file = class_sd + f
				with open(class_gcf_file) as ocgf:
					for i, line in enumerate(ocgf):
						if i == 0: continue
						line = line.strip()
						ls = line.split('\t')
						if ls[0] in bgc_paths and os.path.isfile(bgc_paths[ls[0]][1]):
							gcfs[ls[1]].add(ls[0])
							gcf_classes[ls[1]].add(bgc_class)

		for gcf in gcfs:
			gcf_file = gcf_listings_directory + 'GCF_' + gcf + '.txt'
			gcf_handle = open(gcf_file, 'w')
			samples_with_gcf = defaultdict(int)
			for b in gcfs[gcf]:
				bs = bgc_to_sample[bgc_paths[b][1]]
				samples_with_gcf[bs] += 1
				gcf_handle.write(bgc_paths[b][0] + '\t' + bgc_paths[b][1] + '\n')
			gcf_handle.close()
			samples_with_multiple_bgcs = 0
			for s in samples_with_gcf:
				if samples_with_gcf[s] > 1:
					samples_with_multiple_bgcs += 1
			"""
			'GCF id', 'number of BGCs', 'number of samples', 'samples with multiple BGCs in GCF',
					   'size of the SCC', 'mean number of OGs', 'stdev for number of OGs',
					   'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity',
					   'number of core gene aggregates', 'annotations']) + '\n')
			"""
			sf_handle.write('\t'.join([str(x) for x in ['GCF_' + gcf, len(gcfs[gcf]), len(samples_with_gcf),
														samples_with_multiple_bgcs, 'NA', 'NA', 'NA', 'NA', 'NA',
														'NA', '; '.join(sorted(gcf_classes[gcf]))]]) + '\n')
		sf_handle.close()
	except Exception as e:
		logObject.error("Problem with parsing BiG-SCAPE results directory provided.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())


def updateBGCGenbanksToIncludeAnnotations(protein_annotations, bgc_to_sample, sample_bgc_proteins,
										  antismash_bgcs_directory, antismash_bgcs_directory_updated, logObject):
	sample_bgcs_updated = defaultdict(set)
	bgc_to_sample_updated = {}
	sample_bgc_proteins_update = defaultdict(lambda: defaultdict(set))
	try:
		for s in os.listdir(antismash_bgcs_directory):
			os.system('mkdir %s' % (antismash_bgcs_directory_updated + s))
			for f in os.listdir(antismash_bgcs_directory + s + '/'):
				if not f.endswith('.gbk'): continue
				sample = bgc_to_sample[antismash_bgcs_directory + s + '/' + f]
				update_gbk_file = antismash_bgcs_directory_updated + s + '/' + f
				ugf_handle = open(update_gbk_file, 'w')
				with open(antismash_bgcs_directory + s + '/' + f) as oabduf:
					for rec in SeqIO.parse(oabduf, 'genbank'):
						updated_features = []
						for feature in rec.features:
							if feature.type == "CDS":
								prot_lt = feature.qualifiers.get('locus_tag')[0]
								feature.qualifiers['product'] = [protein_annotations[sample][prot_lt]]
								feature.qualifiers.move_to_end('translation')
								updated_features.append(feature)
							else:
								updated_features.append(feature)
						rec.features = updated_features
						SeqIO.write(rec, ugf_handle, 'genbank')
				ugf_handle.close()
				sample_bgcs_updated[s].add(update_gbk_file)
				bgc_to_sample_updated[update_gbk_file] = s
		for s in sample_bgc_proteins:
			for bgc in sample_bgc_proteins[s]:
				bgc_prots = sample_bgc_proteins[s][bgc]
				new_bgc = antismash_bgcs_directory_updated + '/'.join(bgc.split('/')[-2:])
				sample_bgc_proteins_update[s][new_bgc] = bgc_prots

	except Exception as e:
		logObject.error("Problem with updating AntiSMASH BGC Genbanks to feature KOfam annotations.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return ([sample_bgcs_updated, bgc_to_sample_updated, sample_bgc_proteins_update])


def selectFinalResultsAndCleanUp(outdir, fin_outdir, logObject):
	try:
		delete_set = set(['BLASTing_of_Ortholog_Groups', 'OrthoFinder2_Results', 'KOfam_Annotations',
						  'Prodigal_Gene_Calling_Additional', 'Predicted_Proteomes_Initial',
						  'Prodigal_Gene_Calling', 'Genomic_Genbanks_Initial'])
		for fd in os.listdir(outdir):
			if os.path.isfile(fd): continue
			subdir = outdir + fd + '/'
			if fd in delete_set:
				os.system('rm -rf %s' % subdir)
		if os.path.isdir(outdir + 'BGCs_Retagged_and_Updated'):
			os.system('rm -rf %s' % outdir + 'BGCs_Retagged')
		if os.path.isdir(outdir + 'lsaBGC_AutoExpansion_Results/'):
			os.system('ln -s %s %s' % (
			outdir + 'lsaBGC_AutoExpansion_Results/Updated_GCF_Listings/', fin_outdir + 'Expanded_GCF_Listings'))
			os.system('ln -s %s %s' % (
			outdir + 'lsaBGC_AutoExpansion_Results/Orthogroups.expanded.tsv', fin_outdir + 'Expanded_Orthogroups.tsv'))
			os.system('ln -s %s %s' % (outdir + 'lsaBGC_AutoExpansion_Results/Sample_Annotation_Files.txt',
									   fin_outdir + 'Expanded_Sample_Annotation_Files.txt'))
			if os.path.isfile(outdir + 'Intermediate_Files/GToTree_output.tre'):
				os.system('ln -s %s %s' % (outdir + 'Intermediate_Files/GToTree_output.tre', fin_outdir))
			if os.path.isfile(outdir + 'Intermediate_Files/GToTree_Expected_Similarities.txt'):
				os.system('ln -s %s %s' % (outdir + 'Intermediate_Files/GToTree_Expected_Similarities.txt', fin_outdir))
			if os.path.isfile(outdir + 'Intermediate_Files/Samples_in_GToTree_Tree.txt'):
				os.system('ln -s %s %s' % (outdir + 'Intermediate_Files/Samples_in_GToTree_Tree.txt', fin_outdir))
		else:
			os.system('ln -s %s %s' % (outdir + 'Intermediate_Files/*', fin_outdir))
			if os.path.isdir(outdir + 'BiG_SCAPE_Results_Reformatted/'):
				os.system('ln -s %s %s' % (
				outdir + 'BiG_SCAPE_Results_Reformatted/GCF_Listings/', fin_outdir + 'GCF_Listings'))
				os.system('ln -s %s %s' % (
				outdir + 'BiG_SCAPE_Results_Reformatted/GCF_Details.txt', fin_outdir + 'GCF_Details.txt'))
			elif os.path.isdir(outdir + 'lsaBGC_Cluster_Results/'):
				os.system(
					'ln -s %s %s' % (outdir + 'lsaBGC_Cluster_Results/GCF_Listings/', fin_outdir + 'GCF_Listings'))
				os.system('ln -s %s %s' % (
				outdir + 'lsaBGC_Cluster_Results/GCF_Details_Expanded_Singletons.txt', fin_outdir + 'GCF_Details.txt'))
	except Exception as e:
		raise RuntimeError(traceback.format_exc())


def setupReadyDirectory(directories):
	try:
		assert (type(directories) is list)
		for d in directories:
			if os.path.isdir(d):
				os.system('rm -rf %s' % d)
			os.system('mkdir %s' % d)
	except Exception as e:
		raise RuntimeError(traceback.format_exc())


def p_adjust_bh(p):
	"""
	Benjamini-Hochberg p-value correction for multiple hypothesis testing.
	"""
	p = np.asfarray(p)
	by_descend = p.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	return q[by_orig]


def is_newick(newick):
	"""
	Function to validate if Newick phylogeny file is correctly formatted.
	"""
	try:
		t = Tree(newick)
		return True
	except:
		return False


def is_fastq(fastq):
	"""
	Function to validate if FASTA file is correctly formatted.
	"""
	try:
		with open(fastq) as of:
			SeqIO.parse(of, 'fastq')
		return True
	except:
		return False


def is_fasta(fasta):
	"""
	Function to validate if FASTA file is correctly formatted.
	"""
	try:
		if fasta.endswith('.gz'):
			with gzip.open(fasta, 'rt') as ogf:
				SeqIO.parse(ogf, 'fasta')
		else:
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
		assert (gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.gbk.gz') or gbk.endswith('.gbff.gz'))
		if gbk.endswith('.gz'):
			with gzip.open(gbk, 'rt') as ogf:
				SeqIO.parse(ogf, 'genbank')
		else:
			with open(gbk) as of:
				SeqIO.parse(of, 'genbank')
		return True
	except:
		return False

def parseVersionFromSetupPy():
	"""
	Parses version from setup.py program.
	"""

	setup_py_prog = main_dir + 'setup.py'
	version = 'NA'
	with open(setup_py_prog) as osppf:
		for line in osppf:
			line = line.strip()
			if line.startswith('version='):
				version = line.split('version=')[1]
	return version

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
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	logger.addHandler(fh)
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


def calculateTajimasD(sequences):
	"""
	The code for this functionality was largely taken from Tom Whalley's Tajima's D implementation in Python and further
	modified/corrected based on Wikipedia's page for Tajima's D (Mathematical details).
	"""

	"""Calculate pi"""
	numseqs = len(sequences)
	divisor = math.comb(numseqs, 2)
	combos = itertools.combinations(sequences, 2)
	differences = 0
	for pair in combos:
		seqA = pair[0]
		seqB = pair[1]
		for p, a in enumerate(seqA):
			b = seqB[p]
			if a != b and a != '-' and b != '-':
				differences += 1
	pi = float(differences) / divisor

	"""Calculate s, number of segregation sites)."""
	# Assume if we're in here seqs have already been checked
	combos = itertools.combinations(sequences, 2)
	indexes = set([])
	for pair in combos:
		seqA = pair[0]
		seqB = pair[1]
		for idx, (i, j) in enumerate(zip(seqA, seqB)):
			if i != j and i != '-' and j != '-':
				indexes.add(idx)

	indexes = list(indexes)
	S = len(indexes)

	"""
	Now we have pi (pairwise differences) and s (number
	of segregating sites). This gives us 'little d', so
	now we need to divide it by sqrt of variance.
	"""
	l = len(sequences)

	# calculate D
	a1 = sum([(1.0 / float(i)) for i in range(1, l)])
	a2 = sum([(1.0 / (i ** 2)) for i in range(1, l)])

	b1 = float(l + 1) / (3 * (l - 1))
	b2 = float(2 * ((l ** 2) + l + 3)) / (9 * l * (l - 1))

	c1 = b1 - (1.0 / a1)
	c2 = b2 - (float(l + 2) / (a1 * l)) + (float(a2) / (a1 ** 2.0))

	e1 = float(c1) / a1
	e2 = float(c2) / ((a1 ** 2) + a2)
	if S >= 3:
		D = (float(pi - (float(S) / a1)) / math.sqrt((e1 * S) + ((e2 * S) * (S - 1))))
		return (D)
	else:
		return ("NA")


def is_numeric(x):
	try:
		x = float(x)
		return True
	except:
		return False


def castToNumeric(x):
	try:
		x = float(x)
		return (x)
	except:
		return float('nan')


numeric_columns = set(['GCF Count', 'hg order index', 'hg consensus direction', 'median gene length',
					   'proportion of samples with hg', 'proportion of total populations with hg',
					   'hg median copy count', 'num of hg instances', 'samples with hg', 'ambiguous sites proporition',
					   'Tajimas D', 'proportion variable sites', 'proportion nondominant major allele',
					   'median beta rd',
					   'median dn ds', 'mad dn ds', 'populations with hg', 'proportion of total populations with hg',
					   'most significant Fisher exact pvalues presence absence', 'median Tajimas D per population',
					   'mad Tajimas D per population'])


def loadSampleToGCFIntoPandaDataFrame(gcf_listing_dir):
	import pandas as pd
	panda_df = None
	try:
		data = []
		data.append(['GCF', 'Sample', 'BGC Instances'])
		for f in os.listdir(gcf_listing_dir):
			gcf = f.split('.txt')[0]
			sample_counts = defaultdict(int)
			with open(gcf_listing_dir + f) as ogldf:
				for line in ogldf:
					line = line.strip()
					sample, bgc_path = line.split('\t')
					sample_counts[sample] += 1
			for s in sample_counts:
				data.append([gcf, s, sample_counts[s]])

		panda_dict = {}
		for ls in zip(*data):
			key = ' '.join(ls[0].split('_'))
			vals = ls[1:]
			panda_dict[key] = vals
		panda_df = pd.DataFrame(panda_dict)

	except Exception as e:
		raise RuntimeError(traceback.format_exc())
	return panda_df


def loadSamplesIntoPandaDataFrame(sample_annot_file, pop_spec_file=None):
	import pandas as pd
	panda_df = None
	try:
		if pop_spec_file != None:
			sample_pops = {}
			with open(pop_spec_file) as opsf:
				for line in opsf:
					line = line.strip()
					ls = line.split('\t')
					sample_pops[ls[0]] = ls[1]

		data = []
		if pop_spec_file != None:
			data.append(['Sample', 'Population/Clade'])
		else:
			data.append(['Sample'])
		with open(sample_annot_file) as saf:
			for line in saf:
				line = line.strip('\n')
				ls = line.split('\t')
				if pop_spec_file != None:
					pop = sample_pops[ls[0]]
					data.append([ls[0], pop])
				else:
					data.append([ls[0]])

		panda_dict = {}
		for ls in zip(*data):
			key = ' '.join(ls[0].split('_'))
			vals = ls[1:]
			panda_dict[key] = vals
		panda_df = pd.DataFrame(panda_dict)

	except Exception as e:
		raise RuntimeError(traceback.format_exc())
	return panda_df


def loadTableInPandaDataFrame(input_file):
	import pandas as pd
	panda_df = None
	try:
		data = []
		with open(input_file) as oif:
			for line in oif:
				line = line.strip('\n')
				ls = line.split('\t')
				data.append(ls)

		panda_dict = {}
		for ls in zip(*data):
			key = ' '.join(ls[0].split('_'))
			cast_vals = ls[1:]
			if key in numeric_columns:
				cast_vals = []
				for val in ls[1:]:
					cast_vals.append(castToNumeric(val))
			panda_dict[key] = cast_vals
		panda_df = pd.DataFrame(panda_dict)

	except Exception as e:
		raise RuntimeError(traceback.format_exc())
	return panda_df


def loadCustomPopGeneTableInPandaDataFrame(input_file):
	import pandas as pd
	panda_df = None
	try:
		ignore_data_cats = {'hg_median_copy_count', 'proportion_of_total_populations_with_hg',
							'proportion_variable_sites', 'proportion_nondominant_major_allele', 'median_dn_ds',
							'mad_dn_ds', 'all_domains', 'most_significant_Fisher_exact_pvalues_presence_absence',
							'median_Tajimas_D_per_population', 'mad_Tajimas_D_per_population',
							'most_negative_population_Tajimas_D', 'most_positive_population_Tajimas_D',
							'population_entropy', 'median_fst_like_estimate'}
		data = []
		with open(input_file) as oif:
			for line in oif:
				line = line.strip('\n')
				ls = line.split('\t')
				data.append(ls)

		panda_dict = {}
		for ls in zip(*data):
			if ls[0] in ignore_data_cats: continue
			if ls[0] == 'gcf_annotation':
				updated_ans = []
				for ans in ls[1:]:
					upans = []
					for an in ans.split('; '):
						if not an.startswith('NA:'):
							upans.append(an)
					updated_ans.append('|'.join(upans))
				panda_dict[' '.join(ls[0].split('_'))] = updated_ans
			elif ls[0] == 'annotation':
				key = ' '.join(ls[0].split('_'))
				cleaned_annots = []
				for val in ls[1:]:
					ca = []
					for a in val.split('; '):
						if a != 'hypothetical protein':
							ca.append(a)
					if len(ca) == 0:
						cleaned_annots.append('hypothetical protein')
					else:
						cleaned_annots.append('; '.join(ca))
				panda_dict[key] = cleaned_annots
			else:
				key = ' '.join(ls[0].split('_'))
				cast_vals = ls[1:]
				if key in numeric_columns:
					cast_vals = []
					for val in ls[1:]:
						cast_vals.append(castToNumeric(val))
				panda_dict[key] = cast_vals
		panda_df = pd.DataFrame(panda_dict)

	except Exception as e:
		raise RuntimeError(traceback.format_exc())
	return panda_df


