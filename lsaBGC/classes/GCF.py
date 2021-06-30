import os
import sys
import logging
import traceback
import statistics
import random
import subprocess
import pysam
import multiprocessing
from scipy.stats import f_oneway, fisher_exact, pearsonr, median_absolute_deviation
from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from operator import itemgetter
from collections import defaultdict
from lsaBGC.classes.Pan import Pan
from lsaBGC import util
from pingouin import welch_anova
from pandas import DataFrame
from pomegranate import *
import math

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-3])
RSCRIPT_FOR_BGSEE = lsaBGC_main_directory + '/lsaBGC/Rscripts/bgSee.R'
RSCRIPT_FOR_CLUSTER_ASSESSMENT_PLOTTING = lsaBGC_main_directory + '/lsaBGC/Rscripts/generatePopGenePlots.R'
RSCRIPT_FOR_TAJIMA = lsaBGC_main_directory + '/lsaBGC/Rscripts/calculateTajimasD.R'

class GCF(Pan):
	def __init__(self, bgc_genbanks_listing, gcf_id='GCF_X', logObject=None, lineage_name='Unnamed lineage'):
		super().__init__(bgc_genbanks_listing, lineage_name=lineage_name, logObject=logObject)
		self.gcf_id = gcf_id

		#######
		## Variables not set during initialization
		#######

		# General variables
		self.hg_to_color = None
		self.hg_order_scores = defaultdict(int)
		self.specific_core_homologs =set([])
		self.scc_homologs = set([])
		self.core_homologs = set([])
		self.protocluster_core_homologs = set([])

		# Sequence and alignment directories
		self.nucl_seq_dir = None
		self.prot_seq_dir = None
		self.prot_alg_dir = None
		self.codo_alg_dir = None
		self.nucl_filt_seq_dir = None
		self.prot_filt_seq_dir = None
		self.prot_filt_alg_dir = None
		self.codo_filt_alg_dir = None

		# Concatenated HMMER3 HMM profiles database of homolog groups in GCF
		self.concatenated_profile_HMM = None

		# Dictionary of individual genes to haplotype/allelic representative gene
		self.instance_to_haplotype = {}

		# Set of samples with sequencing reads to avoid for reporting alleles and novel SNVs,
		# these samples do not exhibit enough support for harboring a full BGC for the GC
		self.avoid_samples = set([])

	def identifyKeyHomologGroups(self):
		try:
			initial_samples_with_at_least_one_gcf_hg = set([])
			for hg in self.hg_genes:
				for gene in self.hg_genes[hg]:
					if len(gene.split('_')[0]) == 3:
						gene_info = self.comp_gene_info[gene]
						bgc_id = gene_info['bgc_name']
						sample_id = self.bgc_sample[bgc_id]
						initial_samples_with_at_least_one_gcf_hg.add(sample_id)

			for hg in self.hg_genes:
				sample_counts = defaultdict(int)
				sample_with_hg_as_protocluster_core = 0
				for gene in self.hg_genes[hg]:
					if len(gene.split('_')[0]) == 3:
						gene_info = self.comp_gene_info[gene]
						bgc_id = gene_info['bgc_name']
						sample_id = self.bgc_sample[bgc_id]
						sample_counts[sample_id] += 1
						if gene_info['core_overlap']:
							sample_with_hg_as_protocluster_core += 1

				samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
				samples_with_any_copy = set([s[0] for s in sample_counts.items() if s[1] > 0])

				# check that hg is single-copy-core or just core
				if len(samples_with_single_copy.symmetric_difference(initial_samples_with_at_least_one_gcf_hg)) == 0:
					self.scc_homologs.add(hg)
				if len(samples_with_any_copy.symmetric_difference(initial_samples_with_at_least_one_gcf_hg)) == 0:
					self.core_homologs.add(hg)
				if len(samples_with_any_copy) > 0 and float(sample_with_hg_as_protocluster_core)/len(samples_with_any_copy) >= 0.5:
					self.protocluster_core_homologs.add(hg)

		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with identifying key homolog groups.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def modifyPhylogenyForSamplesWithMultipleBGCs(self, input_phylogeny, result_phylogeny):
		"""
		Function which takes in an input phylogeny and produces a replicate resulting phylogeny with samples/leafs which
		have multiple BGC instances for a GCF expanded.

		:param input_phylogeny: input newick phylogeny file
		:result result_phylogeny: resulting newick phylogeny file
		"""
		try:
			number_of_added_leaves = 0
			t = Tree(input_phylogeny)
			for node in t.traverse('postorder'):
				if node.name in self.sample_bgcs and len(self.sample_bgcs[node.name]) > 1:
					og_node_name = node.name
					node.name = node.name + '_INNERNODE'
					for bgc_id in self.sample_bgcs[og_node_name]:
						# if bgc_id == node.name: continue
						node.add_child(name=bgc_id)
						child_node = t.search_nodes(name=bgc_id)[0]
						child_node.dist = 0
						if bgc_id != og_node_name: number_of_added_leaves += 1

			t.write(format=0, outfile=result_phylogeny)
			if self.logObject:
				self.logObject.info(
					"New phylogeny with an additional %d leafs to reflect samples with multiple BGCs can be found at: %s." % (
						number_of_added_leaves, result_phylogeny))
		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Had difficulties properly editing phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def assignColorsToHGs(self, gene_to_hg, bgc_genes):
		"""
		Simple function to associate each homolog group with a color for consistent coloring.

		:param gene_to_hg: gene to HG relationship.
		:param bgc_genes:  set of genes per HG.
		:return: dictionary mapping each HG to a hex color value.
		"""

		hg_bgc_counts = defaultdict(int)
		for b in bgc_genes:
			for g in bgc_genes[b]:
				if g in gene_to_hg:
					hg_bgc_counts[gene_to_hg[g]] += 1

		hgs = set([])
		for c in hg_bgc_counts:
			if hg_bgc_counts[c] > 1:
				hgs.add(c)

		# read in list of colors
		dir_path = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/'
		colors_file = dir_path + 'other/colors_200.txt'
		colors = []
		with open(colors_file) as ocf:
			colors = [x.strip() for x in ocf.readlines()]
		random.shuffle(colors)

		hg_to_color = {}
		for i, c in enumerate(set(hgs)):
			hg_to_color[c] = colors[i]
		print(hg_to_color); self.hg_to_color = hg_to_color

	def createItolBGCSeeTrack(self, result_track_file):
		"""
		Function to create a track file for visualizing BGC gene architecture across a phylogeny in the interactive tree
		of life (iTol)

		:param result_track_file: The path to the resulting iTol track file for BGC gene visualization.
		"""
		try:
			track_handle = open(result_track_file, 'w')

			if self.logObject:
				self.logObject.info("Writing iTol track file to: %s" % result_track_file)
				self.logObject.info("Track will have label: %s" % self.gcf_id)

			# write header for iTol track file
			track_handle.write('DATASET_DOMAINS\n')
			track_handle.write('SEPARATOR TAB\n')
			track_handle.write('DATASET_LABEL\t%s\n' % self.gcf_id)
			track_handle.write('COLOR\t#000000\n')
			track_handle.write('BORDER_WIDTH\t1\n')
			track_handle.write('BORDER_COLOR\t#000000\n')
			track_handle.write('SHOW_DOMAIN_LABELS\t0\n')
			track_handle.write('DATA\n')

			# write the rest of the iTol track file for illustrating genes across BGC instances
			ref_hg_directions = {}
			bgc_gene_counts = defaultdict(int)
			for bgc in self.bgc_genes:
				bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])

			for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
				bgc = item[0]
				curr_bgc_genes = self.bgc_genes[bgc]
				last_gene_end = max([self.comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
				printlist = [bgc, str(last_gene_end)]
				hg_directions = {}
				hg_lengths = defaultdict(list)
				for lt in curr_bgc_genes:
					ginfo = self.comp_gene_info[lt]
					hg = 'singleton'
					if lt in self.gene_to_hg:
						hg = self.gene_to_hg[lt]
					shape = 'None'
					if ginfo['direction'] == '+':
						shape = 'TR'
					elif ginfo['direction'] == '-':
						shape = 'TL'
					gstart = ginfo['start']
					gend = ginfo['end']
					hg_color = "#dbdbdb"
					if hg in self.hg_to_color:
						hg_color = self.hg_to_color[hg]
					gene_string = '|'.join([str(x) for x in [shape, gstart, gend, hg_color, hg]])
					printlist.append(gene_string)
					if hg != 'singleton':
						hg_directions[hg] = ginfo['direction']
						hg_lengths[hg].append(gend - gstart)
				if i == 0:
					ref_hg_directions = hg_directions
					track_handle.write('\t'.join(printlist) + '\n')
				else:
					flip_support = 0
					keep_support = 0
					for c in ref_hg_directions:
						if not c in hg_directions: continue
						hg_weight = statistics.mean(hg_lengths[c])
						if hg_directions[c] == ref_hg_directions[c]:
							keep_support += hg_weight
						else:
							flip_support += hg_weight

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
		except Exception as e:
			if self.logObject:
				self.logObject.error("Had difficulties creating iTol track for visualization of BGC gene architecture.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def visualizeGCFViaR(self, gggenes_track_file, heatmap_track_file, phylogeny_file, result_pdf_file):
		"""
		Function to create tracks for visualization of gene architecture of BGCs belonging to GCF and run Rscript bgSee.R
		to produce automatic PDFs of plots. In addition, bgSee.R also produces a heatmap to more easily identify homolog
		groups which are conserved across isolates found to feature GCF.

		:param gggenes_track_file: Path to file with gggenes track information (will be created/written to by function, if it doesn't exist!)
		:param heatmap_track_file: Path to file for heatmap visual component (will be created/written to by function, if it doesn't exist!)
		:param phylogeny_file: Phylogeny to use for visualization.
		:param result_pdf_file: Path to PDF file where plots from bgSee.R will be written to.
		"""
		try:
			if os.path.isfile(gggenes_track_file) or os.path.isfile(heatmap_track_file):
				os.system('rm -f %s %s' % (gggenes_track_file, heatmap_track_file))
			gggenes_track_handle = open(gggenes_track_file, 'w')
			heatmap_track_handle = open(heatmap_track_file, 'w')
			if self.logObject:
				self.logObject.info("Writing gggenes input file to: %s" % gggenes_track_file)
				self.logObject.info("Writing heatmap input file to: %s" % heatmap_track_file)
			# write header for track files
			gggenes_track_handle.write('label\tgene\tstart\tend\tforward\tog\tog_color\n')
			heatmap_track_handle.write('label\tog\tog_presence\tog_count\n')

			ref_hg_directions = {}

			bgc_gene_counts = defaultdict(int)
			for bgc in self.bgc_genes:
				bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])

			tree_obj = Tree(phylogeny_file)
			bgc_weights = defaultdict(int)
			all_bgcs_in_tree = set([])
			for leaf in tree_obj:
				all_bgcs_in_tree.add(str(leaf).strip('\n').lstrip('-'))
				bgc_weights[str(leaf).strip('\n').lstrip('-')] += 1

			bgc_hg_presence = defaultdict(lambda: defaultdict(lambda: 'Absent'))
			hg_counts = defaultdict(int)
			for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
				bgc = item[0]
				if not bgc in all_bgcs_in_tree: continue
				curr_bgc_genes = self.bgc_genes[bgc]
				last_gene_end = max([self.comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
				printlist = []
				hg_directions = {}
				hg_lengths = defaultdict(list)
				for lt in curr_bgc_genes:
					ginfo = self.comp_gene_info[lt]
					hg = 'singleton'
					if lt in self.gene_to_hg:
						hg = self.gene_to_hg[lt]

					gstart = ginfo['start']
					gend = ginfo['end']
					forward = "FALSE"
					if ginfo['direction'] == '+': forward = "TRUE"

					hg_color = '"#dbdbdb"'
					if hg in self.hg_to_color:
						hg_color = '"' + self.hg_to_color[hg] + '"'

					gene_string = '\t'.join([str(x) for x in [bgc, lt, gstart, gend, forward, hg, hg_color]])
					printlist.append(gene_string)
					if hg != 'singleton':
						bgc_hg_presence[bgc][hg] = hg
						hg_counts[hg] += bgc_weights[bgc]
						hg_directions[hg] = ginfo['direction']
						hg_lengths[hg].append(gend - gstart)
				if i == 0:
					ref_hg_directions = hg_directions
					gggenes_track_handle.write('\n'.join(printlist) + '\n')
				else:
					flip_support = 0
					keep_support = 0
					for c in ref_hg_directions:
						if not c in hg_directions: continue
						hg_weight = statistics.mean(hg_lengths[c])
						if hg_directions[c] == ref_hg_directions[c]:
							keep_support += hg_weight
						else:
							flip_support += hg_weight

					# flip the genbank visual if necessary, first BGC processed is used as reference guide
					if flip_support > keep_support:
						flip_printlist = []
						for gene_string in printlist:
							gene_info = gene_string.split('\t')
							new_forward = 'TRUE'
							if gene_info[4] == 'TRUE': new_forward = 'FALSE'
							new_gstart = int(last_gene_end) - int(gene_info[3])
							new_gend = int(last_gene_end) - int(gene_info[2])
							new_gene_string = '\t'.join([str(x) for x in
														 [gene_info[0], gene_info[1], new_gstart, new_gend, new_forward,
															gene_info[-2], gene_info[-1]]])
							flip_printlist.append(new_gene_string)
						gggenes_track_handle.write('\n'.join(flip_printlist) + '\n')
					else:
						gggenes_track_handle.write('\n'.join(printlist) + '\n')

			dummy_hg = None
			for bgc in bgc_hg_presence:
				for hg in hg_counts:
					dummy_hg = hg
					heatmap_track_handle.write('\t'.join([bgc, hg, bgc_hg_presence[bgc][hg], str(hg_counts[hg])]) + '\n')

			for bgc in all_bgcs_in_tree:
				if not bgc in bgc_gene_counts.keys():
					gggenes_track_handle.write('\t'.join([bgc] + ['NA']*4 + ['Absent', '"#FFFFFF"']) + '\n')
					heatmap_track_handle.write('\t'.join([bgc, dummy_hg, 'Absent', '1']) + '\n')

			gggenes_track_handle.close()
			heatmap_track_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Had difficulties creating tracks for visualization of BGC gene architecture along phylogeny using R libraries.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		rscript_plot_cmd = ["Rscript", RSCRIPT_FOR_BGSEE, phylogeny_file, gggenes_track_file, heatmap_track_file,
							result_pdf_file]
		if self.logObject:
			self.logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
		try:
			subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			self.logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))

		if self.logObject:
			self.logObject.info('Plotting completed (I think successfully)!')

	def constructCodonAlignments(self, outdir, cores=1, only_scc=False, list_alignments=False, filter_outliers=False):
		"""
		Function to automate construction of codon alignments. This function first extracts protein and nucleotide sequnces
		from BGC Genbanks, then creates protein alignments for each homolog group using MAFFT, and finally converts those
		into codon alignments using PAL2NAL.

		:param outdir: Path to output/workspace directory. Intermediate files (like extracted nucleotide and protein
						 sequences, protein and codon alignments, will be writen to respective subdirectories underneath this
						 one).
		:param cores: Number of cores/threads to use when fake-parallelizing jobs using multiprocessing.
		:param only_scc: Whether to construct codon alignments only for homolog groups which are found to be core and in
						 single copy for samples with the GCF. Note, if working with draft genomes and the BGC is fragmented
						 this should be able to still identify SCC homolog groups across the BGC instances belonging to the
						 GCF.
		"""

		nucl_seq_dir = os.path.abspath(outdir) + '/Nucleotide_Sequences/'
		prot_seq_dir = os.path.abspath(outdir) + '/Protein_Sequences/'
		prot_alg_dir = os.path.abspath(outdir) + '/Protein_Alignments/'
		codo_alg_dir = os.path.abspath(outdir) + '/Codon_Alignments/'
		if filter_outliers:
			nucl_seq_dir = os.path.abspath(outdir) + '/Nucleotide_Sequences_MAD_Refined/'
			prot_seq_dir = os.path.abspath(outdir) + '/Protein_Sequences_MAD_Refined/'
			prot_alg_dir = os.path.abspath(outdir) + '/Protein_Alignments_MAD_Refined/'
			codo_alg_dir = os.path.abspath(outdir) + '/Codon_Alignments_MAD_Refined/'

		if not os.path.isdir(nucl_seq_dir): os.system('mkdir %s' % nucl_seq_dir)
		if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
		if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
		if not os.path.isdir(codo_alg_dir): os.system('mkdir %s' % codo_alg_dir)

		all_samples = set(self.bgc_sample.values())
		try:
			inputs = []
			for hg in self.hg_genes:
				# if len(self.hg_genes[hg]) < 2: continue
				sample_counts = defaultdict(int)
				gene_sequences = {}
				for gene in self.hg_genes[hg]:
					gene_info = self.comp_gene_info[gene]
					bgc_id = gene_info['bgc_name']
					sample_id = self.bgc_sample[bgc_id]
					nucl_seq = gene_info['nucl_seq']
					prot_seq = gene_info['prot_seq']
					sample_counts[sample_id] += 1
					gid = sample_id + '|' + gene
					if only_scc:
						gid = sample_id
					gene_sequences[gid] = tuple([nucl_seq, prot_seq])
				samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
				# check that hg is single-copy-core
				if only_scc and len(samples_with_single_copy.symmetric_difference(all_samples)) > 0:
					continue
				elif only_scc and self.logObject:
					self.logObject.info('Homolog group %s detected as SCC across samples (not individual BGCs).' % hg)
				# check that hg is present in the original instances of GCF
				if len([x for x in gene_sequences.keys() if len(x.split('|')[1].split('_')[0]) == 3]) == 0: continue
				if filter_outliers:
					gene_sequences = util.determineOutliersByGeneLength(gene_sequences)
				inputs.append([hg, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, self.logObject])

			p = multiprocessing.Pool(cores)
			p.map(create_codon_msas, inputs)

			if not filter_outliers:
				self.nucl_seq_dir = nucl_seq_dir
				self.prot_seq_dir = prot_seq_dir
				self.prot_alg_dir = prot_alg_dir
				self.codo_alg_dir = codo_alg_dir
			else:
				self.nucl_filt_seq_dir = nucl_seq_dir
				self.prot_filt_seq_dir = prot_seq_dir
				self.prot_filt_alg_dir = prot_alg_dir
				self.codo_filt_alg_dir = codo_alg_dir

			if list_alignments:
				codon_alg_listings_file = outdir + 'Codon_Alignments_Listings.txt'
				if filter_outliers:
					codon_alg_listings_file = outdir + 'Codon_Alignments_Listings.MAD_Refined.txt'
				codon_alg_listings_handle = open(codon_alg_listings_file, 'w')
				for f in os.listdir(codo_alg_dir):
					codon_alg_listings_handle.write(f.split('.msa.fna')[0] + '\t' + codo_alg_dir + f + '\n')
				codon_alg_listings_handle.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with create protein/codon alignments of SCC homologs for BGC.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def constructGCFPhylogeny(self, output_alignment, output_phylogeny, only_scc=False):
		"""
		Function to create phylogeny based on codon alignments of SCC homolog groups for GCF.

		:param output_alignment: Path to output file for concatenated SCC homolog group alignment.
		:param output_phylogeny: Path to output file for approximate maximum-likelihood phylogeny produced by FastTree2 from
									 concatenated SCC homolog group alignment.
		"""
		try:
			if only_scc:
				bgc_sccs = defaultdict(lambda: "")
				fasta_data = []
				fasta_data_tr = []

				for f in os.listdir(self.codo_alg_dir):
					hg_align_msa = self.codo_alg_dir + f
					# concatenate gene alignments
					with open(hg_align_msa) as opm:
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

				scc_handle = open(output_alignment, 'w')

				for rec in zip(*fasta_data_tr):
					scc_handle.write(rec[0] + '\n' + ''.join(rec[1:]) + '\n')
				scc_handle.close()
			else:
				bgc_sccs = defaultdict(lambda: "")
				fasta_data = []
				fasta_data_tr = []

				for f in os.listdir(self.codo_alg_dir):
					hg_align_msa = self.codo_alg_dir + f
					#print(f)
					# perform consensus calling
					sample_seqs = defaultdict(list)
					with open(hg_align_msa) as opm:
						for rec in SeqIO.parse(opm, 'fasta'):
							sample = rec.id.split('|')[0]
							sample_seqs[sample].append(list(str(rec.seq).upper()))

					for samp in sample_seqs:
						samp_seqs = sample_seqs[samp]
						consensus_seq = []
						for alleles in zip(*samp_seqs):
							valid_alleles = set([a for a in list(alleles) if a in set(['A', 'C', 'G', 'T'])])
							if len(valid_alleles) == 1:
								consensus_seq.append(list(valid_alleles)[0])
							else:
								consensus_seq.append('-')
						bgc_sccs['>' + samp] += "".join(consensus_seq)

				for b in bgc_sccs:
					fasta_data.append([b] + list(bgc_sccs[b]))

				for i, ls in enumerate(zip(*fasta_data)):
					if i == 0:
						fasta_data_tr.append(ls)
					else:
						n_count = len([x for x in ls if x == '-'])
						if (float(n_count) / len(ls)) < 0.1:
							fasta_data_tr.append(list(ls))

				scc_handle = open(output_alignment, 'w')

				for rec in zip(*fasta_data_tr):
					scc_handle.write(rec[0] + '\n' + ''.join(rec[1:]) + '\n')
				scc_handle.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error('Had issues with creating concatenated alignment of the SCC homolog groups.')
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		# use FastTree2 to construct phylogeny
		fasttree_cmd = ['fasttree', '-nt', output_alignment, '>', output_phylogeny]
		if self.logObject:
			self.logObject.info('Running FastTree2 with the following command: %s' % ' '.join(fasttree_cmd))
		try:
			subprocess.call(' '.join(fasttree_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(fasttree_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(fasttree_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(fasttree_cmd))

	def refineBGCGenbanks(self, new_gcf_listing_file, outdir, first_boundary_homolog, second_boundary_homolog):
		"""
		Function to refine BGC Genbanks based on boundaries defined by two single copy core homolog groups. Genbanks
		are filtered to retain only features in between the positions of the two boundary homolog groups. Coordinates
		relevant to the BGC framework are updated (but not all location coordinates!!!).

		:param new_gcf_listing_file: Path to where new GCF listing file will be written.
		:param outdir: Path to workspace directory.
		:param first_boundary_homolog: Identifier of the first boundary homolog group
		:param second_boundary_homolog: Identifier of the second boundary homolog group
		"""
		try:
			refined_gbks_dir = outdir + 'Refined_Genbanks/'
			if not os.path.isdir(refined_gbks_dir): os.system('mkdir %s' % refined_gbks_dir)

			nglf_handle = open(new_gcf_listing_file, 'w')

			first_boundary_homolog_genes = self.hg_genes[first_boundary_homolog]
			second_boundary_homolog_genes = self.hg_genes[second_boundary_homolog]

			for bgc in self.pan_bgcs:
				bgc_genes = set(self.pan_bgcs[bgc].gene_information.keys())
				bgc_fbh_genes = bgc_genes.intersection(first_boundary_homolog_genes)
				bgc_sbh_genes = bgc_genes.intersection(second_boundary_homolog_genes)
				if len(bgc_fbh_genes) == 1 and len(bgc_sbh_genes) == 1:
					refined_gbk = refined_gbks_dir + self.bgc_gbk[bgc].split('/')[-1]
					self.pan_bgcs[bgc].refineGenbank(refined_gbk, list(bgc_fbh_genes)[0], list(bgc_sbh_genes)[0])
					nglf_handle.write('%s\t%s\n' % (self.bgc_sample[bgc], refined_gbk))
				elif self.logObject:
					self.logObject.warning(
						"Dropping the BGC genbank %s from consideration / refinement process because it does not have the boundary homolog groups in single-copy copy." %
						self.bgc_gbk[bgc])
			nglf_handle.close()

		except:
			if self.logObject:
				self.logObject.error('Had an issue refining BGC genbanks associated with GCF')
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def determineHgOrderIndex(self):
		"""
		Function to determine an "ordering" score for homolog groups in GCF. The order score is relative to each GCF,
		even a homolog group has a large or small order score indicates it is on the edges (beginning will be chosen
		arbitrarily).
		"""
		try:
			ref_hg_directions = {}
			bgc_gene_counts = defaultdict(int)
			for bgc in self.bgc_genes:
					bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])

			for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
				bgc = item[0]
				curr_bgc_genes = self.bgc_genes[bgc]
				hg_directions = {}
				hg_lengths = defaultdict(list)
				hg_starts = {}
				for g in curr_bgc_genes:
					ginfo = self.comp_gene_info[g]
					gstart = ginfo['start']
					gend = ginfo['end']
					if g in self.gene_to_hg:
						hg = self.gene_to_hg[g]
						hg_directions[hg] = ginfo['direction']
						hg_lengths[hg].append(gend - gstart)
						hg_starts[hg] = ginfo['start']

				reverse_flag = False
				if i == 0:
					ref_hg_directions = hg_directions
				else:
					flip_support = 0
					keep_support = 0
					for c in ref_hg_directions:
						if not c in hg_directions: continue
						hg_weight = statistics.mean(hg_lengths[c])
						if hg_directions[c] == ref_hg_directions[c]:
							keep_support += hg_weight
						else:
							flip_support += hg_weight

					# reverse ordering
					if flip_support > keep_support:
						reverse_flag = True
				for c in sorted(hg_starts.items(), key=itemgetter(1), reverse=reverse_flag):
						self.hg_order_scores[c[0]] += c[1]

		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues in attempting to calculate order score for each homolog group.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def runPopulationGeneticsAnalysis(self, outdir, cores=1, population=None, filter_outliers=False):
		"""
		Wrapper function which serves to parallelize population genetics analysis.

		:param outdir: The path to the workspace / output directory.
		:param cores: The number of cores (will be used for parallelizing)
		"""

		popgen_dir = outdir + 'Codon_PopGen_Analyses'
		plots_dir = outdir + 'Codon_MSA_Plots'
		final_output_file = outdir + 'Ortholog_Group_Information'

		if filter_outliers:
			final_output_file = outdir + 'Ortholog_Group_Information_MAD_Refined'
			popgen_dir = outdir + 'Codon_PopGen_Analyses_MAD_Refined'
			plots_dir = outdir + 'Codon_MSA_Plots_MAD_Refined'

		if population:
			final_output_file = final_output_file + '_Pop-' + str(population) + '.txt'
			popgen_dir += '_Pop-' + str(population) + '/'
			plots_dir += '_Pop-' + str(population) + '/'
		else:
			final_output_file = final_output_file + '.txt'
			popgen_dir += '/'
			plots_dir += '/'

		if not os.path.isdir(popgen_dir): os.system('mkdir %s' % popgen_dir)
		if not os.path.isdir(plots_dir): os.system('mkdir %s' % plots_dir)

		final_output_handle = open(final_output_file, 'w')
		header = ['gcf_id', 'homolog_group', 'annotation', 'hg_order_index', 'hg_median_copy_count', 'median_gene_length',
					'is_core_to_bgc', 'bgcs_with_hg', 'samples_with_hg', 'proportion_of_samples_with_hg', 'Tajimas_D', 'core_codons',
					'total_variable_codons', 'nonsynonymous_codons', 'synonymous_codons', 'dn_ds', 'all_domains']
		if population:
			header = ['population'] + header
		if self.bgc_population != None:
			header += ['populations_with_hg', 'population_proportion_of_members_with_hg', 'one_way_ANOVA_pvalues_sequence_similarity', 'one_way_ANOVA_pvalues_presence_absence']
		final_output_handle.write('\t'.join(header) + '\n')

		inputs = []
		input_codon_dir = self.codo_alg_dir
		if filter_outliers:
			input_codon_dir = self.codo_filt_alg_dir
		for f in os.listdir(input_codon_dir):
			hg = f.split('.msa.fna')[0]
			codon_alignment_fasta = input_codon_dir + f
			inputs.append([self.gcf_id, hg, codon_alignment_fasta, popgen_dir, plots_dir, self.comp_gene_info, self.hg_genes,
							 self.bgc_sample, self.hg_prop_multi_copy, self.hg_order_scores, dict(self.sample_population),
							 population, self.logObject])

		p = multiprocessing.Pool(cores)
		p.map(popgen_analysis_of_hg, inputs)

		final_output_handle = open(final_output_file, 'a+')
		for f in os.listdir(popgen_dir):
			if not f.endswith('_stats.txt'): continue
			with open(popgen_dir + f) as opf:
				for line in opf:
					line = line
					final_output_handle.write(line)
		final_output_handle.close()

	def identifyGCFInstances(self, outdir, sample_prokka_data, orthofinder_matrix_file, min_size=5, min_core_size=3,
							 gcf_to_gcf_transition_prob=0.9, background_to_background_transition_prob=0.9,
							 syntenic_correlation_threshold=0.8):
		"""
		Function to search for instances of GCF in sample using HMM based approach based on homolog groups as characters,
		"part of GCF" and "not part of GCF" as states - all trained on initial BGCs constituting GCF as identified by
		lsaBGC-Cluster.py. This function utilizes the convenient Python library Pomegranate.

		:param outdir: The path to the workspace / output directory.
		:param sample_prokka_data: Dictionary with keys being sample identifiers and values being paths to genbank and
									 proteome files from Prokka based annotation
		:param cores: The number of cores (will be used for parallelizing)
		"""

		# Estimate HMM parameters
		gcf_hg_probabilities = defaultdict(lambda: 0.0)
		other_hg_probabilities = defaultdict(lambda: 0.0)
		number_gcf_hgs = 0
		number_other_hgs = 0
		specific_hgs = set([])
		with open(orthofinder_matrix_file) as ofmf:
			for i, line in enumerate(ofmf):
				if i == 0: continue
				line = line.strip('\n')
				ls = line.split('\t')
				hg = ls[0]
				if hg in self.hg_genes.keys():
					number_gcf_hgs += 1
					total_genes = 0
					gcf_genes = 0
					for samp_genes in ls[1:]:
						for gene in samp_genes.split(', '):
							if not gene.strip(): continue
							total_genes += 1
							if gene in self.pan_genes:
								gcf_genes += 1

					if float(gcf_genes) == float(total_genes) and self.hg_max_self_evalue[hg][1] == True:
						specific_hgs.add(hg)

					other_hg_probabilities[hg] = 0.01
					gcf_hg_probabilities[hg] = 0.99

					if self.hg_max_self_evalue[hg][1] == False:
						other_hg_probabilities[hg] = min(1.0 - (gcf_genes/float(total_genes)), 0.2)
						gcf_hg_probabilities[hg] = 1.0 - other_hg_probabilities[hg]
				else:
					number_other_hgs += 1

		gcf_hg_probabilities['other'] = 0.2
		other_hg_probabilities['other'] = 0.8

		gcf_distribution = DiscreteDistribution(dict(gcf_hg_probabilities))
		other_distribution = DiscreteDistribution(dict(other_hg_probabilities))

		print(gcf_distribution)
		print(other_distribution)
		print(self.core_homologs)
		print(self.protocluster_core_homologs)
		print(specific_hgs)
		gcf_state = State(gcf_distribution, name='GCF')
		other_state = State(other_distribution, name='Non GCF')

		model = HiddenMarkovModel()
		model.add_states(gcf_state, other_state)

		# estimate transition probabilities
		gcf_to_gcf = gcf_to_gcf_transition_prob #float(number_gcf_hgs - 1) / float(number_gcf_hgs)
		gcf_to_other = 1.0 - gcf_to_gcf_transition_prob #1.0 - gcf_to_gcf
		other_to_other = background_to_background_transition_prob #float(number_other_hgs - 1) / float(number_other_hgs)
		other_to_gcf = 1.0 - background_to_background_transition_prob #1.0 - other_to_other

		start_to_gcf = 0.5 # float(number_gcf_hgs)/float(number_gcf_hgs + number_other_hgs)
		start_to_other = 0.5 #1.0 - start_to_gcf
		gcf_to_end = 0.5 # float(number_gcf_hgs)/float(number_gcf_hgs + number_other_hgs)
		other_to_end = 0.5 # 1.0 - gcf_to_end

		model.add_transition(model.start, gcf_state, start_to_gcf)
		model.add_transition(model.start, other_state, start_to_other)
		model.add_transition(gcf_state, model.end, gcf_to_end)
		model.add_transition(other_state, model.end, other_to_end)
		model.add_transition(gcf_state, gcf_state, gcf_to_gcf)
		model.add_transition(gcf_state, other_state, gcf_to_other)
		model.add_transition(other_state, gcf_state, other_to_gcf)
		model.add_transition(other_state, other_state, other_to_other)

		model.bake()

		# open handle to file where expanded GCF listings will be written
		expanded_gcf_list_file = outdir + 'GCF_Expanded.txt'
		expanded_gcf_list_handle = open(expanded_gcf_list_file, 'w')
		all_samples = set([])
		with open(self.bgc_genbanks_listing) as obglf:
			for line in obglf:
				expanded_gcf_list_handle.write(line)

		sample_bgc_ids = defaultdict(lambda: 1)

		bgc_genbanks_dir = os.path.abspath(outdir + 'BGC_Genbanks') + '/'
		if not os.path.isdir(bgc_genbanks_dir): os.system('mkdir %s' % bgc_genbanks_dir)

		bgc_hmm_evalues_file = outdir + 'GCF_NewInstances_HMMEvalues.txt'
		bgc_hmm_evalues_handle = open(bgc_hmm_evalues_file, 'w')

		# start process of finding new BGCs
		sample_lt_to_hg = defaultdict(dict)
		sample_hgs = defaultdict(set)
		sample_protein_to_hg = defaultdict(dict)
		sample_lt_to_evalue = defaultdict(dict)
		for lt in self.hmmscan_results:
			for i, hits in enumerate(sorted(self.hmmscan_results[lt], key=itemgetter(1))):
				if i == 0:
					sample_lt_to_hg[hits[2]][lt] = hits[0]
					sample_hgs[hits[2]].add(hits[0])
					sample_protein_to_hg[hits[2]][lt] = hits[0]
					sample_lt_to_evalue[hits[2]][lt] = hits[1]

		sample_hg_proteins = defaultdict(lambda: defaultdict(set))
		for sample in sample_hgs:
			if len(sample_hgs[sample]) < 3: continue
			print(sample)

			sample_gcf_predictions = []
			for scaffold in self.scaffold_genes[sample]:
				lts_with_start = []
				for lt in self.scaffold_genes[sample][scaffold]:
					lts_with_start.append([lt, self.gene_location[sample][lt]['start']])

				hgs_ordered = []
				lts_ordered = []
				for lt_start in sorted(lts_with_start, key=itemgetter(1)):
					lt, start = lt_start
					lts_ordered.append(lt)
					if lt in sample_lt_to_hg[sample]:
						hgs_ordered.append(sample_lt_to_hg[sample][lt])
					else:
						hgs_ordered.append('other')

				hg_seq = numpy.array(list(hgs_ordered))
				hmm_predictions = model.predict(hg_seq)

				gcf_state_lts = []
				gcf_state_hgs = []
				for i, hg_state in enumerate(hmm_predictions):
					lt = lts_ordered[i]
					hg = hgs_ordered[i]
					if hg_state == 0:
						gcf_state_lts.append(lt)
						gcf_state_hgs.append(hg)
					if hg_state == 1 or i == (len(hmm_predictions)-1):
						print(gcf_state_lts)
						if len(set(gcf_state_hgs).difference("other")) >= 3:
							print(gcf_state_lts)
							boundary_lt_featured = False
							features_specific_hg = False
							features_protocoluster_hg = False
							if len(self.protocluster_core_homologs.intersection(set(gcf_state_hgs).difference('other'))) > 0: features_protocoluster_hg = True
							if len(self.boundary_genes[sample].intersection(set(gcf_state_lts).difference('other'))) > 0: boundary_lt_featured = True
							if len(specific_hgs.intersection(set(gcf_state_hgs).difference('other'))) > 0: features_specific_hg = True
							sample_gcf_predictions.append([gcf_state_lts, gcf_state_hgs, len(gcf_state_lts), len(set(gcf_state_hgs).difference("other")), len(set(gcf_state_hgs).difference("other").intersection(self.core_homologs)), scaffold, boundary_lt_featured, features_specific_hg, features_protocoluster_hg])
							print(sample_gcf_predictions[-1])
						gcf_state_lts = []
						gcf_state_hgs = []

			if len(sample_gcf_predictions) == 0: continue

			sorted_sample_gcf_predictions = [x for x in sorted(sample_gcf_predictions, key=itemgetter(3), reverse=True)]

			sample_gcf_predictions_filtered = []
			sample_edge_gcf_predictions_filtered = []
			cumulative_edge_hgs = set([])
			visited_scaffolds_with_edge_gcf_segment = set([])
			for gcf_segment in sorted_sample_gcf_predictions:
				if (gcf_segment[3] >= min_size and gcf_segment[4] >= min_core_size) or (gcf_segment[-1]) or (gcf_segment[-2]) or (gcf_segment[3] >= 3 and gcf_segment[-3] and not gcf_segment[5] in visited_scaffolds_with_edge_gcf_segment):
					# code to determine whether syntenically, the considered segment aligns with what is expected.
					segment_hg_order = []
					bgc_hg_orders = defaultdict(list)

					copy_count_of_hgs_in_segment = defaultdict(int)
					for hg in gcf_segment[1]:
						copy_count_of_hgs_in_segment[hg] += 1

					for gi, g in enumerate(gcf_segment[0]):
						hg = gcf_segment[1][gi]
						if copy_count_of_hgs_in_segment[hg] != 1: continue
						gene_midpoint = (self.gene_location[sample][g]['start'] + self.gene_location[sample][g]['end'])/2.0
						segment_hg_order.append(gene_midpoint)

						for bgc in self.bgc_genes:
							bg_matching = []
							for bg in self.bgc_genes[bgc]:
								if bg in self.gene_to_hg:
									hg_of_bg = self.gene_to_hg[bg]
									if hg_of_bg == hg: bg_matching.append(bg)
							if len(bg_matching) == 1:
								bgc_gene_midpoint = (self.comp_gene_info[bg_matching[0]]['start'] + self.comp_gene_info[bg_matching[0]]['end']) / 2.0
								bgc_hg_orders[bgc].append(bgc_gene_midpoint)
							else:
								bgc_hg_orders[bgc].append(None)

					best_corr = None
					for bgc in self.bgc_genes:
						try:
							assert(len(segment_hg_order) == len(bgc_hg_orders[bgc]))

							list1 = []
							list2 = []
							for iter, hgval1 in enumerate(segment_hg_order):
								hgval2 = bgc_hg_orders[bgc][iter]
								if hgval1 == None or hgval2 == None: continue
								list1.append(hgval1)
								list2.append(hgval2)

							corr, pval = pearsonr(list1, list2)
							corr = abs(corr)
							print(corr)
							if (best_corr and best_corr < corr) or (not best_corr):
								best_corr = corr
						except:
							pass

					if not best_corr or best_corr < syntenic_correlation_threshold: continue

					if (gcf_segment[3] >= min_size and gcf_segment[4] >= min_core_size) or (gcf_segment[-1]) or (gcf_segment[-2]):
						sample_gcf_predictions_filtered.append(gcf_segment)
					elif gcf_segment[3] >= 3 and gcf_segment[-3] and not gcf_segment[5] in visited_scaffolds_with_edge_gcf_segment:
						sample_edge_gcf_predictions_filtered.append(gcf_segment)
						visited_scaffolds_with_edge_gcf_segment.add(gcf_segment[5])
						cumulative_edge_hgs = cumulative_edge_hgs.union(set(gcf_segment[1]))

			if len(sample_edge_gcf_predictions_filtered) >= 2:
				print(sample_edge_gcf_predictions_filtered)
				if len(cumulative_edge_hgs) >= min_size and len(cumulative_edge_hgs.intersection(self.core_homologs)) >= min_core_size:
					sample_gcf_predictions_filtered += sample_edge_gcf_predictions_filtered

			print(sample)
			for gcf_state_lts in sample_gcf_predictions_filtered:
				print(gcf_state_lts)
				clean_sample_name = util.cleanUpSampleName(sample)
				bgc_genbank_file = bgc_genbanks_dir + clean_sample_name + '_BGC-' + str(sample_bgc_ids[sample]) + '.gbk'
				sample_bgc_ids[sample] += 1

				min_bgc_pos = min([self.gene_location[sample][g]['start'] for g in gcf_state_lts[0]])
				max_bgc_pos = max([self.gene_location[sample][g]['end'] for g in gcf_state_lts[0]])

				util.createBGCGenbank(sample_prokka_data[sample]['genbank'], bgc_genbank_file, gcf_state_lts[5], min_bgc_pos, max_bgc_pos)
				expanded_gcf_list_handle.write('\t'.join([clean_sample_name, bgc_genbank_file]) + '\n')

				for i, lt in enumerate(gcf_state_lts[0]):
					hg = gcf_state_lts[1][i]
					evalue = 10.0
					if lt in sample_lt_to_evalue[sample]: evalue = sample_lt_to_evalue[sample][lt]
					bgc_hmm_evalues_handle.write('\t'.join([bgc_genbank_file, sample, lt, hg, str(evalue)]) + '\n')

				for lt in gcf_state_lts[0]:
					if lt in sample_protein_to_hg[sample].keys():
						hg = sample_protein_to_hg[sample][lt]
						sample_hg_proteins[clean_sample_name][hg].add(lt)
				all_samples.add(clean_sample_name)

		expanded_gcf_list_handle.close()
		bgc_hmm_evalues_handle.close()

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
						sample_hg_proteins[original_samples[j]][hg] = sample_hg_proteins[original_samples[j]][hg].union(set(prot.split(', ')))

		expanded_orthofinder_matrix_file = outdir + 'Orthogroups.expanded.csv'
		expanded_orthofinder_matrix_handle = open(expanded_orthofinder_matrix_file, 'w')

		header = [''] + [s for s in sorted(all_samples)]
		expanded_orthofinder_matrix_handle.write('\t'.join(header) + '\n')
		for hg in sorted(all_hgs):
			printlist = [hg]
			for s in sorted(all_samples):
				printlist.append(', '.join(sample_hg_proteins[s][hg]))
			expanded_orthofinder_matrix_handle.write('\t'.join(printlist) + '\n')
		expanded_orthofinder_matrix_handle.close()

	def extractGenesAndCluster(self, genes_representative_fasta, genes_fasta, codon_alignments_file, bowtie2_db_prefix):
		"""
		Function to cluster gene sequences for each homolog group using codon alignments and then select representative
		sequences which will be used to build a Bowtie2 reference database.

		:param genes_representative_fasta: Path to FASTA file which will harbor gene + flanks
		:param codon_alignments_file: File listing paths to codon alignments (2nd column) for each homolog group (1st
		                              column).
		:param bowtei2_db_prefix: Path to prefix of Bowtie2 refernece database/index to be used for aligning downstream
															in the lsaBGC-DiscoVary workflow
		"""

		try:
			grf_handle = open(genes_representative_fasta, 'w')
			gf_handle = open(genes_fasta, 'w')
			with open(codon_alignments_file) as ocaf:
				for line in ocaf:
					line = line.strip()
					hg, cod_alignment_fasta = line.split('\t')
					alleles_clustered, pair_matching = util.determineAllelesFromCodonAlignment(cod_alignment_fasta)
					all_genes_in_codon_alignments = set([item.split('|')[-1] for sublist in alleles_clustered for item in sublist])
					try:
						assert(len(self.hg_genes[hg].symmetric_difference(all_genes_in_codon_alignments)) == 0)
					except:
						self.logObject.warning("Not all genes featured in codon alignment for homolog group %s, these will be excluded." % hg)

					for allele_cluster in alleles_clustered:
						best_rep_score = defaultdict(float)
						for ag1 in alleles_clustered[allele_cluster]:
							gf_handle.write('>' + hg + '|' + allele_cluster + '|' + ag1 + '\n' + str(self.comp_gene_info[ag1.split('|')[-1]]['nucl_seq']) + '\n')
							for ag2 in alleles_clustered[allele_cluster]:
								best_rep_score[ag1] += pair_matching[ag1][ag2]
						representative_gene = sorted([x for x in best_rep_score if best_rep_score[x] == max(best_rep_score.values())])[0]
						for ag in alleles_clustered[allele_cluster]:
							self.instance_to_haplotype[hg + '|' + allele_cluster + '|' + ag] = hg + '|' + allele_cluster + '|' + representative_gene
						grf_handle.write('>' + hg + '|' + allele_cluster + '|' + representative_gene + '\n' + str(self.comp_gene_info[representative_gene.split('|')[-1]]['nucl_seq']) + '\n')
			grf_handle.close()
			gf_handle.close()

			bowtie2_build = ['bowtie2-build', genes_representative_fasta, bowtie2_db_prefix]
			if self.logObject:
				self.logObject.info('Running the following command: %s' % ' '.join(bowtie2_build))
			try:
				subprocess.call(' '.join(bowtie2_build), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				if self.logObject:
					self.logObject.info('Successfully ran: %s' % ' '.join(bowtie2_build))
			except:
				if self.logObject:
					self.logObject.error('Had an issue running: %s' % ' '.join(bowtie2_build))
				raise RuntimeError('Had an issue running: %s' % ' '.join(bowtie2_build))
			if self.logObject:
				self.logObject.info('Build Bowtie2 database/index for %s' % genes_representative_fasta)

		except Exception as e:
			if self.logObject:
				self.logObject.error("Unable to create FASTA of genes and/or subsequent bowtie2 database.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def extractGeneWithFlanksAndCluster(self, genes_with_flanks_fasta, cd_hit_clusters_fasta_file, cd_hit_nr_fasta_file, bowtie2_db_prefix):
		"""
		Function to extract gene sequences and surrounding flanking sequences into a FASTA file, which will then be
		clustered using CD-HIT at both a coarse and very granular (just remove 100% redundancy) and to construct
		a Bowtie2 reference database of the granular clustering. From the coarse clustering, representative genes
		will be selected to depict different alleles of the gene and stored as a dictionary.

		:param genes_with_flanks_fasta: Path to FASTA file which will harbor gene + flanks
		:param cd_hit_clusters_fasta_file: Path to FASTA file which will be used to output coarse level clustering by
											 CD-HIT
		:param cd_hit_nr_fasta_file: Path to FASTA file which will be used to output granular level clustering by CD-HIT
		:param bowtei2_db_prefix: Path to prefix of Bowtie2 refernece database/index to be used for aligning downstream
															in the lsaBGC-DiscoVary workflow
		"""
		try:
			gwff_handle = open(genes_with_flanks_fasta, 'w')
			for bgc in self.bgc_genes:
				for gene in self.bgc_genes[bgc]:
					if gene in self.gene_to_hg:
						gwff_handle.write('>' + gene + '|' + bgc + '|' + self.gene_to_hg[gene] + '\n' + self.comp_gene_info[gene]['nucl_seq_with_flanks'] + '\n')
			gwff_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Unable to extract flanking sequences of gene into FASTA file.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		cd_hit_nr = ['cd-hit-est', '-i', genes_with_flanks_fasta, '-o', cd_hit_nr_fasta_file, '-G', '1', '-g',
					 '1', '-d', '0', '-n', '10', '-M', '2000', '-c', '1.0', '-aL', '1.0', '-aS', '1.0', '-T', '1']
		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(cd_hit_nr))
		try:
			subprocess.call(' '.join(cd_hit_nr), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(cd_hit_nr))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(cd_hit_nr))
			raise RuntimeError('Had an issue running: %s' % ' '.join(cd_hit_nr))
		if self.logObject:
			self.logObject.info('Ran CD-HIT for collapsing redundancy.')

		cd_hit_cluster = ['cd-hit-est', '-i', genes_with_flanks_fasta, '-o', cd_hit_clusters_fasta_file, '-G', '1',
							'-g',
							'1', '-d', '0', '-n', '10', '-M', '2000', '-c', '0.98', '-aL', '0.95', '-aS', '0.95', '-T',
							'1']
		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(cd_hit_cluster))
		try:
			subprocess.call(' '.join(cd_hit_cluster), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(cd_hit_cluster))
		except:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(cd_hit_cluster))
			raise RuntimeError('Had an issue running: %s' % ' '.join(cd_hit_cluster))
		if self.logObject:
			self.logObject.info('Ran CD-HIT for clustering genes, with their flanks, into haplotype groups.')

		bowtie2_build = ['bowtie2-build', cd_hit_nr_fasta_file, bowtie2_db_prefix]
		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(bowtie2_build))
		try:
			subprocess.call(' '.join(bowtie2_build), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(bowtie2_build))
		except:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(bowtie2_build))
			raise RuntimeError('Had an issue running: %s' % ' '.join(bowtie2_build))
		if self.logObject:
			self.logObject.info('Build Bowtie2 database/index for %s' % cd_hit_nr_fasta_file)

		try:
			cd_hit_clusters_cltr_file = cd_hit_clusters_fasta_file + '.clstr'
			assert (os.path.isfile(cd_hit_clusters_cltr_file))

			cluster = []
			with open(cd_hit_clusters_cltr_file) as off:
				for line in off:
					line = line.strip()
					ls = line.split()
					if line.startswith('>'):
						if len(cluster) > 0:
							for g in cluster:
								self.instance_to_haplotype[g] = rep
						cluster = []
						rep = None
					else:
						gene_id = ls[2][1:-3]
						cluster.append(gene_id)
						if line.endswith('*'): rep = gene_id
			if len(cluster) > 0:
				if len(cluster) > 0:
					for g in cluster:
						self.instance_to_haplotype[g] = rep

		except Exception as e:
			if self.logObject:
				self.logObject.error("Unable to parse CD-HIT clustering of gene sequences (with flanks) to obtain representative sequence per cluster.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def runSNVMining(self, paired_end_sequencing_file, bowtie2_ref_fasta, codon_alignment_file, bowtie2_alignment_dir, results_dir, cores=1):
		"""
		Wrapper function for mining for novel SNVs across genes of GCF.

		:param paired_end_sequencing_file: tab delimited file with three columns: (1) sample name (2) path to forward
											 reads and (3) path to reverse reads
		:param bowtie2_ref_fasta: FASTA file corresponding
		:param bowtie2_alignment_dir: Path to directory where Bowtie2 alignments were written. This directory should
										include BAM files ending in *.filtered.sorted.bam which are sorted and indexed.
		:param results_dir: Path to directory where results of SNV mining will be written.
		:param cores: The number of processes to be run in parallel.
		"""

		try:

			gene_pos_to_msa_pos = defaultdict(lambda: defaultdict(dict))
			gene_pos_to_allele = defaultdict(lambda: defaultdict(dict))
			codon_alignment_lengths = defaultdict(int)
			with open(codon_alignment_file) as ocaf:
				for line in ocaf:
					line = line.strip()
					hg, cod_alignment = line.split('\t')
					seq_count = 0
					with open(cod_alignment) as oca:
						for j, rec in enumerate(SeqIO.parse(oca, 'fasta')):
							sample_id, gene_id = rec.id.split('|')
							real_pos = 1
							for msa_pos, bp in enumerate(str(rec.seq)):
								if j == 0:
									codon_alignment_lengths[hg] += 1
								if bp != '-':
									gene_pos_to_allele[hg][gene_id][real_pos] = bp.upper()
									gene_pos_to_msa_pos[hg][gene_id][real_pos] = msa_pos+1
									real_pos += 1

							seq_count += 1

			process_args = []
			with open(paired_end_sequencing_file) as opesf:
				for line in opesf:
					line = line.strip()
					sample, frw_read, rev_read = line.split('\t')
					process_args.append([sample, bowtie2_alignment_dir + sample + '.sorted.bam',
										 bowtie2_ref_fasta, self.instance_to_haplotype, results_dir, self.hg_genes,
										 self.comp_gene_info, dict(gene_pos_to_msa_pos), dict(gene_pos_to_allele),
										 dict(codon_alignment_lengths), self.logObject])

			p = multiprocessing.Pool(cores)
			p.map(snv_miner_single, process_args)
			p.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def phaseAndSummarize(self, paired_end_sequencing_file, codon_alignment_file, snv_mining_outdir, outdir):
		try:
			novelty_report_file = outdir + 'Novelty_Report.txt'
			homolog_presence_report_file = outdir + 'Sample_Homolog_Group_Presence.txt'
			sample_bgc_coverage_file = outdir + 'Sample_BGC_Overall_Depth.txt'
			no_handle = open(novelty_report_file, 'w')
			hpr_handle = open(homolog_presence_report_file, 'w')
			sbc_handle = open(sample_bgc_coverage_file, 'w')
			no_handle.write('\t'.join(['gcf_id', 'sample', 'homolog_group', 'position_along_msa', 'alternate_allele',
									   'codon_position', 'alternate_codon', 'alternate_aa', 'dn_or_ds', 'ts_or_tv',
									   'reference_allele', 'reference_sample', 'reference_gene', 'reference_position',
									   'ref_codon', 'ref_aa', 'snv_support']) + '\n')
			mges = set(['transp', 'integrase'])
			purine_alleles = set(['A', 'G'])
			specific_homolog_groups = set([])
			for hg in self.hg_differentiation_stats:
				if self.hg_differentiation_stats[hg]['able_to_differentiate']:
					specific_homolog_groups.add(hg)

			gene_ignore_positions = defaultdict(set)
			gene_edgy_positions = defaultdict(set)
			gene_exclusion = defaultdict(lambda: [0, 1000000])
			gene_core_positions = defaultdict(set)
			gene_pos_to_msa_pos = defaultdict(lambda: defaultdict(dict))
			gene_pos_to_allele = defaultdict(lambda: defaultdict(dict))
			msa_pos_alleles = defaultdict(lambda: defaultdict(set))
			msa_pos_ambiguous_freqs = defaultdict(lambda: defaultdict(float))
			with open(codon_alignment_file) as ocaf:
				for line in ocaf:
					line = line.strip()
					hg, cod_alignment = line.split('\t')
					seq_count = 0
					msa_pos_ambiguous_counts = defaultdict(int)
					msa_pos_non_ambiguous_counts = defaultdict(int)
					seqlen_information = {}
					msa_positions = set([])
					with open(cod_alignment) as oca:
						for rec in SeqIO.parse(oca, 'fasta'):
							sequence_without_gaps = str(rec.seq).upper().replace('-', '')
							sample_id, gene_id = rec.id.split('|')
							seqlen = len(sequence_without_gaps)
							seqlen_lower = 50
							seqlen_upper = seqlen - seqlen_lower
							seqlen_information[gene_id] = [seqlen, seqlen_lower, seqlen_upper]

							real_pos = 1
							for msa_pos, bp in enumerate(str(rec.seq)):
								msa_positions.add(msa_pos+1)
								if bp != '-':
									msa_pos_non_ambiguous_counts[msa_pos + 1] += 1
									gene_pos_to_allele[hg][gene_id][real_pos] = bp.upper()
									gene_pos_to_msa_pos[hg][gene_id][real_pos] = msa_pos + 1
									real_pos += 1
									msa_pos_alleles[hg][msa_pos + 1].add(bp.upper())
								else:
									msa_pos_ambiguous_counts[msa_pos + 1] += 1
							seq_count += 1

					for pos, seqs_with_al in msa_pos_non_ambiguous_counts.items():
						if seqs_with_al/float(seq_count) >= 0.8:
							gene_core_positions[hg].add(pos)

					cod_algn_len = max(msa_positions)
					for p in (list(range(1, 51)) + list(range(cod_algn_len-50, cod_algn_len+1))):
						gene_ignore_positions[hg].add(p)

					for pos in msa_pos_ambiguous_counts:
						msa_pos_ambiguous_freqs[hg][pos] = msa_pos_ambiguous_counts[pos] / float(seq_count)
						if (msa_pos_ambiguous_counts[pos]/float(seq_count)) >= 0.1:
							for p in range(pos - 50, pos + 51):
								gene_ignore_positions[hg].add(p)

			with open(paired_end_sequencing_file) as ossf:
				for line in ossf:
					pe_sample = line.strip().split('\t')[0]
					result_file = snv_mining_outdir + pe_sample + '.txt'
					snv_file = snv_mining_outdir + pe_sample + '.snvs'
					homolog_group_positions = defaultdict(list)
					homolog_group_depths = defaultdict(list)
					previous_positions = []
					hg_median_depths = {}
					hg_first_position_of_stop_codon = defaultdict(lambda: None)
					with open(result_file) as orf:
						for i, line in enumerate(orf):
							if i == 0: continue
							line = line.strip()
							homolog_group = line.split(',')[0]
							position, a, c, g, t = [int(x) for x in line.split(',')[1:]]
							min_stop_codon_depth = 5
							total_depth = sum([a,c,g,t])

							major_alleles = []
							if a >= min_stop_codon_depth:
								major_alleles.append('A')
							if c >= min_stop_codon_depth:
								major_alleles.append('C')
							if g >= min_stop_codon_depth:
								major_alleles.append('G')
							if t >= min_stop_codon_depth:
								major_alleles.append('T')

							if position % 3 == 0:
								if len(previous_positions) == 2:
									for al1 in previous_positions[0]:
										for al2 in previous_positions[1]:
											for al3 in major_alleles:
												codon = al1 + al2 + al3
												if codon in set(['TAG', 'TGA', 'TAA']):
													if not homolog_group in hg_first_position_of_stop_codon:
														hg_first_position_of_stop_codon[homolog_group] = position-2
								previous_positions = []
							else:
								previous_positions.append(major_alleles)

							homolog_group_depths[homolog_group].append(total_depth)
							homolog_group_positions[homolog_group].append(position)

					present_homolog_groups = set([])
					mid_depths_all_present_hgs = []
					for hg in homolog_group_depths:
						hg_depths = homolog_group_depths[hg]
						hg_positions = homolog_group_positions[hg]

						core_positions_covered = 0
						for i, pos in enumerate(hg_positions):
							dp = hg_depths[i]
							if pos in gene_core_positions[hg] and dp >= 1:
								core_positions_covered += 1

						if float(core_positions_covered)/len(gene_core_positions[hg]) >= 0.9:
							present_homolog_groups.add(hg)

							hg10per = int(np.percentile(list(range(0, len(hg_depths) + 1)), 10))
							hg_depths_trimmed = hg_depths[hg10per:-hg10per]
							hg_trimmed_depth_median = statistics.median(hg_depths_trimmed)
							hg_median_depths[hg] = hg_trimmed_depth_median

					if len(hg_median_depths) < 5: continue
					median_of_medians = statistics.median(list(hg_median_depths.values()))
					mad_of_medians = median_absolute_deviation(list(hg_median_depths.values()))

					outlier_homolog_groups = set([])
					for hg in present_homolog_groups:
						hg_median = hg_median_depths[hg]
						if hg_median > (median_of_medians + (2*mad_of_medians)) or hg_median < (median_of_medians - (2*mad_of_medians)):
							outlier_homolog_groups.add(hg)

					present_homolog_groups = present_homolog_groups.difference(outlier_homolog_groups)
					if len(present_homolog_groups.intersection(self.core_homologs))/float(len(self.core_homologs)) < 0.7 and len(present_homolog_groups.intersection(specific_homolog_groups)) == 0: continue

					for hg in present_homolog_groups:
						hg_depths = homolog_group_depths[hg]
						hg10per = int(np.percentile(list(range(0, len(hg_depths) + 1)), 10))
						hg_depths_trimmed = hg_depths[hg10per:-hg10per]
						mid_depths_all_present_hgs += hg_depths_trimmed
						product_has_mge_term = False
						for gene in self.hg_genes[hg]:
							for word in mges:
								if word in self.comp_gene_info[gene]['product'].lower():
									product_has_mge_term = True
						hpr_handle.write('\t'.join([str(x) for x in [pe_sample, hg, (hg in specific_homolog_groups),
																	 self.hg_prop_multi_copy[hg], product_has_mge_term,
																	 statistics.median(hg_depths_trimmed),
																	 hg_first_position_of_stop_codon[hg],
																	 ','.join([str(x) for x in sorted(gene_ignore_positions[hg])])]]) + '\n')

					trimmed_depth_median = statistics.median(mid_depths_all_present_hgs)
					trimmed_depth_mad = median_absolute_deviation(mid_depths_all_present_hgs)
					sbc_handle.write(pe_sample + '\t' + str(trimmed_depth_median) + '\t' + str(trimmed_depth_mad) + '\n')

					with open(snv_file) as of:
						for i, line in enumerate(of):
							line = line.strip()
							ls = line.split('\t')
							snv_count = ls[1]
							gsh, ref_pos, ref_al, alt_al = ls[0].split('_|_')
							hg, allele_cluster, sample, gene = gsh.split('|')

							ref_pos = int(ref_pos)

							# check whether to regard this SNV as potentially legit
							if not hg in present_homolog_groups: continue
							if self.hg_prop_multi_copy[hg] >= 0.05: continue
							if any(word in self.comp_gene_info[gene]['product'].lower() for word in mges): continue
							if not ref_pos in gene_pos_to_msa_pos[hg][gene]: continue
							msa_pos = gene_pos_to_msa_pos[hg][gene][ref_pos]
							msa_pos_als = msa_pos_alleles[hg][msa_pos]
							if hg_first_position_of_stop_codon[hg] != None and msa_pos >= hg_first_position_of_stop_codon[hg]: continue
							if int(snv_count) >= 5 and int(homolog_group_depths[hg][msa_pos-1]) <= (trimmed_depth_median+(2*trimmed_depth_mad)) and int(homolog_group_depths[hg][msa_pos-1]) >= trimmed_depth_median-(3*trimmed_depth_mad):
								assert (ref_al in msa_pos_alleles[hg][msa_pos] and ref_al == gene_pos_to_allele[hg][gene][ref_pos])
								if not alt_al in msa_pos_als and not msa_pos in gene_ignore_positions[hg]:
									#if not alt_al in msa_pos_als and not msa_pos in gene_edgy_positions[hg] and msa_pos_ambiguous_freqs[hg][msa_pos] <= 0.1:
									codon_position = None
									ref_codon = None
									ref_aa = None
									alt_codon = None
									alt_aa = None
									dn_or_ds = None
									ts_or_tv = "transition"
									if (ref_al in purine_alleles) != (alt_al in purine_alleles):
										ts_or_tv = "transversion"
									if ref_pos % 3 == 1:
										codon_position = 1
										ref_codon = ref_al + gene_pos_to_allele[hg][gene][ref_pos + 1] + \
													gene_pos_to_allele[hg][gene][ref_pos + 2]
										alt_codon = alt_al + gene_pos_to_allele[hg][gene][ref_pos + 1] + \
													gene_pos_to_allele[hg][gene][ref_pos + 2]
										ref_aa = str(Seq(ref_codon).translate())
										alt_aa = str(Seq(alt_codon).translate())
									elif ref_pos % 3 == 2:
										codon_position = 2
										ref_codon = gene_pos_to_allele[hg][gene][ref_pos - 1] + ref_al + \
													gene_pos_to_allele[hg][gene][ref_pos + 1]
										alt_codon = gene_pos_to_allele[hg][gene][ref_pos - 1] + alt_al + \
													gene_pos_to_allele[hg][gene][ref_pos + 1]
										ref_aa = str(Seq(ref_codon).translate())
										alt_aa = str(Seq(alt_codon).translate())
									elif ref_pos % 3 == 0:
										codon_position = 3
										ref_codon = gene_pos_to_allele[hg][gene][ref_pos - 2] + \
													gene_pos_to_allele[hg][gene][ref_pos - 1] + ref_al
										alt_codon = gene_pos_to_allele[hg][gene][ref_pos - 2] + \
													gene_pos_to_allele[hg][gene][ref_pos - 1] + alt_al
										ref_aa = str(Seq(ref_codon).translate())
										alt_aa = str(Seq(alt_codon).translate())
									if ref_aa != alt_aa:
										dn_or_ds = "non-synonymous"
									else:
										dn_or_ds = "synonymous"
									no_handle.write('\t'.join([str(x) for x in [self.gcf_id, pe_sample, hg, msa_pos, alt_al,
																		codon_position, alt_codon, alt_aa, dn_or_ds,
																		ts_or_tv, ref_al, sample, gene, ref_pos,
																		ref_codon, ref_aa, snv_count]]) + '\n')


			no_handle.close()
			hpr_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with generating matrices showcasing allele presence across samples.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

def snv_miner_single(input_args):
	"""
	Function to mine for novel SNVs and identify alleles of homolog groups for GCF present in paired-end sequencing
	dataset.
	"""
	sample, bam_alignment, ref_fasta, hg_gene_to_rep, res_dir, bgc_hg_genes, comp_gene_info, gene_pos_to_msa_pos, gene_pos_to_allele, codon_alignment_lengths, logObject = input_args
	try:
		hg_rep_genes = defaultdict(set)
		for g, r in hg_gene_to_rep.items():
			hg_rep_genes[r].add(g)

		if not os.path.isfile(bam_alignment): return
		snvs_file = res_dir + sample + '.snvs'
		result_file = res_dir + sample + '.txt'
		details_file = res_dir + sample + '.full.txt'
		snv_outf = open(snvs_file, 'w')
		res_outf = open(result_file, 'w')
		det_outf = open(details_file, 'w')

		res_outf.write('Contig,Position,Sample-A,Sample-C,Sample-G,Sample-T\n')

		bam_handle = pysam.AlignmentFile(bam_alignment, 'rb')

		topaligns_file = res_dir + sample + '_topaligns.bam'
		topaligns_file_sorted = res_dir + sample + '_topaligns.sorted.bam'
		topaligns_handle = pysam.AlignmentFile(topaligns_file, "wb", template=bam_handle)

		for hg, hg_genes in bgc_hg_genes.items():
			read_ascores_per_allele = defaultdict(list)
			hg_genes_covered = 0
			gene_sequence = {}
			total_reads = set([])
			with open(ref_fasta) as opff:
				for rec in SeqIO.parse(opff, 'fasta'):
					if rec.id.split('|')[0] != hg: continue
					_, allele_cluster, _, g = rec.id.split('|')
					ginfo = comp_gene_info[g]
					gene_sequence[g] = str(rec.seq)
					gstart = ginfo['start']
					gend = ginfo['end']

					gene_length = gend - gstart + 1

					gene_covered_1 = 0

					try:
						for pileupcolumn in bam_handle.pileup(contig=rec.id, stepper="nofilter"):
							pos_depth = 0
							for pileupread in pileupcolumn.pileups:
								if pileupread.is_del or pileupread.is_refskip: continue
								read = pileupread.alignment
								if read.query_qualities[pileupread.query_position] < 30: continue
								pos_depth += 1
							if pos_depth >= 1:
								gene_covered_1 += 1
					except:
						pass

					gene_coverage_1 = gene_covered_1 / float(gene_length)
					if gene_coverage_1 < 0.90: continue
					hg_genes_covered += 1
					#print('\t'.join([sample, hg, rec.id, str(gene_coverage_1), str(gene_coverage_3)]))

					for read_alignment in bam_handle.fetch(rec.id):
						read_name = read_alignment.query_name
						total_reads.add(read_name)
						read_ascore = read_alignment.tags[0][1]

						read_ref_positions = set(read_alignment.get_reference_positions())

						first_real_alignment_pos = None
						last_real_alignment_pos = None
						indel_positions = set([])
						matches = set([])
						for b in read_alignment.get_aligned_pairs(with_seq=True):
							if not (b[0] == None or b[1] == None):
								if first_real_alignment_pos == None:
									first_real_alignment_pos = b[1]
								last_real_alignment_pos = b[1]
								if b[2].isupper():
									matches.add(b[1])
							else:
								indel_positions.add(b[1])

						main_alignment_positions = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
						read_has_indel = len(main_alignment_positions.intersection(indel_positions)) > 0
						matching_percentage = float(len(matches))/float(len(main_alignment_positions))

						read_ascores_per_allele[read_name].append([g, read_ascore, matching_percentage, len(main_alignment_positions), read_has_indel, read_alignment])

			accounted_reads = set([])
			hg_align_pos_alleles = defaultdict(lambda: defaultdict(set))
			supported_snvs = defaultdict(lambda: defaultdict(set))
			for read in read_ascores_per_allele:
				top_score = -1000000
				score_sorted_alignments = sorted(read_ascores_per_allele[read], key=itemgetter(1), reverse=True)
				for i, align in enumerate(score_sorted_alignments):
					if i == 0: top_score = align[1]
					if align[1] == top_score and ((align[2] >= 0.99 and align[3] >= 60) or (align[2] >= 0.9 and align[3] >= 100)):# and align[4] == False:
						read_alignment = align[-1]
						topaligns_handle.write(read_alignment)

						min_read_ref_pos = min(read_alignment.get_reference_positions())
						read_referseq = read_alignment.get_reference_sequence().upper()
						read_queryseq = read_alignment.query_sequence
						read_queryqua = read_alignment.query_qualities

						for b in read_alignment.get_aligned_pairs(with_seq=True):
							if b[0] == None or b[1] == None: continue
							ref_pos = b[1]+1
							alt_al = read_queryseq[b[0]].upper()
							ref_al = read_referseq[b[1] - min_read_ref_pos].upper()
							assert (ref_al == str(gene_sequence[align[0]]).upper()[b[1]])
							if b[2] == 'n' or ref_al == 'N' or alt_al == 'N': continue
							que_qual = read_queryqua[b[0]]
							if (que_qual >= 30) and ((ref_pos+3) < len(gene_sequence[align[0]])):
								cod_pos = gene_pos_to_msa_pos[hg][align[0]][ref_pos]
								hg_align_pos_alleles[cod_pos][alt_al].add(read)
								accounted_reads.add(read)
								det_outf.write('\t'.join([str(x) for x in [sample, hg, align[0], ref_pos, cod_pos, ref_al, alt_al, read, align[1], align[2], align[3], align[4]]]) + '\n')
								if b[2].islower():
									assert (alt_al != ref_al)
									snv_id = str(read_alignment.reference_name) + '_|_' + str(ref_pos) + '_|_' + ref_al + '_|_' + alt_al
									supported_snvs[snv_id][read].add(align[1])

			for pos in range(1, codon_alignment_lengths[hg]+1):
				printlist = [hg, str(pos)]
				for al in ['A', 'C', 'G', 'T']:
					printlist.append(str(len(hg_align_pos_alleles[pos][al])))
				res_outf.write(','.join(printlist) + '\n')

			for snv in supported_snvs:
				support_info = []
				for read in supported_snvs[snv]:
					support_info.append(read + '_|_' + str(max(supported_snvs[snv][read])))
				snv_outf.write('\t'.join([snv, str(len(support_info))] + support_info) + '\n')

			print('\t'.join([hg, str(len(total_reads)), str(len(accounted_reads))]))

		snv_outf.close()
		res_outf.close()
		det_outf.close()
		topaligns_handle.close()
		bam_handle.close()

		os.system("samtools sort -@ %d %s -o %s" % (1, topaligns_file, topaligns_file_sorted))
		os.system("samtools index %s" % topaligns_file_sorted)

	except Exception as e:
		if logObject:
			logObject.error("Issues with mining for SNVs/parsing alignment file.")
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())


def paired_snv_miner(input_args):
	"""
	Function to mine for novel SNVs and identify alleles of homolog groups for GCF present in paired-end sequencing
	dataset.
	"""
	sample, bam_alignment, ref_fasta, hg_gene_to_rep, res_dir, bgc_hg_genes, comp_gene_info, gene_pos_to_msa_pos, gene_pos_to_allele, codon_alignment_lengths, logObject = input_args
	try:
		hg_rep_genes = defaultdict(set)
		for g, r in hg_gene_to_rep.items():
			hg_rep_genes[r].add(g)

		if not os.path.isfile(bam_alignment): return
		snvs_file = res_dir + sample + '.snvs'
		result_file = res_dir + sample + '.txt'
		details_file = res_dir + sample + '.full.txt'
		snv_outf = open(snvs_file, 'w')
		res_outf = open(result_file, 'w')
		det_outf = open(details_file, 'w')

		res_outf.write('Contig,Position,Sample-A,Sample-C,Sample-G,Sample-T\n')

		bam_handle = pysam.AlignmentFile(bam_alignment, 'rb')

		topaligns_file = res_dir + sample + '_topaligns.bam'
		topaligns_file_sorted = res_dir + sample + '_topaligns.sorted.bam'
		topaligns_handle = pysam.AlignmentFile(topaligns_file, "wb", template=bam_handle)

		for hg, hg_genes in bgc_hg_genes.items():
			read_ascores_per_allele = defaultdict(list)
			hg_genes_covered = 0
			gene_sequence = {}
			total_reads = set([])
			with open(ref_fasta) as opff:
				for rec in SeqIO.parse(opff, 'fasta'):
					if rec.id.split('|')[0] != hg: continue
					_, allele_cluster, _, g = rec.id.split('|')
					ginfo = comp_gene_info[g]
					gene_sequence[g] = str(rec.seq)
					gstart = ginfo['start']
					gend = ginfo['end']

					gene_length = gend - gstart + 1

					gene_covered_1 = 0

					try:
						for pileupcolumn in bam_handle.pileup(contig=rec.id, stepper="nofilter"):
							pos_depth = 0
							for pileupread in pileupcolumn.pileups:
								if pileupread.is_del or pileupread.is_refskip: continue
								read = pileupread.alignment
								if read.query_qualities[pileupread.query_position] < 30: continue
								pos_depth += 1
							if pos_depth >= 1:
								gene_covered_1 += 1
					except:
						pass

					gene_coverage_1 = gene_covered_1 / float(gene_length)
					if gene_coverage_1 < 0.90: continue
					hg_genes_covered += 1
					#print('\t'.join([sample, hg, rec.id, str(gene_coverage_1), str(gene_coverage_3)]))

					for read1_alignment, read2_alignment in util.read_pair_generator(bam_handle, rec.id, gene_length):
						if read1_alignment and read2_alignment:
							read_name = read1_alignment.query_name
							total_reads.add(read_name)
							combined_ascore = read1_alignment.tags[0][1] + read2_alignment.tags[0][1]
							read1_ref_positions = set(read1_alignment.get_reference_positions())
							read2_ref_positions = set(read2_alignment.get_reference_positions())

							align_union = len(read1_ref_positions.union(read2_ref_positions))
							align_intersect = len(read1_ref_positions.intersection(read2_ref_positions))
							min_align_length = min(len(read1_ref_positions), len(read2_ref_positions))
							align_overlap_prop = float(align_intersect) / float(align_union) # / min_align_length

							matches_1 = set([])
							first_real_alignment_pos = None
							last_real_alignment_pos = None
							indel_positions = set([])
							for b in read1_alignment.get_aligned_pairs(with_seq=True):
								if not (b[0] == None or b[1] == None):
									if first_real_alignment_pos == None:
										first_real_alignment_pos = b[1]
									last_real_alignment_pos = b[1]
									if b[2].isupper():
										matches_1.add(b[1])
								else:
									indel_positions.add(b[1])

							main_alignment_positions_1 = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
							read1_has_indel = len(main_alignment_positions_1.intersection(indel_positions)) > 0

							matches_2 = set([])
							first_real_alignment_pos = None
							last_real_alignment_pos = None
							indel_positions = set([])
							for b in read2_alignment.get_aligned_pairs(with_seq=True):
								if not (b[0] == None or b[1] == None):
									if first_real_alignment_pos == None:
										first_real_alignment_pos = b[1]
									last_real_alignment_pos = b[1]
									if b[2].isupper():
										matches_2.add(b[1])
								else:
									indel_positions.add(b[1])


							main_alignment_positions_2 = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
							read2_has_indel = len(main_alignment_positions_2.intersection(indel_positions)) > 0

							matches = matches_1.union(matches_2)
							main_alignment_positions = main_alignment_positions_1.union(main_alignment_positions_2)

							matching_percentage = float(len(matches))/float(len(main_alignment_positions))
							matching_percentage_read1 = float(len(matches_1))/float(len(main_alignment_positions_1))
							matching_percentage_read2 = float(len(matches_2))/float(len(main_alignment_positions_2))

							#if matching_percentage < 0.95: continue
							#if read1_has_indel or read2_has_indel or matching_percentage < 0.95 or align_overlap_prop > 0.75: continue

							read_ascores_per_allele[read_name].append(
								[g, combined_ascore, matching_percentage, # max([matching_percentage, matching_percentage_read1, matching_percentage_read2]) ,
								 len(main_alignment_positions),
								 (read1_has_indel or read2_has_indel),
								 abs(matching_percentage_read1 - matching_percentage_read2), abs(len(main_alignment_positions_1) - len(main_alignment_positions_2))/float(len(main_alignment_positions)),
								 [read1_alignment, read2_alignment]])
							"""
							if max([matching_percentage, matching_percentage_read1, matching_percentage_read2]) == matching_percentage:
								read_ascores_per_allele[read_name].append([g, combined_ascore, matching_percentage, len(main_alignment_positions), (read1_has_indel or read2_has_indel), abs(matching_percentage_read1-matching_percentage_read2), [read1_alignment, read2_alignment]])
							elif max([matching_percentage, matching_percentage_read1, matching_percentage_read2]) == matching_percentage_read1:
								read_ascores_per_allele[read_name].append([g, combined_ascore, matching_percentage_read1, len(main_alignment_positions), (read1_has_indel or read2_has_indel), 0.0, [read1_alignment]])
							elif max([matching_percentage, matching_percentage_read1, matching_percentage_read2]) == matching_percentage_read2:
								read_ascores_per_allele[read_name].append([g, combined_ascore, matching_percentage_read2, len(main_alignment_positions), (read1_has_indel or read2_has_indel), 0.0, [read2_alignment]])
							"""
						elif read1_alignment or read2_alignment:
							read_alignment = None
							if read1_alignment != None:
								read_alignment = read1_alignment
							else:
								read_alignment = read2_alignment

							read_name = read_alignment.query_name
							total_reads.add(read_name)

							matches = set([])
							first_real_alignment_pos = None
							last_real_alignment_pos = None
							indel_positions = set([])
							for b in read.get_aligned_pairs(with_seq=True):
								if b[0] != None and b[1] != None and b[2] != None:
									if first_real_alignment_pos == None:
										first_real_alignment_pos = b[1]
									last_real_alignment_pos = b[1]
									if b[2].isupper():
										matches.add(b[1])
								else:
									indel_positions.add(b[1])
							main_alignment_positions = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
							has_indel = len(main_alignment_positions.intersection(indel_positions)) > 0

							read_length = len(set(read_alignment.get_reference_positions()))
							matching_percentage = float(len(matches))/float(len(main_alignment_positions))

							# 0 added just to signify that it is a single mate contributing to the paired end combined ascore
							combined_ascore = 0 + read_alignment.tags[0][1]

							read_ascores_per_allele[read_name].append([g, combined_ascore, matching_percentage, len(main_alignment_positions), has_indel, 0.0, 0.0, [read_alignment]])

			accounted_reads = set([])
			hg_align_pos_alleles = defaultdict(lambda: defaultdict(set))
			supported_snvs = defaultdict(lambda: defaultdict(set))
			for read in read_ascores_per_allele:
				top_score = -1000000
				score_sorted_alignments = sorted(read_ascores_per_allele[read], key=itemgetter(1), reverse=True)
				for i, align in enumerate(score_sorted_alignments):
					if i == 0: top_score = align[1]
					if align[1] == top_score and ((align[2] >= 0.99 and align[3] >= 100) or (align[2] >= 0.95 and align[3] >= 180)) and align[4] == False and align[5] < 0.02: # and align[6] < 0.25:
						for read_alignment in align[-1]:
							topaligns_handle.write(read_alignment)

							min_read_ref_pos = min(read_alignment.get_reference_positions())
							read_referseq = read_alignment.get_reference_sequence().upper()
							read_queryseq = read_alignment.query_sequence
							read_queryqua = read_alignment.query_qualities

							for b in read_alignment.get_aligned_pairs(with_seq=True):
								if b[0] == None or b[1] == None: continue
								ref_pos = b[1]+1
								alt_al = read_queryseq[b[0]].upper()
								ref_al = read_referseq[b[1] - min_read_ref_pos].upper()
								assert (ref_al == str(gene_sequence[align[0]]).upper()[b[1]])
								if b[2] == 'n' or ref_al == 'N' or alt_al == 'N': continue
								que_qual = read_queryqua[b[0]]
								if (que_qual >= 30) and ((ref_pos+3) < len(gene_sequence[align[0]])):
									cod_pos = gene_pos_to_msa_pos[hg][align[0]][ref_pos]
									hg_align_pos_alleles[cod_pos][alt_al].add(read)
									accounted_reads.add(read)
									det_outf.write('\t'.join([str(x) for x in [sample, hg, align[0], ref_pos, cod_pos, ref_al, alt_al, read, align[1], align[2], align[3], align[4]]]) + '\n')
									if b[2].islower():
										assert (alt_al != ref_al)
										snv_id = str(read_alignment.reference_name) + '_|_' + str(ref_pos) + '_|_' + ref_al + '_|_' + alt_al
										supported_snvs[snv_id][read].add(align[1])

			for pos in range(1, codon_alignment_lengths[hg]+1):
				printlist = [hg, str(pos)]
				for al in ['A', 'C', 'G', 'T']:
					printlist.append(str(len(hg_align_pos_alleles[pos][al])))
				res_outf.write(','.join(printlist) + '\n')

			for snv in supported_snvs:
				support_info = []
				for read in supported_snvs[snv]:
					support_info.append(read + '_|_' + str(max(supported_snvs[snv][read])))
				snv_outf.write('\t'.join([snv, str(len(support_info))] + support_info) + '\n')

			print('\t'.join([hg, str(len(total_reads)), str(len(accounted_reads))]))

		snv_outf.close()
		res_outf.close()
		det_outf.close()
		topaligns_handle.close()
		bam_handle.close()

		os.system("samtools sort -@ %d %s -o %s" % (1, topaligns_file, topaligns_file_sorted))
		os.system("samtools index %s" % topaligns_file_sorted)

	except Exception as e:
		if logObject:
			logObject.error("Issues with mining for SNVs/parsing alignment file.")
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def popgen_analysis_of_hg(inputs):
	"""
	Helper function which is to be called from the runPopulationGeneticsAnalysis() function to parallelize population
	genetics analysis of each homolog group.

	:param inputs: list of inputs passed in by GCF.runPopulationGeneticsAnalysis().
	"""
	gcf_id, hg, codon_alignment_fasta, popgen_dir, plots_dir, comp_gene_info, hg_genes, bgc_sample, hg_prop_multi_copy, hg_order_scores, sample_population, population, logObject = inputs
	domain_plot_file = plots_dir + hg + '_domain.txt'
	position_plot_file = plots_dir + hg + '_position.txt'
	popgen_plot_file = plots_dir + hg + '_popgen.txt'
	plot_pdf_file = plots_dir + hg + '.pdf'

	domain_plot_handle = open(domain_plot_file, 'w')
	position_plot_handle = open(position_plot_file, 'w')
	popgen_plot_handle = open(popgen_plot_file, 'w')

	seqs = []
	samples_ordered = []
	genes_ordered = []
	bgc_codons = defaultdict(list)
	num_codons = None
	samples = set([])
	gene_lengths = []
	gene_locs = defaultdict(dict)
	core_counts = defaultdict(int)
	products = set([])

	updated_codon_alignment_fasta = popgen_dir + codon_alignment_fasta.split('/')[-1]
	updated_codon_alignment_handle = open(updated_codon_alignment_fasta, 'w')

	with open(codon_alignment_fasta) as ocaf:
		for rec in SeqIO.parse(ocaf, 'fasta'):
			sample_id, gene_id = rec.id.split('|')
			if len(gene_id.split('_')[0]) == 3:
				if comp_gene_info[gene_id]['core_overlap']:
					core_counts['core'] += 1
				else:
					core_counts['auxiliary'] += 1
			if population != sample_population[sample_id] and population != None: continue
			updated_codon_alignment_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
			products.add(comp_gene_info[gene_id]['product'])
			real_pos = 1
			seqs.append(list(str(rec.seq)))
			codons = [str(rec.seq)[i:i + 3] for i in range(0, len(str(rec.seq)), 3)]
			num_codons = len(codons)
			bgc_codons[rec.id] = codons
			samples.add(sample_id)
			samples_ordered.append(sample_id)
			genes_ordered.append(gene_id)
			for msa_pos, bp in enumerate(str(rec.seq)):
				if bp != '-':
					gene_locs[gene_id][real_pos] = msa_pos + 1
					real_pos += 1
			gene_lengths.append(len(str(rec.seq).replace('-', '')))

	updated_codon_alignment_handle.close()
	codon_alignment_fasta = updated_codon_alignment_fasta
	if len(seqs) == 0: return

	is_core = False
	if float(core_counts['core']) / sum(core_counts.values()) >= 0.8: is_core = True

	median_gene_length = statistics.median(gene_lengths)

	variable_sites = set([])
	conserved_sites = set([])
	position_plot_handle.write('\t'.join(['pos', 'num_seqs', 'num_alleles', 'num_gaps', 'maj_allele_freq']) + '\n')

	sample_differences_to_consensus = defaultdict(lambda: defaultdict(int))
	for i, ls in enumerate(zip(*seqs)):
		al_counts = defaultdict(int)
		for al in ls:
			if al != '-': al_counts[al] += 1
		if sum(al_counts.values()) == 0: continue
		maj_allele_count = max(al_counts.values())
		tot_count = sum(al_counts.values())
		num_seqs = len(ls)
		num_alleles = len(al_counts.keys())
		num_gaps = num_seqs - tot_count
		maj_allele_freq = float(maj_allele_count) / tot_count
		maj_allele = [a[0] for a in al_counts.items() if a[0] != '-' and maj_allele_count == a[1]][0]
		for j, al in enumerate(ls):
			sid = samples_ordered[j]
			gid = genes_ordered[j]
			if al != maj_allele:
				sample_differences_to_consensus[sid][gid] += 1
			else:
				sample_differences_to_consensus[sid][gid] += 0
		position_plot_handle.write('\t'.join([str(x) for x in [i + 1, num_seqs, num_alleles, num_gaps, maj_allele_freq]]) + '\n')
		if maj_allele_freq <= 0.90:
			variable_sites.add(i)
		else:
			conserved_sites.add(i)
	position_plot_handle.close()

	differential_domains = set([])
	domain_positions_msa = defaultdict(set)
	domain_min_position_msa = defaultdict(lambda: 1e8)
	all_domains = set([])
	for gene in gene_locs:
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
					domain_plot_handle.write(
						'\t'.join([str(x) for x in [dom[0], i, min_pos, max_pos]]) + '\n')
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
		cods = set([])
		for bgc in bgc_codons:
			cod = bgc_codons[bgc][cod_index]
			cods.add(cod)
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

		if len(cod_count) == 0: continue
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

	rscript_plot_cmd = ["Rscript", RSCRIPT_FOR_CLUSTER_ASSESSMENT_PLOTTING, domain_plot_file, position_plot_file,
						popgen_plot_file,
						plot_pdf_file]
	if logObject:
		logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
	try:
		subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
	except Exception as e:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	tajima_results = popgen_dir + hg + '.tajima.txt'
	if len(seqs) >= 3:
		rscript_tajimaD_cmd = ["Rscript", RSCRIPT_FOR_TAJIMA, codon_alignment_fasta, tajima_results]
		if logObject:
			logObject.info('Running R pegas for calculating Tajima\'s D from codon alignment with the following command: %s' % ' '.join(
				rscript_tajimaD_cmd))
		try:
			subprocess.call(' '.join(rscript_tajimaD_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if logObject:
				logObject.info('Successfully ran: %s' % ' '.join(rscript_tajimaD_cmd))
		except Exception as e:
			if logObject:
				logObject.error('Had an issue running: %s' % ' '.join(rscript_tajimaD_cmd))
				logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	tajimas_d = "NA"
	if os.path.isfile(tajima_results):
		with open(tajima_results) as otrf:
			for i, line in enumerate(otrf):
				if i == 4:
					try:
						tajimas_d = float(line.strip())
					except:
						pass

	prop_samples_with_hg = len(samples) / float(len(set(bgc_sample.values())))

	hg_info = []
	if population:
		hg_info = [population]
	hg_info += [gcf_id, hg, '; '.join(products), hg_order_scores[hg], hg_prop_multi_copy[hg],
				median_gene_length, is_core, len(seqs), len(samples), prop_samples_with_hg, tajimas_d, total_core_codons,
				total_variable_codons, nonsynonymous_sites, synonymous_sites, dn_ds, '; '.join(all_domains)]

	if sample_population and not population:
		input_anova_data_seqsim = []
		input_anova_header_seqsim = ['sample', 'population', 'differences_to_consensus']

		population_counts = defaultdict(int)
		for s, p in sample_population.items():
			population_counts[p] += 1

		pops_with_hg = set([])
		pop_count_with_hg = defaultdict(int)
		for s in sample_differences_to_consensus:
			min_diff_to_consensus = 1e100
			for g in sample_differences_to_consensus[s]:
				if min_diff_to_consensus > sample_differences_to_consensus[s][g]:
					min_diff_to_consensus = sample_differences_to_consensus[s][g]
			if min_diff_to_consensus < 1e100:
				pop_count_with_hg[sample_population[s]] += 1
				data_row = [s, sample_population[s], float(min_diff_to_consensus)]
				input_anova_data_seqsim.append(data_row)
				pops_with_hg.add(sample_population[s])

		anova_pval_seqsim = "NA"
		if len(pops_with_hg) >= 2:
			anova_input_df = DataFrame(np.array(input_anova_data_seqsim), columns=input_anova_header_seqsim)
			anova_input_df['differences_to_consensus'] = anova_input_df['differences_to_consensus'].astype(float)
			aov = welch_anova(dv='differences_to_consensus', between='population', data=anova_input_df)
			anova_pval_seqsim = aov.iloc[0, 4]

		fishers_pvals = []
		for pop in pop_count_with_hg:
			other_count = sum([x[1] for x in pop_count_with_hg.items() if x[0] != pop])
			other_total = sum([x[1] for x in population_counts.items() if x[0] != pop])
			odds, pval = fisher_exact([[pop_count_with_hg[pop], population_counts[pop]], [other_count, other_total]])
			fishers_pvals.append(pval)

		fisher_pval = "NA"
		if len(fishers_pvals) > 0: fisher_pval = min(fishers_pvals)

		hg_population_info = [len(pops_with_hg), '|'.join([str(x[0]) + '=' + str(float(x[1])/population_counts[x[0]]) for x in pop_count_with_hg.items()]),
								 fisher_pval, str(anova_pval_seqsim).replace('nan', 'NA')]
		hg_info += hg_population_info

	hg_stats_handle = open(popgen_dir + hg + '_stats.txt', 'w')
	hg_stats_handle.write('\t'.join([str(x) for x in hg_info]) + '\n')
	hg_stats_handle.close()

def create_codon_msas(inputs):
	"""
	Helper function which is to be called from the constructCodonAlignments() function to parallelize construction
	of codon alignments for each homolog group of interest in the GCF.
	:param inputs: list of inputs passed in by GCF.constructCodonAlignments().
	"""
	hg, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, logObject = inputs

	hg_nucl_fasta = nucl_seq_dir + '/' + hg + '.fna'
	hg_prot_fasta = prot_seq_dir + '/' + hg + '.faa'
	hg_prot_msa = prot_alg_dir + '/' + hg + '.msa.faa'
	hg_codo_msa = codo_alg_dir + '/' + hg + '.msa.fna'

	hg_nucl_handle = open(hg_nucl_fasta, 'w')
	hg_prot_handle = open(hg_prot_fasta, 'w')
	for s in gene_sequences:
		hg_nucl_handle.write('>' + s + '\n' + str(gene_sequences[s][0]) + '\n')
		hg_prot_handle.write('>' + s + '\n' + str(gene_sequences[s][1]) + '\n')
	hg_nucl_handle.close()
	hg_prot_handle.close()

	mafft_cmd = ['mafft', '--maxiterate', '1000', '--localpair', hg_prot_fasta, '>', hg_prot_msa]
	pal2nal_cmd = ['pal2nal.pl', hg_prot_msa, hg_nucl_fasta, '-output', 'fasta', '>', hg_codo_msa]

	if logObject:
		logObject.info('Running mafft with the following command: %s' % ' '.join(mafft_cmd))
	try:
		subprocess.call(' '.join(mafft_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(mafft_cmd))
	except Exception as e:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(mafft_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(mafft_cmd))

	if logObject:
		logObject.info('Running PAL2NAL with the following command: %s' % ' '.join(pal2nal_cmd))
	try:
		subprocess.call(' '.join(pal2nal_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(pal2nal_cmd))
	except Exception as e:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(pal2nal_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(pal2nal_cmd))

	if logObject:
		logObject.info('Achieved codon alignment for homolog group %s' % hg)
