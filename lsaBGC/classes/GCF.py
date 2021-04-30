import os
import sys
import logging
import traceback
import statistics
import random
import subprocess
import pysam
import multiprocessing
from scipy.stats import f_oneway
from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter
from collections import defaultdict
from lsaBGC.classes.Pan import Pan
from lsaBGC import util

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

		# Sequence and alignment directories
		self.nucl_seq_dir = None
		self.prot_seq_dir = None
		self.prot_alg_dir = None
		self.codo_alg_dir = None

		# Concatenated HMMER3 HMM profiles database of homolog groups in GCF
		self.concatenated_profile_HMM = None

		# Dictionary of individual genes to haplotype/allelic representative gene
		self.instance_to_haplotype = {}

		# Set of samples with sequencing reads to avoid for reporting alleles and novel SNVs,
		# these samples do not exhibit enough support for harboring a full BGC for the GC
		self.avoid_samples = set([])

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
		self.hg_to_color = hg_to_color

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
			for leaf in tree_obj:
				bgc_weights[str(leaf).strip('\n').lstrip('-')] += 1

			bgc_hg_presence = defaultdict(lambda: defaultdict(lambda: 'Absent'))
			hg_counts = defaultdict(int)
			for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
				bgc = item[0]
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
			gggenes_track_handle.close()

			for bgc in bgc_hg_presence:
				for hg in hg_counts:
					heatmap_track_handle.write(
						'\t'.join([bgc, hg, bgc_hg_presence[bgc][hg], str(hg_counts[hg])]) + '\n')
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

	def constructCodonAlignments(self, outdir, cores=1, only_scc=False, list_alignments=False):
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

		nucl_seq_dir = os.path.abspath(outdir + 'Nucleotide_Sequences') + '/'
		prot_seq_dir = os.path.abspath(outdir + 'Protein_Sequences') + '/'
		prot_alg_dir = os.path.abspath(outdir + 'Protein_Alignments') + '/'
		codo_alg_dir = os.path.abspath(outdir + 'Codon_Alignments') + '/'

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
				inputs.append(
					[hg, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, self.logObject])

			p = multiprocessing.Pool(cores)
			p.map(create_codon_msas, inputs)

			self.nucl_seq_dir = nucl_seq_dir
			self.prot_seq_dir = prot_seq_dir
			self.prot_alg_dir = prot_alg_dir
			self.codo_alg_dir = codo_alg_dir

			if list_alignments:
				codon_alg_listings_file = outdir + 'Codon_Alignments_Listings.txt'
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
					print(f)
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

	def runPopulationGeneticsAnalysis(self, outdir, cores=1):
		"""
		Wrapper function which serves to parallelize population genetics analysis.

		:param outdir: The path to the workspace / output directory.
		:param cores: The number of cores (will be used for parallelizing)
		"""

		popgen_dir = outdir + 'Codon_PopGen_Analyses/'
		plots_dir = outdir + 'Codon_MSA_Plots/'
		if not os.path.isdir(popgen_dir): os.system('mkdir %s' % popgen_dir)
		if not os.path.isdir(plots_dir): os.system('mkdir %s' % plots_dir)

		final_output_handle = open(outdir + 'Ortholog_Group_Information.txt', 'w')
		header = ['gcf_id', 'homolog_group', 'annotation', 'hg_order_index', 'hg_median_copy_count', 'median_gene_length',
				  'is_core_to_bgc', 'bgcs_with_hg', 'proportion_of_samples_with_hg', 'Tajimas_D', 'core_codons',
				  'total_variable_codons', 'nonsynonymous_codons', 'synonymous_codons', 'dn_ds', 'all_domains']
		if self.bgc_population != None:
			header += ['populations_with_hg', 'population_proportion_of_members_with_hg', 'one_way_ANOVA_pvalues']

		final_output_handle.write('\t'.join(header) + '\n')

		inputs = []
		for f in os.listdir(self.codo_alg_dir):
			hg = f.split('.msa.fna')[0]
			codon_alignment_fasta = self.codo_alg_dir + f
			inputs.append([self.gcf_id, hg, codon_alignment_fasta, popgen_dir, plots_dir, self.comp_gene_info, self.hg_genes,
						   self.bgc_sample, self.hg_prop_multi_copy, self.hg_order_scores, dict(self.sample_population),
						   self.logObject])

		p = multiprocessing.Pool(cores)
		p.map(popgen_analysis_of_hg, inputs)

		for f in os.listdir(popgen_dir):
			if not f.endswith('_stats.txt'): continue
			with open(popgen_dir + f) as opf:
				for line in opf:
					line = line
					final_output_handle.write(line)
		final_output_handle.close()

	def constructHMMProfiles(self, outdir, cores=1):
		"""
		Wrapper function to construct Hmmer3 HMMs for each of the homolog groups.

		:param outdir: The path to the workspace / output directory.
		:param cores: The number of cores (will be used for parallelizing)
		"""

		prot_seq_dir = os.path.abspath(outdir + 'Protein_Sequences') + '/'
		prot_alg_dir = os.path.abspath(outdir + 'Protein_Alignments') + '/'
		prot_hmm_dir = os.path.abspath(outdir + 'Profile_HMMs') + '/'
		if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
		if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
		if not os.path.isdir(prot_hmm_dir): os.system('mkdir %s' % prot_hmm_dir)

		all_samples = set(self.bgc_sample.values())
		try:
			inputs = []
			for hg in self.hg_genes:
				sample_counts = defaultdict(int)
				sample_sequences = {}
				for gene in self.hg_genes[hg]:
					gene_info = self.comp_gene_info[gene]
					bgc_id = gene_info['bgc_name']
					sample_id = self.bgc_sample[bgc_id]
					prot_seq = gene_info['prot_seq']
					sample_counts[sample_id] += 1
					sample_sequences[sample_id] = prot_seq
				samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
				samples_with_any_copy = set([s[0] for s in sample_counts.items() if s[1] > 0])

				# check that hg is single-copy-core or just core
				if len(samples_with_single_copy.symmetric_difference(all_samples)) == 0:
					self.scc_homologs.add(hg)
				if len(samples_with_any_copy.symmetric_difference(all_samples)) == 0:
					self.core_homologs.add(hg)
					if self.hg_prop_multi_copy[hg] <= 0.05:
						self.specific_core_homologs.add(hg)
				inputs.append([hg, sample_sequences, prot_seq_dir, prot_alg_dir, prot_hmm_dir, self.logObject])

			p = multiprocessing.Pool(cores)
			p.map(create_hmm_profiles, inputs)
			p.close()

			if self.logObject:
				self.logObject.info(
					"Successfully created profile HMMs for each homolog group. Now beginning concatenation into single file.")
			self.concatenated_profile_HMM = outdir + 'All_GCF_Homologs.hmm'
			os.system('rm -f %s %s.h3*' % (self.concatenated_profile_HMM, self.concatenated_profile_HMM))
			for f in os.listdir(prot_hmm_dir):
				os.system('cat %s >> %s' % (prot_hmm_dir + f, self.concatenated_profile_HMM))

			hmmpress_cmd = ['hmmpress', self.concatenated_profile_HMM]
			if self.logObject:
				self.logObject.info(
				'Running hmmpress on concatenated profiles with the following command: %s' % ' '.join(hmmpress_cmd))
			try:
				subprocess.call(' '.join(hmmpress_cmd), shell=True, stdout=subprocess.DEVNULL,
								stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				if self.logObject:
					self.logObject.info('Successfully ran: %s' % ' '.join(hmmpress_cmd))
			except:
				if self.logObject:
					self.logObject.error('Had an issue running: %s' % ' '.join(hmmpress_cmd))
					self.logObject.error(traceback.format_exc())
				raise RuntimeError('Had an issue running: %s' % ' '.join(hmmpress_cmd))

		except:
			if self.logObject:
				self.logObject.error("Issues with running hmmpress on profile HMMs.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def runHMMScanAndAssignBGCsToGCF(self, outdir, sample_prokka_data, orthofinder_matrix_file, cores=1):
		"""

		"""
		search_res_dir = os.path.abspath(outdir + 'HMMScan_Results') + '/'
		bgc_genbanks_dir = os.path.abspath(outdir + 'BGC_Genbanks') + '/'
		if not os.path.isdir(search_res_dir): os.system('mkdir %s' % search_res_dir)
		if not os.path.isdir(bgc_genbanks_dir): os.system('mkdir %s' % bgc_genbanks_dir)

		sample_bgc_ids = defaultdict(lambda: 1)
		boundary_genes = defaultdict(set)
		scaffold_genes = defaultdict(lambda: defaultdict(set))
		gene_location = defaultdict(dict)
		gene_id_to_order = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
		gene_order_to_id = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

		tot_bgc_proteins = defaultdict(int)
		hmmscan_cmds = []
		for sample in sample_prokka_data:
			sample_genbank = sample_prokka_data[sample]['genbank']
			sample_proteome = sample_prokka_data[sample]['predicted_proteome']

			gene_to_scaff, scaff_genes, bound_genes, gito, goti = util.parseGenbankAndFindBoundaryGenes(sample_genbank)

			boundary_genes[sample] = bound_genes
			scaffold_genes[sample] = scaff_genes
			gene_location[sample] = gene_to_scaff
			gene_id_to_order[sample] = gito
			gene_order_to_id[sample] = goti

			result_file = search_res_dir + sample + '.txt'
			hmmscan_cmd = ['hmmscan', '--max', '--cpu', '1', '--tblout', result_file, self.concatenated_profile_HMM, sample_proteome, self.logObject]
			hmmscan_cmds.append(hmmscan_cmd)

		p = multiprocessing.Pool(cores)
		p.map(util.multiProcess, hmmscan_cmds)
		p.close()

		protein_hits = defaultdict(list)
		sample_hg_hits_with_evalue = defaultdict(lambda: defaultdict(list))
		scaffold_proteins = defaultdict(lambda: defaultdict(set))
		sample_protein_to_hg = defaultdict(dict)
		sample_hgs = defaultdict(set)
		for sample in sample_prokka_data:
			result_file = search_res_dir + sample + '.txt'
			assert (os.path.isfile(result_file))

			with open(result_file) as orf:
				for line in orf:
					if line.startswith("#"): continue
					line = line.strip()
					ls = line.split()
					hg = ls[0]
					gene_id = ls[2]
					scaffold = gene_location[sample][gene_id]['scaffold']
					eval = float(ls[4])
					if eval <= 1e-5:
						protein_hits[gene_id].append([hg, eval, sample, scaffold])

		for p in protein_hits:
			for i, hits in enumerate(sorted(protein_hits[p], key=itemgetter(1))):
				if i == 0:
					sample_protein_to_hg[hits[2]][p] = hits[0]
					sample_hgs[hits[2]].add(hits[0])
					scaffold_proteins[hits[2]][hits[3]].add(p)
					sample_hg_hits_with_evalue[hits[2]][hits[0]].append([p, hits[1], hits[3]])

		expanded_gcf_list_file = outdir + 'GCF_Expanded.txt'
		expanded_gcf_list_handle = open(expanded_gcf_list_file, 'w')

		with open(self.bgc_genbanks_listing) as obglf:
			for line in obglf:
				expanded_gcf_list_handle.write(line)

		all_samples = set([])
		sample_hg_proteins = defaultdict(lambda: defaultdict(set))
		for sample in sample_hg_hits_with_evalue:
			sample_gcf_anchor_proteins = set([])
			scaffold_anchor_proteins = defaultdict(set)
			for hg in self.specific_core_homologs:
				for i, hits in enumerate(sorted(sample_hg_hits_with_evalue[sample][hg], key=itemgetter(1))):
					if i == 0:
						sample_gcf_anchor_proteins.add(hits[0])
						if hits[1] <= 1e-20:
							scaffold_anchor_proteins[hits[2]].add(hits[0])

			if len(sample_gcf_anchor_proteins) > 0:
				number_of_sample_scaffolds_with_anchor_hgs = len(set(scaffold_anchor_proteins.keys()))

				for scaffold in scaffold_anchor_proteins:
					max_anchor_pos = [0, None]
					min_anchor_pos = [100000000, None]

					for protein in list(scaffold_anchor_proteins[scaffold]):
						anchor_protein_start = gene_location[sample][protein]['start']
						anchor_protein_end = gene_location[sample][protein]['end']
						if anchor_protein_start < min_anchor_pos[0]: min_anchor_pos = [anchor_protein_start, protein]
						if anchor_protein_end > max_anchor_pos[0]: max_anchor_pos = [anchor_protein_end, protein]

					range_between_anchor_proteins = max_anchor_pos[0] - min_anchor_pos[0]
					if range_between_anchor_proteins <= 50000 and range_between_anchor_proteins > 0:
						start_anchor_gene_id = min_anchor_pos[1]
						end_anchor_gene_id = max_anchor_pos[1]
						start_anchor_gene_rank = gene_id_to_order[sample][scaffold][start_anchor_gene_id]
						end_anchor_gene_rank = gene_id_to_order[sample][scaffold][end_anchor_gene_id]

						scaffold_gcf_gene_ids = set([])
						for gene_rank in range(start_anchor_gene_rank, end_anchor_gene_rank + 1):
							scaffold_gcf_gene_ids.add(gene_order_to_id[sample][scaffold][gene_rank])

						# try expanding upstream
						strikes = 0
						curr_rank = start_anchor_gene_rank - 1
						while strikes < 3:
							if not curr_rank in gene_order_to_id[sample][scaffold]: break
							curr_rank_gene = gene_order_to_id[sample][scaffold][curr_rank]
							if curr_rank_gene in scaffold_proteins[sample][scaffold]:
								scaffold_gcf_gene_ids.add(curr_rank_gene)
								strikes = 0
							else:
								strikes += 1
							curr_rank -= 1

						# try expanding downstream
						strikes = 0
						curr_rank = end_anchor_gene_rank + 1
						while strikes < 3:
							if not curr_rank in gene_order_to_id[sample][scaffold]: break
							curr_rank_gene = gene_order_to_id[sample][scaffold][curr_rank]
							if curr_rank_gene in scaffold_proteins[sample][scaffold]:
								scaffold_gcf_gene_ids.add(curr_rank_gene)
								strikes = 0
							else:
								strikes += 1
							curr_rank += 1

						scaffold_hgs = set([])
						for g in scaffold_gcf_gene_ids:
							if g in sample_protein_to_hg[sample].keys():
								hg = sample_protein_to_hg[sample][g]
								scaffold_hgs.add(hg)

						if (number_of_sample_scaffolds_with_anchor_hgs == 1 or \
								(number_of_sample_scaffolds_with_anchor_hgs > 1 and \
								 len(boundary_genes[sample].intersection(scaffold_gcf_gene_ids)) > 0)) and \
								(len(sample_hgs[sample].intersection(self.core_homologs))/len(self.core_homologs) >= 0.9) and \
								(len(scaffold_hgs.intersection(self.core_homologs))/len(self.core_homologs) >= 0.2):
							clean_sample_name = sample.replace('-', '_').replace(':', '_').replace('.', '_').replace('=', '_')
							bgc_genbank_file = bgc_genbanks_dir + clean_sample_name + '_BGC-' + str(sample_bgc_ids[sample]) + '.gbk'
							sample_bgc_ids[sample] += 1

							min_bgc_pos = min([gene_location[sample][g]['start'] for g in scaffold_gcf_gene_ids])
							max_bgc_pos = max([gene_location[sample][g]['end'] for g in scaffold_gcf_gene_ids])
							util.createBGCGenbank(sample_prokka_data[sample]['genbank'], bgc_genbank_file, scaffold, min_bgc_pos, max_bgc_pos)
							expanded_gcf_list_handle.write('\t'.join([clean_sample_name, bgc_genbank_file]) + '\n')

							for g in scaffold_gcf_gene_ids:
								if g in sample_protein_to_hg[sample].keys():
									hg = sample_protein_to_hg[sample][g]
									sample_hg_proteins[clean_sample_name][hg].add(g)
							all_samples.add(clean_sample_name)
		expanded_gcf_list_handle.close()

		original_samples = []
		all_hgs = set([])
		with open(orthofinder_matrix_file) as omf:
			for i, line in enumerate(omf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					original_samples = ls[1:]
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
						#print('>' + gene + '|' + bgc + '|' + self.gene_to_hg[gene] + '\n' + self.comp_gene_info[gene]['nucl_seq_with_flanks'])
						gwff_handle.write('>' + gene + '|' + bgc + '|' + self.gene_to_hg[gene] + '\n' + self.comp_gene_info[gene]['nucl_seq_with_flanks'] + '\n')
			gwff_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Unable to extract flanking sequences of gene into FASTA file.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		cd_hit_nr = ['cd-hit-est', '-i', genes_with_flanks_fasta, '-o', cd_hit_nr_fasta_file, '-G', '1', '-g',
					 '1', '-d', '0', '-n', '10', '-M', '2000', '-c', '1.0', '-aL', '0.0', '-aS', '1.0', '-T', '1']
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

	def runSNVMining(self, paired_end_sequencing_file, bowtie2_ref_fasta, bowtie2_alignment_dir, results_dir, cores=1):
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
			process_args = []
			with open(paired_end_sequencing_file) as opesf:
				for line in opesf:
					line = line.strip()
					sample, frw_read, rev_read = line.split('\t')
					process_args.append([sample, bowtie2_alignment_dir + sample + '.filtered.sorted.bam',
										 bowtie2_ref_fasta, self.instance_to_haplotype, results_dir, self.hg_genes,
										 self.comp_gene_info, self.logObject])

			p = multiprocessing.Pool(cores)
			p.map(snv_miner, process_args)
			p.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def createSummaryMatricesForMetaNovelty(self, paired_end_sequencing_file, snv_mining_outdir, outdir):
		try:

			all_samples = set([])
			for hg in self.hg_genes:
				sample_counts = defaultdict(int)
				sample_sequences = {}
				for gene in self.hg_genes[hg]:
					if len(gene.split('_')[0]) == 3:
						gene_info = self.comp_gene_info[gene]
						bgc_id = gene_info['bgc_name']
						sample_id = self.bgc_sample[bgc_id]
						sample_counts[sample_id] += 1
						all_samples.add(sample_id)
				samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
				samples_with_any_copy = set([s[0] for s in sample_counts.items() if s[1] > 0])

				# check that hg is single-copy-core or just core
				if len(samples_with_single_copy.symmetric_difference(all_samples)) == 0:
					self.scc_homologs.add(hg)
				if len(samples_with_any_copy.symmetric_difference(all_samples)) == 0:
					self.core_homologs.add(hg)
					if self.hg_prop_multi_copy[hg] <= 0.05:
						self.specific_core_homologs.add(hg)

			sample_homolog_groups = defaultdict(set)
			samples = set([])
			hg_allele_representatives = set([])
			sample_allele_reads = defaultdict(lambda: defaultdict(int))
			sample_allele_unique_reads = defaultdict(lambda: defaultdict(int))
			sample_allele_novelty_reads = defaultdict(lambda: defaultdict(int))
			sample_allele_unique_novelty_reads = defaultdict(lambda: defaultdict(int))
			with open(paired_end_sequencing_file) as ossf:
				for line in ossf:
					sample = line.strip().split('\t')[0]
					result_file = snv_mining_outdir + sample + '.txt'
					samples.add(sample)
					if not os.path.isfile(result_file): continue
					with open(result_file) as orf:
						for i, hg_al in enumerate(orf):
							if i == 0: continue
							hg_al = hg_al.strip()
							hg, allele_representative, reads, reads_with_novelty, reads_uniquely_mapping, reads_uniquely_mapping_with_novelty = hg_al.split('\t')
							car = hg + '|' + allele_representative
							sample_homolog_groups[sample].add(hg)
							hg_allele_representatives.add(car)
							sample_allele_reads[sample][car] = int(reads)
							sample_allele_unique_reads[sample][car] = int(reads_uniquely_mapping)
							sample_allele_novelty_reads[sample][car] = int(reads_with_novelty)
							sample_allele_unique_novelty_reads[sample][car] = int(reads_uniquely_mapping_with_novelty)

			for s in samples:
				samp_hgs = sample_homolog_groups[s]
				print(self.core_homologs)
				if len(samp_hgs.intersection(self.core_homologs))/float(len(self.core_homologs)) < 0.9:
					self.avoid_samples.add(s)

			final_matrix_reads_file = outdir + 'Sample_by_OG_Allele_Read_Counts.matrix.txt'
			final_matrix_novelty_reads_file = outdir + 'Sample_by_OG_Allele_Novelty_Read_Counts.matrix.txt'
			final_matrix_unique_reads_file = outdir + 'Final_OG_Allele_Unique_Read_Counts.matrix.txt'
			final_matrix_unique_and_novelty_reads_file = outdir + 'Final_OG_Allele_Unique_and_Novelty_Read_Counts.matrix.txt'

			final_matrix_reads_handle = open(final_matrix_reads_file, 'w')
			final_matrix_novelty_reads_handle = open(final_matrix_novelty_reads_file, 'w')
			final_matrix_unique_reads_handle = open(final_matrix_unique_reads_file, 'w')
			final_matrix_unique_and_novelty_reads_handle = open(final_matrix_unique_and_novelty_reads_file, 'w')

			final_matrix_reads_handle.write(
				'\t'.join(['Sample/OG_Allele'] + list(sorted(hg_allele_representatives))) + '\n')
			final_matrix_novelty_reads_handle.write(
				'\t'.join(['Sample/OG_Allele'] + list(sorted(hg_allele_representatives))) + '\n')
			final_matrix_unique_reads_handle.write(
				'\t'.join(['Sample/OG_Allele'] + list(sorted(hg_allele_representatives))) + '\n')
			final_matrix_unique_and_novelty_reads_handle.write(
				'\t'.join(['Sample/OG_Allele'] + list(sorted(hg_allele_representatives))) + '\n')

			for s in samples:
				if s in self.avoid_samples: continue
				printlist_all = [s]
				printlist_uni = [s]
				printlist_nov = [s]
				printlist_uni_nov = [s]
				for h in hg_allele_representatives:
					printlist_all.append(str(sample_allele_reads[s][h]))
					printlist_uni.append(str(sample_allele_unique_reads[s][h]))
					printlist_nov.append(str(sample_allele_novelty_reads[s][h]))
					printlist_uni_nov.append(str(sample_allele_unique_novelty_reads[s][h]))
				final_matrix_reads_handle.write('\t'.join(printlist_all) + '\n')
				final_matrix_novelty_reads_handle.write('\t'.join(printlist_nov) + '\n')
				final_matrix_unique_reads_handle.write('\t'.join(printlist_uni) + '\n')
				final_matrix_unique_and_novelty_reads_handle.write('\t'.join(printlist_uni_nov) + '\n')

			final_matrix_reads_handle.close()
			final_matrix_novelty_reads_handle.close()
			final_matrix_unique_reads_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with generating matrices showcasing allele presence across samples.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def generateNoveltyReport(self, codon_alignment_file, snv_mining_outdir, outdir):
		try:
			novelty_report_file = outdir + 'Novelty_Report.txt'
			no_handle = open(novelty_report_file, 'w')
			no_handle.write('\t'.join(['gcf_id', 'sample', 'homolog_group', 'position_along_msa', 'alternate_allele',
									   'snv_count', 'reference_sample', 'reference_gene', 'reference_position',
									   'reference_allele']) + '\n')

			gene_pos_to_msa_pos = defaultdict(lambda: defaultdict(dict))
			msa_pos_alleles = defaultdict(lambda: defaultdict(set))
			with open(codon_alignment_file) as ocaf:
				for line in ocaf:
					line = line.strip()
					hg, cod_alignment = line.split('\t')
					with open(cod_alignment) as oca:
						for rec in SeqIO.parse(oca, 'fasta'):
							if hg == 'OG0001886':
								print(rec.id)
								print(str(rec.seq))
							sample_id, gene_id = rec.id.split('|')
							real_pos = 1
							for msa_pos, bp in enumerate(str(rec.seq)):
								if bp != '-':
									gene_pos_to_msa_pos[hg][gene_id][real_pos] = msa_pos+1
									real_pos += 1
									msa_pos_alleles[hg][msa_pos+1].add(bp.upper())

			for f in os.listdir(snv_mining_outdir):
				if not f.endswith('.snvs'): continue
				pe_sample = f.split('.snvs')[0]
				if pe_sample in self.avoid_samples: continue
				mges = set(['transp', 'integrase'])

				with open(snv_mining_outdir + f) as of:
					for i, line in enumerate(of):
						if i == 0: continue
						line = line.strip()
						ls = line.split('\t')
						snv_count = ls[1]
						gsh, ref_pos, ref_al, alt_al = ls[0].split('_|_')
						gene, sample, hg = gsh.split('|')
						if self.hg_prop_multi_copy[hg] >= 0.05: continue
						if any(word in self.comp_gene_info[gene]['product'].lower() for word in mges): continue
						ref_pos = int(ref_pos)+1

						#if self.comp_gene_info[gene]['direction'] == '+': ref_pos = ref_pos
						#else: ref_pos += 1
						print(ref_pos)
						print(line)
						print(self.comp_gene_info[gene])
						print(self.comp_gene_info[gene]['nucl_seq_with_flanks'][ref_pos+self.comp_gene_info[gene]['relative_start']-3:ref_pos+self.comp_gene_info[gene]['relative_start']+4])
						if not ref_pos in gene_pos_to_msa_pos[hg][gene]: continue
						msa_pos = gene_pos_to_msa_pos[hg][gene][int(ref_pos)]
						msa_pos_als = msa_pos_alleles[hg][msa_pos]
						print(msa_pos)
						print(msa_pos_alleles[hg][msa_pos - 2])
						print(msa_pos_alleles[hg][msa_pos - 1])
						print('=>\t' + str(msa_pos_als))
						print(msa_pos_alleles[hg][msa_pos + 1])
						print(msa_pos_alleles[hg][msa_pos + 2])
						print('\t'.join([pe_sample, hg, str(msa_pos), sample, gene, str(ref_pos), ref_al, alt_al, snv_count]))
						print('-'*80)
						assert (ref_al in msa_pos_alleles[hg][msa_pos])
						if not alt_al in msa_pos_als:
							no_handle.write('\t'.join([self.gcf_id, pe_sample, hg, str(msa_pos), alt_al, snv_count, sample, gene, str(ref_pos), ref_al]) + '\n')

			no_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues generating reports of alleles found per sample.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

def snv_miner(input_args):
	"""
	Function to mine for novel SNVs and identify alleles of homolog groups for GCF present in paired-end sequencing
	dataset.
	"""
	sample, bam_alignment, ref_fasta, hg_gene_to_rep, res_dir, bgc_hg_genes, comp_gene_info, logObject = input_args
	try:
		hg_rep_genes = defaultdict(set)
		for g, r in hg_gene_to_rep.items():
			hg_rep_genes[r].add(g)

		if not os.path.isfile(bam_alignment): return
		snvs_file = res_dir + sample + '.snvs'
		snv_outf = open(snvs_file, 'w')
		result_file = res_dir + sample + '.txt'
		outf = open(result_file, 'w')
		outf.write('\t'.join(['# hg', 'allele_representative', 'reads', 'reads_with_novelty', 'reads_uniquely_mapping',
							  'reads_uniquely_mapping_with_novelty']) + '\n')

		bam_handle = pysam.AlignmentFile(bam_alignment, 'rb')

		topaligns_file = res_dir + sample + '_topaligns.bam'
		topaligns_file_sorted = res_dir + sample + '_topaligns.sorted.bam'
		topaligns_handle = pysam.AlignmentFile(topaligns_file, "wb", template=bam_handle)

		unialigns_file = res_dir + sample + '_unialigns.bam'
		unialigns_file_sorted = res_dir + sample + '_unialigns.sorted.bam'
		unialigns_handle = pysam.AlignmentFile(unialigns_file, "wb", template=bam_handle)

		for hg, hg_genes in bgc_hg_genes.items():
			read_ascores_per_allele = defaultdict(list)
			read_genes_mapped = defaultdict(set)
			snv_counts = defaultdict(set)
			hg_genes_covered = 0
			rep_alignments = defaultdict(lambda: defaultdict(set))
			with open(ref_fasta) as opff:
				for rec in SeqIO.parse(opff, 'fasta'):
					if rec.id.split('|')[-1] != hg: continue
					g, sample, _ = rec.id.split('|')
					ginfo = comp_gene_info[g]

					gstart = ginfo['relative_start']
					gend = ginfo['relative_end']
					offset = gstart

					gene_length = gend - gstart

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
					if gene_coverage_1 < 0.90: continue
					hg_genes_covered += 1

					for read1_alignment, read2_alignment in util.read_pair_generator(bam_handle, region_string=rec.id,
																				start=gstart, stop=gend+1):
						if read1_alignment is None or read2_alignment is None: continue
						read_name = read1_alignment.query_name
						read1_ascore = read1_alignment.tags[0][1]
						read2_ascore = read2_alignment.tags[0][1]
						combined_ascore = read1_ascore + read2_ascore

						snvs = set([])
						g_rep = hg_gene_to_rep[rec.id]

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

						n_found = False
						for b in read1_alignment.get_aligned_pairs(with_seq=True):
							if b[0] == None or b[1] == None: continue
							if not b[2].islower(): continue
							ref_pos = b[1]
							que_qual = read1_queryqua[b[0]]
							alt_al = read1_queryseq[b[0]].upper()
							ref_al = read1_referseq[b[1] - min_read1_ref_pos].upper()
							if b[2] == 'n' or ref_al == 'N' or alt_al == 'N': n_found = True; break
							assert (ref_al == str(rec.seq).upper()[b[1]])
							assert (alt_al != ref_al)
							if que_qual >= 30 and not alignment_has_indel and mismatch_count <= 5 and read_overlap_prop <= 0.25 and min_read_length >= 75:
								snvs.add(str(rec.id) + '_|_' + str(ref_pos - offset) + '_|_' + ref_al + '_|_' + alt_al)
								snv_counts[
									str(rec.id) + '_|_' + str(ref_pos - offset) + '_|_' + ref_al + '_|_' + alt_al].add(
									read_name)

						for b in read2_alignment.get_aligned_pairs(with_seq=True):
							if b[0] == None or b[1] == None: continue
							if not b[2].islower(): continue
							ref_pos = b[1]
							que_qual = read2_queryqua[b[0]]
							alt_al = read2_queryseq[b[0]].upper()
							ref_al = read2_referseq[b[1] - min_read2_ref_pos].upper()
							if b[2] == 'n' or ref_al == 'N' or alt_al == 'N': n_found = True; break
							assert (ref_al == str(rec.seq).upper()[b[1]])
							assert (alt_al != ref_al)
							if que_qual >= 30 and not alignment_has_indel and mismatch_count <= 5 and read_overlap_prop <= 0.25 and min_read_length >= 75:
								snvs.add(str(rec.id) + '_|_' + str(ref_pos - offset) + '_|_' + ref_al + '_|_' + alt_al)
								snv_counts[
									str(rec.id) + '_|_' + str(ref_pos - offset) + '_|_' + ref_al + '_|_' + alt_al].add(
									read_name)

						if n_found: continue
						read_genes_mapped[read1_alignment.query_name].add(rec.id)
						rep_alignments[rec.id][read1_alignment.query_name].add(tuple([read1_alignment, read2_alignment]))
						read_ascores_per_allele[read1_alignment.query_name].append(
							[g_rep.split('|')[0], g_rep.split('|')[1], combined_ascore, snvs, g])

			#if hg_genes_covered / float(len(hg_genes)) < 0.80: continue

			supported_snvs = set([])
			allele_reads = defaultdict(set)
			allele_reads_with_mismatch = defaultdict(set)
			multi_partitioned_reads = set([])
			for read in read_ascores_per_allele:
				top_score = -1000000
				top_score_grep = None
				score_sorted_alignments = sorted(read_ascores_per_allele[read], key=itemgetter(2), reverse=True)
				for i, align in enumerate(score_sorted_alignments):
					g_rep = align[0] + '|' + align[1] + '|' + hg
					if i == 0: top_score = align[2]; top_score_grep = g_rep
					if (i == 0 and align[2] == top_score) and (
							len(score_sorted_alignments) == 1 or align[2] > score_sorted_alignments[i + 1][2]):
						for snv in align[3]:
							if len(snv_counts[snv]) >= 5:
								supported_snvs.add(snv)
					if align[2] == top_score:
						#g_map = read_genes_mapped[read].intersection(hg_rep_genes[g_rep])
						#g_map_prop = len(g_map) / float(len(hg_rep_genes[g_rep]))
						#if g_map_prop < 0.80: continue

						allele_reads[g_rep].add(read)
						for snv in align[3]:
							if len(snv_counts[snv]) >= 5:
								allele_reads_with_mismatch[g_rep].add(read)
								break

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
					'\t'.join([str(x) for x in [hg, al, len(allele_reads[al]), len(allele_reads_with_mismatch[al]),
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
	except Exception as e:
		if logObject:
			logObject.error("Issues with mining for SNVs/parsing alignment file.")
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def create_hmm_profiles(inputs):
	"""

	"""
	hg, sample_sequences, prot_seq_dir, prot_alg_dir, prot_hmm_dir, logObject = inputs

	hg_prot_fasta = prot_seq_dir + '/' + hg + '.faa'
	hg_prot_msa = prot_alg_dir + '/' + hg + '.msa.faa'
	hg_prot_hmm = prot_hmm_dir + '/' + hg + '.hmm'

	hg_prot_handle = open(hg_prot_fasta, 'w')
	for s in sample_sequences:
		hg_prot_handle.write('>' + s + '\n' + str(sample_sequences[s]) + '\n')
	hg_prot_handle.close()

	mafft_cmd = ['mafft', '--maxiterate', '1000', '--localpair', hg_prot_fasta, '>', hg_prot_msa]
	if logObject:
		logObject.info('Running mafft with the following command: %s' % ' '.join(mafft_cmd))
	try:
		subprocess.call(' '.join(mafft_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(mafft_cmd))
	except:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(mafft_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(mafft_cmd))

	hmmbuild_cmd = ['hmmbuild', '--amino', '-n', hg, hg_prot_hmm, hg_prot_msa]
	if logObject:
		logObject.info('Running hmmbuild (from HMMER3) with the following command: %s' % ' '.join(hmmbuild_cmd))
	try:
		subprocess.call(' '.join(hmmbuild_cmd), shell=True, stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(hmmbuild_cmd))
	except:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(hmmbuild_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(hmmbuild_cmd))

	if logObject:
		logObject.info('Constructed profile HMM for homolog group %s' % hg)

def popgen_analysis_of_hg(inputs):
	"""
	Helper function which is to be called from the runPopulationGeneticsAnalysis() function to parallelize population
	genetics analysis of each homolog group.

	:param inputs: list of inputs passed in by GCF.runPopulationGeneticsAnalysis().
	"""
	gcf_id, hg, codon_alignment_fasta, popgen_dir, plots_dir, comp_gene_info, hg_genes, bgc_sample, hg_prop_multi_copy, hg_order_scores, sample_population, logObject = inputs
	domain_plot_file = plots_dir + hg + '_domain.txt'
	position_plot_file = plots_dir + hg + '_position.txt'
	popgen_plot_file = plots_dir + hg + '_popgen.txt'
	plot_pdf_file = plots_dir + hg + '.pdf'

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
	position_plot_handle.write(
		'\t'.join(['pos', 'num_seqs', 'num_alleles', 'num_gaps', 'maj_allele_freq']) + '\n')
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
	for gene in hg_genes[hg]:
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

	hg_info = [gcf_id, hg, '; '.join(products), hg_order_scores[hg], hg_prop_multi_copy[hg],
				median_gene_length, is_core, len(seqs), prop_samples_with_hg, tajimas_d, total_core_codons,
				total_variable_codons, nonsynonymous_sites, synonymous_sites, dn_ds, '; '.join(all_domains)]

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

		pop_prop_with_hg = defaultdict(float)
		pops_with_hg = 0
		populations_order = []
		within_population_differences = []
		for pop in population_samples:
			pop_prop_with_hg[pop] = len(population_samples[pop].intersection(samples)) / len(
				population_samples[pop])
			if pop_prop_with_hg[pop] > 0: pops_with_hg += 1
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
		hg_population_info = [pops_with_hg,
							   '|'.join([str(x[0]) + '=' + str(x[1]) for x in pop_prop_with_hg.items()]),
							   anova_pval]
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
