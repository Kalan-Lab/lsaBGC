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

FLANK_SIZE = 500

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