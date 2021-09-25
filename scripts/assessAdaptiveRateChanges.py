#!/usr/bin/env python

import os
import sys
import argparse
from scipy import stats
from collections import defaultdict
from ete3 import Tree
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.codonalign.codonalignment import CodonAlignment
from Bio import SeqIO

from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio import BiopythonWarning
from Bio.codonalign.codonseq import _get_codon_list, CodonSeq, cal_dn_ds
from Bio.codonalign.chisq import chisqprob
from Bio.codonalign.codonalignment import _get_codon2codon_matrix, _get_subgraph, _count_replacement
from scipy.stats import chisquare
GCFID = None
def mktest(codon_alns, codon_table=None):
	"""McDonald-Kreitman test for neutrality.
	Implement the McDonald-Kreitman test for neutrality (PMID: 1904993)
	This method counts changes rather than sites
	(http://mkt.uab.es/mkt/help_mkt.asp).
	Arguments:
	 - codon_alns  - list of CodonAlignment to compare (each
	   CodonAlignment object corresponds to gene sampled from a species)
	Return the p-value of test result.
	"""
	import copy

	if codon_table is None:
		codon_table = CodonTable.generic_by_id[1]
	if not all(isinstance(i, CodonAlignment) for i in codon_alns):
		raise TypeError("mktest accepts CodonAlignment list.")
	codon_aln_len = [i.get_alignment_length() for i in codon_alns]
	if len(set(codon_aln_len)) != 1:
		raise RuntimeError( "CodonAlignment object for mktest should be of equal length." )
	codon_num = codon_aln_len[0] // 3

	# prepare codon_dict (taking stop codon as an extra amino acid)
	codon_dict = copy.deepcopy(codon_table.forward_table)
	for stop in codon_table.stop_codons:
		codon_dict[stop] = "stop"

	# prepare codon_lst
	codon_lst = []
	codpos_allele_count = defaultdict(lambda: defaultdict(int))
	for codon_aln in codon_alns:
		codon_lst.append([])
		for i in codon_aln:
			cod_lst = [str(i.seq)[p:p + 3] for p in range(0, len(str(i.seq)), 3)]
			codon_lst[-1].append(cod_lst)
			for p in range(codon_num):
				codpos_allele_count[p][cod_lst[p]] += 1

	codon_set = []
	for i in range(codon_num):
		gap_count = codpos_allele_count[i]['---']
		tot_count = sum(codpos_allele_count[i].values())
		assert(tot_count > 0)
		gap_prop = float(gap_count)/float(tot_count)
		if gap_prop >= 0.1: continue

		uniq_codons = []
		for j in codon_lst:
			uniq_codon = set([k[i] for k in j])
			uniq_codon_filt = set([])
			for cod in uniq_codon:
				if not '-' in cod and codpos_allele_count[i][cod] >= 2:
					uniq_codon_filt.add(cod)
			uniq_codons.append(uniq_codon_filt)
		codon_set.append(uniq_codons)

	syn_fix, nonsyn_fix, syn_poly, nonsyn_poly = 0, 0, 0, 0
	G, nonsyn_G = _get_codon2codon_matrix(codon_table=codon_table)
	for i in codon_set:
		#if len(i[0]) == 0 or len(i[1]) == 0: continue
		all_codon = sorted(i[0].union(i[1]))
		if len(all_codon) == 1: continue

		fix_or_not = ((len(i[0]) == 1) and (len(i[1]) == 1))
		if fix_or_not:
			# B and A fixed on different alleles
			all_codon = [i[1][0], i[0][0]]
			nonsyn_subgraph = _get_subgraph(all_codon, nonsyn_G)
			subgraph = _get_subgraph(all_codon, G)
			this_non = _count_replacement(all_codon, nonsyn_subgraph)
			this_syn = _count_replacement(all_codon, subgraph) - this_non
			nonsyn_fix += this_non
			syn_fix += this_syn
		elif len(i[1]) == 1 and len(i[0]) > 1 and list(i[1])[0] in i[0] and len(i[0]) == 2:
			# only species B fixed, species A is bi-allelic at codon site with one allele matching the fixed
			# allele in species B.
			all_codon = [i[1][0], [x for x in i[0] if x != i[1][0]][0]]
			nonsyn_subgraph = _get_subgraph(all_codon, nonsyn_G)
			subgraph = _get_subgraph(all_codon, G)
			this_non = _count_replacement(all_codon, nonsyn_subgraph)
			this_syn = _count_replacement(all_codon, subgraph) - this_non
			nonsyn_poly += this_non
			syn_poly += this_syn

	pval = np.nan
	if syn_fix >= 5 and nonsyn_fix >= 5 and syn_poly >= 5 and nonsyn_poly >= 5:
		try:
			chisquare2, pval = chisquare([syn_fix, nonsyn_fix, syn_poly, nonsyn_poly])
		except:
			pass
	return([pval, syn_fix, nonsyn_fix, syn_poly, nonsyn_poly])

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

def read_codalign_file(codon_alignment_listing_file):
	try:
		sample_seqs = defaultdict(lambda: defaultdict(list))
		for i, hg_line in enumerate(open(codon_alignment_listing_file)):
			hg, hg_codon_alignment_file = hg_line.strip().split('\t')

			with open(hg_codon_alignment_file) as ohcaf:
				for rec in SeqIO.parse(ohcaf, 'fasta'):
					sample = rec.id.split('|')[0]
					sample_seqs[hg][sample].append(str(rec.seq).strip().upper().replace('N', '-'))
		return sample_seqs
	except:
		sys.stderr.write("Problem parsing codon alignment. Exiting now ...")
		raise RuntimeError


def comp_species(sp1, sp2, skin_species_samples, sample_seqs):
	comp_hg_info = []
	pvalues = []

	for hg in sample_seqs:
		sp1_seqs = []
		sp2_seqs = []
		sp1_samples_with_hg = set([])
		sp2_samples_with_hg = set([])

		for sample in skin_species_samples[sp1]:
			for seq in sample_seqs[hg][sample]:
				sp1_samples_with_hg.add(sample)
				sp1_seqs.append(SeqRecord(CodonSeq(seq)))

		for sample in skin_species_samples[sp2]:
			for seq in sample_seqs[hg][sample]:
				sp2_samples_with_hg.add(sample)
				sp2_seqs.append(SeqRecord(CodonSeq(seq)))

		if len(skin_species_samples[sp1]) == 0 or len(skin_species_samples[sp2]) == 0: continue
		sp1_prop_with_hg = len(sp1_samples_with_hg)/float(len(skin_species_samples[sp1]))
		sp2_prop_with_hg = len(sp2_samples_with_hg)/float(len(skin_species_samples[sp2]))

		if len(sp1_samples_with_hg) >= 3 and len(sp2_samples_with_hg) >= 3 and sp1_prop_with_hg >= 0.25 and sp2_prop_with_hg >= 0.25:
			print(sp1 + '\t' + sp2 + '\t' + hg)
			sp1_cod_alg_obj = CodonAlignment(sp1_seqs)
			sp2_cod_alg_obj = CodonAlignment(sp2_seqs)
			list_of_cod_algns = [sp1_cod_alg_obj, sp2_cod_alg_obj]
			pval, syn_fix, nonsyn_fix, syn_poly, nonsyn_poly = mktest(list_of_cod_algns)
			pvalues.append(pval)
			comp_hg_info.append([hg, sp1, sp2, sp1_prop_with_hg, sp2_prop_with_hg, syn_fix, nonsyn_fix, syn_poly, nonsyn_poly])
			print([pval, syn_fix, nonsyn_fix, syn_poly, nonsyn_poly])
	return([comp_hg_info, pvalues])

def speciesComparisonMKTest(skin_associated, gcf_id, codon_alignment_file, output):
	global GCFID
	GCFID = gcf_id
	try:
		assert (os.path.isfile(skin_associated) and os.path.isfile(codon_alignment_file))
	except:
		sys.stderr.write("An input file does not exist. Exiting now ..."); raise RuntimeError

	sample_seqs = read_codalign_file(codon_alignment_file)
	output = os.path.abspath(output)

	try:
		assert (not os.path.isfile(output))
	except:
		sys.stderr.write("Output file already exists. Please remove/rename."); raise RuntimeError

	out = open(output, 'w')
	out.write('gcf_id\thg\tspecies_1\tspecies_2\tpval\tadj_pvalue\tprop_sp1_with_hg\tprop_sp2_with_hg\tsyn_fix\tnonsyn_fix\tsyn_poly\tnonsyn_poly\n')

	all_pvalues = []
	all_comp_hgs = []

	skin_species = set([])
	with open(skin_associated) as osa:
		for line in osa:
			line = line.strip()
			skin_species.add(line)

	skin_species_samples = defaultdict(set)
	with open(codon_alignment_file) as ocaf:
		for line in ocaf:
			line = line.strip()
			hg, hg_cod_alg_file = line.split('\t')
			with open(hg_cod_alg_file) as ohcaf:
				for rec in SeqIO.parse(ohcaf, 'fasta'):
					sample = rec.id.split('|')[0]
					spec = ' '.join(sample.split('_')[:2])
					if spec in skin_species:
						skin_species_samples[spec].add(sample)

	for i, sp1 in enumerate(sorted(skin_species)):
		for j, sp2 in enumerate(sorted(skin_species)):
			if sp1 == sp2: continue
			comp_hg_info, pvalues = comp_species(sp1, sp2, skin_species_samples, sample_seqs)
			all_comp_hgs += comp_hg_info
			all_pvalues += pvalues

	adj_pvalues = p_adjust_bh(all_pvalues)

	for i, data in enumerate(all_comp_hgs):
		adj_pvalue = adj_pvalues[i]
		sf = float(data[-4])
		nsf = float(data[-3])
		sp = float(data[-2])
		nsp = float(data[-1])
		rate_fold_change = np.nan
		if sf > 0 and sp > 0:
			dn_ds = nsf/sf
			pn_ps = nsp/sp
			if dn_ds > 0:
				rate_fold_change = pn_ps/dn_ds
			else:
				rate_fold_change = 'infinity'
		out.write('\t'.join([str(x) for x in ([gcf_id] + data[:3] + [all_pvalues[i], adj_pvalue, rate_fold_change] + data[3:])]) + '\n')
	out.close()

if __name__ == '__main__':
	# Pull out the arguments.
	parser = argparse.ArgumentParser(description="""This program assesses directional selection based on the 
	McDonald-Kreitman test.""")

	parser.add_argument('-s', '--skin_species', help='List of species associated with the skin.', required=True)
	parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
	parser.add_argument('-a', '--codon_alignments', help="File listing the codon alignments for each homolog group in the GCF. Can be found as part of PopGene output.", required=True)
	parser.add_argument('-o', '--output', help="Output iTol dataset file.", required=True)

	args = parser.parse_args()
	speciesComparisonMKTest(args.skin_species, args.gcf_id, args.codon_alignments, args.output)