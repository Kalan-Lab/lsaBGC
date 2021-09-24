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
from Bio.codonalign.codonalignment import _get_codon2codon_matrix, _get_subgraph, _count_replacement, _G_test

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
		raise RuntimeError(
			"CodonAlignment object for mktest should be of equal length."
		)
	codon_num = codon_aln_len[0] // 3
	# prepare codon_dict (taking stop codon as an extra amino acid)
	codon_dict = copy.deepcopy(codon_table.forward_table)
	for stop in codon_table.stop_codons:
		codon_dict[stop] = "stop"
	# prepare codon_lst
	codon_lst = []
	for codon_aln in codon_alns:
		codon_lst.append([])
		for i in codon_aln:
			codon_lst[-1].append(_get_codon_list(i.seq))
	codon_set = []
	for i in range(codon_num):
		uniq_codons = []
		for j in codon_lst:
			uniq_codon = {k[i] for k in j}
			uniq_codons.append(uniq_codon)
		codon_set.append(uniq_codons)
	syn_fix, nonsyn_fix, syn_poly, nonsyn_poly = 0, 0, 0, 0
	G, nonsyn_G = _get_codon2codon_matrix(codon_table=codon_table)
	for i in codon_set:
		all_codon = set(i[0].union(*i[1:]))
		all_codon = all_codon.difference(set(['---']))
		if len(all_codon) <= 1: continue
		fix_or_not = all(len(k) == 1 for k in i)
		if fix_or_not:
			# fixed
			nonsyn_subgraph = _get_subgraph(all_codon, nonsyn_G)
			subgraph = _get_subgraph(all_codon, G)
			this_non = _count_replacement(all_codon, nonsyn_subgraph)
			this_syn = _count_replacement(all_codon, subgraph) - this_non
			nonsyn_fix += this_non
			syn_fix += this_syn
		else:
			# not fixed
			all_codon = i[0].difference(set(['---']))
			if len(all_codon) <= 1: continue
			nonsyn_subgraph = _get_subgraph(all_codon, nonsyn_G)
			subgraph = _get_subgraph(all_codon, G)
			this_non = _count_replacement(all_codon, nonsyn_subgraph)
			this_syn = _count_replacement(all_codon, subgraph) - this_non
			nonsyn_poly += this_non
			syn_poly += this_syn
	pval = "NA"
	try:
		pval = _G_test([syn_fix, nonsyn_fix, syn_poly, nonsyn_poly])
	except:
		pass
	print([pval, syn_fix, nonsyn_fix, syn_poly, nonsyn_poly])
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
					sample_seqs[hg][sample].append(str(rec.seq).upper().replace('N', '-'))
		return sample_seqs
	except:
		sys.stderr.write("Problem parsing codon alignment. Exiting now ...")
		raise RuntimeError


def node_mktesting(node_id, sample_seqs, all_children, all_tree_samples):
	node_hg_info = []
	pvalues = []

	for hg in sample_seqs:
		print(hg)
		node_cod_seqs = []
		other_cod_seqs = []
		for sample in all_children:
			for seq in sample_seqs[hg][sample]:
				node_cod_seqs.append(SeqRecord(CodonSeq(seq)))
		for sample in all_tree_samples.difference(all_children):
			for seq in sample_seqs[hg][sample]:
				other_cod_seqs.append(SeqRecord(CodonSeq(seq)))

		if len(node_cod_seqs) >= 3 and len(other_cod_seqs) >= 3:
			node_cod_alg_obj = CodonAlignment(node_cod_seqs)
			other_cod_alg_obj = CodonAlignment(other_cod_seqs)
			list_of_cod_algns = [node_cod_alg_obj, other_cod_alg_obj]

			pval, syn_fix, nonsyn_fix, syn_poly, nonsyn_poly = mktest(list_of_cod_algns)
			pvalues.append(pval)
			node_hg_info.append([node_id, hg, '; '.join(all_children), syn_fix, nonsyn_fix, syn_poly, nonsyn_poly])

	return([node_hg_info, pvalues])

LEAF_NAMES = set([])

def is_innernode(i):
	try:
		assert(not i  in LEAF_NAMES)
		return True
	except:
		return False


def recursively_get_children(direct_children_map, curr_node):
	""" feeling like fibonacci """
	children = set([])
	direct_children = direct_children_map[curr_node]
	for child in direct_children:
		if is_innernode(child):
			children = children.union(recursively_get_children(direct_children_map, child))
		else:
			children.add(child)
	return children


def parse_phylogeny(tree):
	try:
		t = Tree(tree)
		counter = 1
		for node in t.traverse("postorder"):
			if not node.is_leaf() and node.name.strip() == '':
				node.name = str(counter)
				counter += 1

		direct_children = defaultdict(set)
		for node in t.traverse("postorder"):
			try:
				parent_name = node.up.name
				direct_children[parent_name].add(node.name)
			except: pass
		return direct_children
	except:
		sys.stderr.write("Problem parsing phylogeny! Please check input newick file. Exiting now ...")
		raise RuntimeError


def crawlingMKTest(tree, gcf_id, codon_alignment_file, output):
	try:
		assert (os.path.isfile(tree) and os.path.isfile(codon_alignment_file))
	except:
		sys.stderr.write("Either phylogeny or codon alignments listing file does not exist. Exiting now ..."); raise RuntimeError

	direct_children = parse_phylogeny(tree)
	#print(direct_children)
	sample_seqs = read_codalign_file(codon_alignment_file)
	print(sample_seqs['OG0000157'])
	output = os.path.abspath(output)

	try:
		assert (not os.path.isfile(output))
	except:
		sys.stderr.write("Output file already exists. Please remove/rename."); raise RuntimeError

	out = open(output, 'w')
	out.write('hg\tnode\tchildren\tadj_pvalue\tsyn_fix\tnonsyn_fix\tsyn_poly\tnonsyn_poly\n')

	all_pvalues = []
	all_node_hgs = []

	all_tree_samples = set([])
	for leaf in Tree(tree):
		all_tree_samples.add(str(leaf).strip('\n').lstrip('-'))
	global LEAF_NAMES
	LEAF_NAMES = all_tree_samples

	#print(all_tree_samples)
	for par in direct_children:
		#print(par)
		all_children = recursively_get_children(direct_children, par)
		if len(all_children) >= 5:
			#print('----------------')
			#print(all_children)
			node_hgs, pvalues = node_mktesting(par, sample_seqs, all_children, all_tree_samples)
			print(pvalues)
			print(node_hgs)
			all_node_hgs += node_hgs
			all_pvalues += pvalues

	adj_pvalues = p_adjust_bh(all_pvalues)
	for i, data in enumerate(all_node_hgs):
		if adj_pvalues[i] < 0.05:
			out.write('\t'.join([str(x) for x in ([gcf_id, data[1], data[0]] + data[2:])]) + '\n')
	out.close()

if __name__ == '__main__':
	# Pull out the arguments.
	parser = argparse.ArgumentParser(description=""" This program crawls up a phylogenetic tree and runs the 
	McDonald-Kreitman Test between the children of each innernode to the remaining sequneces/samples.""")

	parser.add_argument('-t', '--tree', help='Phylogenetic tree in Newick format. Inner nodes must be named!', required=True)
	parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
	parser.add_argument('-a', '--codon_alignments', help="File listing the codon alignments for each homolog group in the GCF. Can be found as part of PopGene output.", required=True)
	parser.add_argument('-o', '--output', help="Output iTol dataset file.", required=True)

	args = parser.parse_args()

	crawlingMKTest(args.tree, args.gcf_id, args.codon_alignments, args.output)