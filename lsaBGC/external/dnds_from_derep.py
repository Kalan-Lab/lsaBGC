#!/usr/bin/env python

# Version 0.2.2 - 7.14.20
# Tried to add handling of empty files
# Fixed stupid nsmiscan version

# Version 0.2.1 - 7.14.20
# Updated the substitution_matrix that seems to be different now?
# Added some additional debugging

# Version 0.1.3
# Matt Olm
# mattolm@berkeley.edu
# 02.20.19
# https://biotite.berkeley.edu/j/user/mattolm/notebooks/OtherProjects/RefSeq/goANI_2_development_3_dRep.ipynb

# CHANGED IN VERSION 0.1.3:
# Fixed a too many alignments bug with a workaround https://biotite.berkeley.edu/j/user/mattolm/notebooks/OtherProjects/RefSeq/goANI_2_development_4_debugging.ipynb

# CHANGED IN VERSION 0.2.0:
# updated to match https://biotite.berkeley.edu/j/user/mattolm/notebooks/OtherProjects/RefSeq/_SupplementalNotebooks_1_dnds_HGT_calculations.ipynb


"""
__modification_other__ = "Rauf Salamzade"

__author__ = "Matt Olm"
__version__ = "0.2.2"
__license__ = "MIT"

Copyright (c) <year> <copyright holders>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


"""
Modifications made for lsaBGC framework, as well as removing of unused functions for simplicity.
"""

import os
import sys
import glob
import scipy
import pickle
import shutil
import argparse
import textwrap
import datetime
import traceback
import numpy as np
import pandas as pd

import warnings
warnings.simplefilter("ignore")

import numpy as np
from math import log
import concurrent.futures
from itertools import permutations
from collections import defaultdict

import Bio
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SubsMat import MatrixInfo
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.codonalign.codonalphabet import default_codon_table

from Bio.Align import substitution_matrices


def goANI_dnds_calculation(fna1, faa1, fna2, faa2, gedb, debug=False):
    '''
    This is a threadable command to determine the dn/ds of two genomes
    based on a list of genes

    Arguments:
        fna1 : .fna file of genome1
        faa1 : .faa file of genome1
        fna2 : .fna file of genome2
        faa2 : .faa file of genome2
        gedb : datatable listing the genes to align and calculate dn/ds for

    Returns:
        dndb : data-table containing raw dn/ds information
    '''
    # load .fasta files
    g1n = SeqIO.to_dict(SeqIO.parse(fna1, 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA()))
    g1a = SeqIO.to_dict(SeqIO.parse(faa1, 'fasta', alphabet=IUPAC.protein))
    g2n = SeqIO.to_dict(SeqIO.parse(fna2, 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA()))
    g2a = SeqIO.to_dict(SeqIO.parse(faa2, 'fasta', alphabet=IUPAC.protein))

    # set up aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    #print(MatrixInfo.blosum62)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -12
    aligner.extend_gap_score = -3

    # set up table
    table = defaultdict(list)

    # for every gene-pair to align
    j = 0
    for i, row in gedb.iterrows():
        try:
            # get the sequences
            a1 = g1a[row['qry_id']]
            if a1[-1] == '*':
                a1 = a1[:-1]
            a2 = g2a[row['sbj_id']]
            if a2[-1] == '*':
                a2 = a2[:-1]

            s1 = g1n[row['qry_id']]
            s2 = g2n[row['sbj_id']]

            # alingn them
            alignments = aligner.align(a1.seq, a2.seq)

            # Arbitrary cutoff to make sure this doesn't bug out
            if len(alignments) > 10000:
                print("Ahhh! {0} vs {1} has {2} alignments".format(row['qry_id'], row['sbj_id'],len(alignments)))
                raise Exception('Too many alignments exception')

            # convert to multi-sequence alignment
            alignment = min(alignments)
            ass = str(alignment).split('\n')
            msa = MultipleSeqAlignment([SeqRecord(Seq(ass[0], alphabet=IUPAC.protein)),
                                       SeqRecord(Seq(ass[-2], alphabet=IUPAC.protein))])

            # convert to codon alignment
            codon_aln = Bio.codonalign.build(msa, [s1, s2])

            # calculate dn/ds on the codon alignment
            dS, S, dN, N = custom_dn_ds(codon_aln._records[0], codon_aln._records[1])

            # save
            table['qry_id'].append(row['qry_id'])
            table['sbj_id'].append(row['sbj_id'])
            table['S_changed'].append(dS)
            table['S_sites'].append(S)
            table['N_changed'].append(dN)
            table['N_sites'].append(N)

            j += 1
            if debug:
                if j >= 10:
                    break

        except Exception as e:
            print("Alignment exception- {0}".format(e))
            table['qry_id'].append(row['qry_id'])
            table['sbj_id'].append(row['sbj_id'])
            table['S_changed'].append(0)
            table['S_sites'].append(0)
            table['N_changed'].append(0)
            table['N_sites'].append(0)

    dnDb = pd.DataFrame(table)

    return dnDb

def custom_dn_ds(codon_seq1, codon_seq2, method="NG86", codon_table=default_codon_table, k=1, cfreq=None):
    """
    http://biopython.org/DIST/docs/api/Bio.codonalign.codonseq-pysrc.html#cal_dn_ds

    Calculate dN and dS of the given two sequences.

    Available methods:
      - NG86  - `Nei and Gojobori (1986)`_ (PMID 3444411).
      - LWL85 - `Li et al. (1985)`_ (PMID 3916709).
      - ML    - `Goldman and Yang (1994)`_ (PMID 7968486).
      - YN00  - `Yang and Nielsen (2000)`_ (PMID 10666704).

    .. _`Nei and Gojobori (1986)`: http://www.ncbi.nlm.nih.gov/pubmed/3444411
    .. _`Li et al. (1985)`: http://www.ncbi.nlm.nih.gov/pubmed/3916709
    .. _`Goldman and Yang (1994)`: http://mbe.oxfordjournals.org/content/11/5/725
    .. _`Yang and Nielsen (2000)`: https://doi.org/10.1093/oxfordjournals.molbev.a026236

    Arguments:
    - codon_seq1 - CodonSeq or or SeqRecord that contains a CodonSeq
    - codon_seq2 - CodonSeq or or SeqRecord that contains a CodonSeq
    - w  - transition/transversion ratio
    - cfreq - Current codon frequency vector can only be specified
     when you are using ML method. Possible ways of
     getting cfreq are: F1x4, F3x4 and F61.

    """
    # Convert to codonseq
    if isinstance(codon_seq1, Bio.codonalign.codonseq.CodonSeq) \
            and isinstance(codon_seq2, Bio.codonalign.codonseq.CodonSeq):
        pass
    elif isinstance(codon_seq1, SeqRecord) and isinstance(codon_seq2, SeqRecord):
        codon_seq1 = codon_seq1.seq
        codon_seq2 = codon_seq2.seq
    else:
        raise TypeError("cal_dn_ds accepts two CodonSeq objects or SeqRecord "
                      "that contains CodonSeq as its seq!")
    if len(codon_seq1.get_full_rf_table()) != len(codon_seq2.get_full_rf_table()):
        raise RuntimeError("full_rf_table length of seq1 ({0}) and seq2 ({1}) "
                         "are not the same".format(
                             len(codon_seq1.get_full_rf_table()),
                             len(codon_seq2.get_full_rf_table()))
                         )
    if cfreq is None:
        cfreq = 'F3x4'
    elif cfreq is not None and method != 'ML':
        raise RuntimeError("cfreq can only be specified when you "
                         "are using ML method")
    if cfreq not in ('F1x4', 'F3x4', 'F61'):
        import warnings
        warnings.warn("Unknown cfreq ({0}). Only F1x4, F3x4 and F61 are "
                    "acceptable. Use F3x4 in the following.".format(cfreq))
        cfreq = 'F3x4'

    # get codon lists
    seq1_codon_lst = get_codon_list(codon_seq1)
    seq2_codon_lst = get_codon_list(codon_seq2)

    # filter for stop codons
    #if ('TGA' in seq1_codon_lst) or ('TGA' in seq2_codon_lst):
    #    return 0,0,0,0

    # remove gaps in seq_codon_lst
    seq1 = []
    seq2 = []
    for i, j in zip(seq1_codon_lst, seq2_codon_lst):
        if ('-' not in i) and ('-' not in j):
            seq1.append(i)
            seq2.append(j)

    dS, S, dN, N = _ng86_custom(seq1, seq2, k, codon_table)
    return dS, S, dN, N

def _ng86_custom(seq1, seq2, k, codon_table):
    """
    NG86 method main function (PRIVATE).
    http://biopython.org/DIST/docs/api/Bio.codonalign.codonseq-pysrc.html#cal_dn_ds

    MATT - You changed this to return details instead of dN, dS.

    It now returns changed S, total S, changed N, total N

    """
    S_sites1, N_sites1 = count_sites(seq1,
                                            codon_table=codon_table, k=k)
    S_sites2, N_sites2 = count_sites(seq2,
                                            codon_table=codon_table, k=k)
    S_sites = (S_sites1 + S_sites2) / 2.0
    N_sites = (N_sites1 + N_sites2) / 2.0
    SN = [0, 0]
    for i, j in zip(seq1, seq2):
        SN = [m + n for m, n in zip(SN, _count_diff_NG86(i, j,
                                                           codon_table=codon_table))]
    #print(SN)
    ps = SN[0] / S_sites
    pn = SN[1] / N_sites
    if ps < 3 / 4:
        dS = abs(-3.0 / 4 * log(1 - 4.0 / 3 * ps))
    else:
        dS = -1
    if pn < 3 / 4:
        dN = abs(-3.0 / 4 * log(1 - 4.0 / 3 * pn))
    else:
        dN = -1
    return SN[0], S_sites, SN[1], N_sites


def get_codon_list(codonseq):
    """List of codons according to full_rf_table for counting (PRIVATE)."""
    full_rf_table = codonseq.get_full_rf_table()
    codon_lst = []
    for i, k in enumerate(full_rf_table):
        if isinstance(k, int):
            start = k
            try:
                end = int(full_rf_table[i+1])
            except IndexError:
                end = start+3
            this_codon = str(codonseq[start:end])
            if len(this_codon) == 3:
                codon_lst.append(this_codon)
            else:
                codon_lst.append(str(this_codon.ungap()))
        elif str(codonseq[int(k):int(k)+3]) == "---":
            codon_lst.append("---")
        else:
          # this may be problematic, as normally no codon shoud
          # fall into this condition
            codon_lst.append(codonseq[int(k):int(k)+3])
    return codon_lst


def count_sites(codon_lst, k=1, codon_table=default_codon_table):
    S_site = 0.0  # synonymous sites
    N_site = 0.0  # non-synonymous sites
    purine = ('A', 'G')
    pyrimidine = ('T', 'C')
    base_tuple = ('A', 'T', 'C', 'G')
    for codon in codon_lst:
        neighbor_codon = {'transition': [], 'transversion': []}
        # classify neighbor codons
        codon = codon.replace('U', 'T')
        if codon == '---':
            continue
        for n, i in enumerate(codon):
            for j in base_tuple:
                if i == j:
                    pass
                elif i in purine and j in purine:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transition'].append(this_codon)
                elif i in pyrimidine and j in pyrimidine:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transition'].append(this_codon)
                else:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transversion'].append(this_codon)

        if codon in codon_table.stop_codons:
            print("STOP DETECTED")
            continue

        aa = codon_table.forward_table[codon]
        this_codon_N_site = this_codon_S_site = 0
        for neighbor in neighbor_codon['transition']:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += 1
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += 1
            else:
                this_codon_N_site += 1
        for neighbor in neighbor_codon['transversion']:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += k
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += k
            else:
                this_codon_N_site += k
        norm_const = (this_codon_N_site + this_codon_S_site)/3
        S_site += float(this_codon_S_site) / float(norm_const)
        N_site += float(this_codon_N_site) / float(norm_const)
    return (S_site, N_site)

def _count_diff_NG86(codon1, codon2, codon_table=default_codon_table):
    """Count differences between two codons (three-letter string; PRIVATE).

      The function will take multiple pathways from codon1 to codon2
      into account.
      """

    if not isinstance(codon1, str) or not isinstance(codon2, str):
        raise TypeError('_count_diff_NG86 accepts string object to represent codon ({0}, {1} detected)'.format(type(codon1),
                        type(codon2)))
    if len(codon1) != 3 or len(codon2) != 3:
        raise RuntimeError('codon should be three letter string ({0}, {1} detected)'.format(len(codon1),
                           len(codon2)))
    SN = [0, 0]  # synonymous and nonsynonymous counts
    if codon1 == '---' or codon2 == '---':
        return SN
    base_tuple = ('A', 'C', 'G', 'T')
    if not all(i in base_tuple for i in codon1):
        raise RuntimeError('Unrecognized character detected in codon1 {0} (Codons consist of A, T, C or G)'.format(codon1))
    if not all(i in base_tuple for i in codon2):
        raise RuntimeError('Unrecognized character detected in codon2 {0} (Codons consist of A, T, C or G)'.format(codon2))
    if codon1 == codon2:
        return SN
    else:
        diff_pos = []
        for (i, k) in enumerate(zip(codon1, codon2)):
            if k[0] != k[1]:
                diff_pos.append(i)

        def compare_codon(
            codon1,
            codon2,
            codon_table=default_codon_table,
            weight=1,
            ):
            """Compare two codon accounting for different pathways."""

            sd = nd = 0
            if len(set(map(codon_table.forward_table.get, [codon1,
                   codon2]))) == 1:
                # If both codons end up with the same aa
#                 for x in map(codon_table.forward_table.get, [codon1,
#                    codon2]):
#                     print(x)
#                 print(codon1, codon2)
                sd += weight
            else:
                nd += weight

            return (sd, nd)

        if len(diff_pos) == 1:
            SN = [i + j for (i, j) in zip(SN, compare_codon(codon1,
                  codon2, codon_table=codon_table))]
        elif len(diff_pos) == 2:
            codon2_aa = codon_table.forward_table[codon2]
            for i in diff_pos:
                temp_codon = codon1[:i] + codon2[i] + codon1[i + 1:]
                SN = [i + j for (i, j) in zip(SN, compare_codon(codon1,
                      temp_codon, codon_table=codon_table, weight=0.5))]
                SN = [i + j for (i, j) in zip(SN,
                      compare_codon(temp_codon, codon2,
                      codon_table=codon_table, weight=0.5))]
        elif len(diff_pos) == 3:
            codon2_aa = codon_table.forward_table[codon2]
            paths = list(permutations([0, 1, 2], 3))
            tmp_codon = []
            for p in paths:
                tmp1 = codon1[:p[0]] + codon2[p[0]] + codon1[p[0] + 1:]
                tmp2 = tmp1[:p[1]] + codon2[p[1]] + tmp1[p[1] + 1:]
                tmp_codon.append((tmp1, tmp2))
                SN = [i + j for (i, j) in zip(SN, compare_codon(codon1,
                      tmp1, codon_table, weight=0.5 / 3))]
                SN = [i + j for (i, j) in zip(SN, compare_codon(tmp1,
                      tmp2, codon_table, weight=0.5 / 3))]
                SN = [i + j for (i, j) in zip(SN, compare_codon(tmp2,
                      codon2, codon_table, weight=0.5 / 3))]
    return SN