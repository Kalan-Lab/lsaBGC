"""
Author: Rauf Salamzade
Kalan Lab - UW Madison - MMI
01/08/2021
"""

import os
import sys
from Bio import SeqIO
from collections import defaultdict
import argparse
import random
from ete3 import Tree

# read in list of colors
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
colors_file = dir_path + 'colors_200.txt'
COLORS = []
with open(colors_file) as ocf:
    COLORS = [x.strip() for x in ocf.readlines()]

def assignColorsToCOGs(cogs):
    """
    :param cogs: set of CoGs.
    :return: dictionary mapping each CoG to a hex color value.
    """
    random.shuffle(COLORS)
    cog_to_color = {}
    for i, c in enumerate(set(cogs)):
        cog_to_color[c] = COLORS[i]
    return(cog_to_color)

def parseGenbanks(gbk):
    """
    :param gbk: AntiSMASH Genbank file
    :return: list of lists (genes) where each sub-list corresponding to a gene contains:
        1. locus tag
        2. start (smaller) coordinate
        3. end (larger) coordinate
        4. gene direction
        5. protein sequence
    each gene is presented in the same order it appears in the AntiSMASH genbank file
    """
    genes = []
    with open(gbk) as ogbk:
        for rec in SeqIO.parse(ogbk, 'genbank'):
            for feature in rec.features:
                if feature.type == "CDS":
                    lt = feature.qualifiers.get('locus_tag')[0]
                    prot_seq = feature.qualifiers.get('translation')[0]
                    start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    direction = str(feature.location).split('(')[1].split(')')[0]
                    genes.append([lt, start, end, direction, prot_seq])
    return(genes)

def parseOrthoFinderMatrix(orthofinder_matrix, all_gene_lts):
    """
    :param orthofinder_matrix: OrthoFinderV2 matrix Orthogroups.csv file, should also include singleton orthologroups
    :param all_gene_lts: Set of all the relevant gene locus tag identifiers found in BGC Genbanks
    :return: dictionary mapping gene locus tags to CoGs
    """
    gene_to_cog = {}
    with open(orthofinder_matrix) as ofm:
        for i, line in enumerate(ofm):
            if i == 0: continue
            line = line.strip()
            ls = line.split('\t')
            cog = ls[0]
            for sgs in ls[1:]:
                for g in sgs.split(', '):
                    if g in all_gene_lts:
                        gene_to_cog[g] = cog
    return(gene_to_cog)

def bgsee(bgc_specs_file, orthofinder_matrix, dataset_label, output_prefix, core_phylogeny, species_phylogeny):

    ### vet input files quickly
    sys.stderr.write('Checking if input files exist and are in valid formatting...\n')
    try:
        if species_phylogeny:
            assert(os.path.isfile(species_phylogeny))
        assert(os.path.isfile(orthofinder_matrix))
        assert(os.path.isfile(bgc_specs_file))
        with open(bgc_specs_file) as obsf:
            for line in obsf:
                line = line.strip()
                gbk, sample = line.split('\t')
                assert(os.path.isfile(gbk))
                with open(gbk) as ogbk:
                    for rec in SeqIO.parse(ogbk, 'genbank'):
                        pass
    except:
        sys.stderr.write('Error validating input files. Exiting now...\n')
        sys.exit(1)

    sys.stderr.write('Input files check out, at least for now...\n')

    sample_index = defaultdict(int)
    edited_sample_names = []
    antismash_genbanks = []
    with open(bgc_specs_file) as obsf:
        for line in obsf:
            line = line.strip()
            gbk, sample = line.split('\t')
            try:
                assert (os.path.isfile(gbk))
                antismash_genbanks.append(gbk)
                if sample_index[sample] == 0:
                    edited_sample_names.append(sample)
                else:
                    edited_sample_names.append(sample + '_' + str(sample_index[sample]+1))
                sample_index[sample] += 1
            except:
                pass

    if species_phylogeny:
        t = Tree(species_phylogeny)
        for node in t.traverse('postorder'):
            if node.name in edited_sample_names and sample_index[node.name] > 0:
                for i in range(2, sample_index[node.name]+1):
                    node.add_sister(name=(node.name + '_' + str(i)))
                    sister_node = t.search_nodes(name=(node.name + '_' + str(i)))[0]
                    sister_node.dist = node.dist
        t.write(format=1, outfile=output_prefix + '.expanded.tre')

    # parse and process genbanks
    bgc_genes = {}
    all_gene_lts = set([])
    for i, gbk in enumerate(antismash_genbanks):
        bgc_name = edited_sample_names[i]
        genes = parseGenbanks(gbk)
        bgc_genes[bgc_name] = genes
        all_gene_lts = all_gene_lts.union(set([x[0] for x in genes]))

    # parse orthofinder matrix
    gene_to_cog = parseOrthoFinderMatrix(orthofinder_matrix, all_gene_lts)

    # assign each relevant CoG to a hexagonal color value
    cog_to_color = assignColorsToCOGs(gene_to_cog.values())

    # write header for iTol track file
    outf = open(output_prefix + '.iTol.txt', 'w')
    outf.write('DATASET_DOMAINS\n')
    outf.write('SEPARATOR TAB\n')
    outf.write('DATASET_LABEL\t%s\n' % dataset_label)
    outf.write('COLOR\t#000000\n')
    outf.write('#BORDER_WIDTH\t1\n')
    outf.write('BORDER_COLOR\t#000000\n')
    outf.write('SHOW_DOMAIN_LABELS\t0\n')
    outf.write('DATA\n')

    # write the rest of the iTol track file for illustrating genes across BGC instances
    ref_cog_directions = {}
    for i, s in enumerate(bgc_genes):
        samp_genes = bgc_genes[s]
        last_gene_end = samp_genes[-1][2]
        end = int(round(float(last_gene_end)/100.0)) + 1
        printlist = [s, str(end)]
        cog_directions = {}
        for ginfo in samp_genes:
            cog = 'singleton'
            if ginfo[0] in gene_to_cog:
                cog = gene_to_cog[ginfo[0]]
            shape = 'None'
            if ginfo[-2] == '+': shape = 'TR'
            elif ginfo[-2] == '-': shape = 'TL'
            gstart = int(round(float(ginfo[1])/100.0))+1
            gend = int(round(float(ginfo[2])/100.0))+1
            cog_color = "#dbdbdb"
            if cog in cog_to_color:
                cog_color = cog_to_color[cog]
            gene_string = '|'.join([str(x) for x in [shape, gstart, gend, cog_color, cog]])
            printlist.append(gene_string)
            if cog != 'singleton':
                cog_directions[cog] = ginfo[-2]
        if i == 0:
            ref_cog_directions = cog_directions
            outf.write('\t'.join(printlist) + '\n')
        else:
            flip_support = 0
            keep_support = 0
            for c in ref_cog_directions:
                if not c in cog_directions: continue
                if cog_directions[c] == ref_cog_directions[c]: keep_support += 1
                else: flip_support += 1

            # flip the genbank visual if necessary, first BGC processed is used as reference guide
            if flip_support > keep_support:
                flip_printlist = printlist[:2]
                bgc_end = int(printlist[1])
                for gene_string in printlist[2:]:
                    gene_info = gene_string.split('|')
                    new_shape = None
                    if gene_info[0] == 'TR': new_shape = 'TL'
                    elif gene_info[0] == 'TL': new_shape = 'TR'
                    new_gstart = int(bgc_end) - int(gene_info[2])
                    new_gend = int(bgc_end) - int(gene_info[1])
                    new_gene_info = '|'.join([new_shape, str(new_gstart), str(new_gend)] + gene_info[-2:])
                    flip_printlist.append(new_gene_info)
                outf.write('\t'.join(flip_printlist) + '\n')
            else:
                outf.write('\t'.join(printlist) + '\n')
    outf.close()

    if core_phylogeny:
        tmp_dir = output_prefix + '_tmp/'
        os.system('mkdir %s' % tmp_dir)

        # gather list of core CoGs
        all_bgcs = set([])
        bgc_cog_counts = defaultdict(lambda: defaultdict(list))
        bgc_sccs = defaultdict(lambda: "")
        for i, b in enumerate(bgc_genes):
            all_bgcs.add(b)
            samp_genes = bgc_genes[b]
            for ginfo in samp_genes:
                if not ginfo[0] in gene_to_cog: continue
                cog = gene_to_cog[ginfo[0]]
                prot_seq = ginfo[-1]
                bgc_cog_counts[cog][b].append(prot_seq)
        for c in bgc_cog_counts:
            flag_bgc_scc = True
            seqs = []
            for b in all_bgcs:
                if len(bgc_cog_counts[c][b]) != 1:
                    flag_bgc_scc = False
                    break
                else:
                    seqs.append([b, bgc_cog_counts[c][b][0]])
            if flag_bgc_scc:
                c = '_'.join(c.split())
                # write protein fasta for cog
                prot_faa = tmp_dir + c + '.faa'
                prot_handle = open(prot_faa, 'w')
                for g in seqs:
                    prot_handle.write('>' + g[0] + '\n' + g[1] + '\n')
                prot_handle.close()

                # align protein alignment using muscle
                prot_msa = tmp_dir + c + '.msa.faa'
                os.system('muscle < %s > %s' % (prot_faa, prot_msa))

                # concatenate gene alignments
                with open(prot_msa) as opm:
                    for rec in SeqIO.parse(opm, 'fasta'):
                        bgc_sccs['>' + rec.id] += str(rec.seq).upper()

        os.system('rm -rf %s' % tmp_dir)

        fasta_data = []
        for b in bgc_sccs:
            fasta_data.append([b] + list(bgc_sccs[b]))

        fasta_data_tr = []
        for i, ls in enumerate(zip(*fasta_data)):
            if i == 0:
                fasta_data_tr.append(ls)
            else:
                n_count = len([x for x in ls if x == '-'])
                if (float(n_count) / len(ls)) < 0.1:
                    fasta_data_tr.append(list(ls))

        scc_faa = output_prefix + '.msa.faa'
        scc_tre = output_prefix + '.msa.tre'
        scc_handle = open(scc_faa, 'w')

        for rec in zip(*fasta_data_tr):
            scc_handle.write(rec[0] + '\n' + ''.join(rec[1:]) + '\n')
        scc_handle.close()

        # use FastTree2 to construct phylogeny
        os.system('fasttree %s > %s' % (scc_faa, scc_tre))

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Program to visualize BGCs using the iTol Infrastructure.
    """)

    parser.add_argument('-b', '--bgc_specs_file', help='BGC specifications file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-l', '--dataset_label', help="Dataset label for iTol track", required=False, default='BGSee')
    parser.add_argument('-o', '--output_prefix', help="Prefix for output files.", required=False, default='iTol_BGSee.iTol.txt')
    parser.add_argument('-p', '--create_core_gcf_phylogeny', action='store_true', help="Create phylogeny from core COGs.", required=False, default=False)
    parser.add_argument('-s', '--species_phylogeny', help="The species phylogeny in Newick format.", required=False, default=None)
    args = parser.parse_args()

    bgsee(args.bgc_specs_file, args.orthofinder_matrix, args.dataset_label, args.output_prefix, args.create_core_gcf_phylogeny, args.species_phylogeny)