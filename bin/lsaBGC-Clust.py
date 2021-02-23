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
import itertools
import statistics

def gatherAnnotation(gbk):
    product = "NA"
    with open(gbk) as ogbk:
        for rec in SeqIO.parse(ogbk, 'genbank'):
            for feature in rec.features:
                if feature.type == "region":
                    product = feature.qualifiers.get('product')[0]
    return product

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
    avg_nz_cog_copycounts = defaultdict(int)
    with open(orthofinder_matrix) as ofm:
        for i, line in enumerate(ofm):
            if i == 0: continue
            line = line.strip('\n')
            ls = line.split('\t')
            cog = ls[0]
            flag = False
            counts = []
            for sgs in ls[1:]:
                for g in sgs.split(', '):
                    if g in all_gene_lts:
                        gene_to_cog[g] = cog
                        flag = True
                if len(sgs.split(', ')) > 0:
                    counts.append(len(sgs.split(', ')))
            avg_nz_cog_copycounts[cog] = float(sum(counts))/float(len(counts))
    return(gene_to_cog, avg_nz_cog_copycounts)

def bgclust(bgc_specs_file, orthofinder_matrix, outdir, threads, mcl_inflation, run_inflation_tests):
    ### vet input files quickly
    sys.stderr.write('Checking if input files exist and are in valid formatting...\n')
    try:
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

    outdir = os.path.abspath(outdir) + '/'
    if not os.path.isdir(outdir):
        sys.stderr.write("Output directory does not exist. Creating it ...\n")
        os.system('mkdir %s' % outdir)

    sys.stderr.write('Reading in and processing GBKs ...\n')
    antismash_genbank_to_sample = {}
    with open(bgc_specs_file) as obsf:
        for line in obsf:
            line = line.strip()
            gbk, sample = line.split('\t')
            antismash_genbank_to_sample[gbk] = sample

    # parse and process genbanks
    bgc_genes = {}
    bgc_annot = {}
    all_gene_lts = set([])
    sample_bgc_indices = defaultdict(lambda: 1)
    bgc_gbks = {}
    name_to_sample = {}
    for gbk in antismash_genbank_to_sample:
        samp_name = antismash_genbank_to_sample[gbk]
        bgc_name = samp_name + '_' + str(sample_bgc_indices[samp_name])
        name_to_sample[bgc_name] = samp_name
        genes = parseGenbanks(gbk)
        annotation = gatherAnnotation(gbk)
        bgc_genes[bgc_name] = genes
        bgc_annot[bgc_name] = annotation
        bgc_gbks[bgc_name] = gbk
        all_gene_lts = all_gene_lts.union(set([x[0] for x in genes]))
        sample_bgc_indices[samp_name] += 1

    # parse orthofinder matrix
    gene_to_cog, avg_nz_cog_counts = parseOrthoFinderMatrix(orthofinder_matrix, all_gene_lts)

    # write the rest of the iTol track file for illustrating genes across BGC instances
    bgc_cogs = defaultdict(set)
    auto_singleton_bgcs = set([])
    for i, b in enumerate(bgc_genes):
        bgs = bgc_genes[b]
        cogs = set([])
        for ginfo in bgs:
            if ginfo[0] in gene_to_cog:
                cog = gene_to_cog[ginfo[0]]
                cogs.add(cog)
        if len(cogs) > 0:
            bgc_cogs[b] = cogs
        else:
            auto_singleton_bgcs.add(b)

    sys.stderr.write('Calculating overlap in Ortholog Groups between BGC GBKs ...\n')
    pair_relations_file = outdir + 'bgc_pair_relationships.txt'
    prf_handle = open(pair_relations_file, 'w')
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

    mci_file = outdir + 'bgc_pair_relationships.mci'
    tab_file = outdir + 'bgc_pair_relationships.tab'
    mcxload_cmd = 'mcxload -abc %s --stream-mirror -write-tab %s -o %s' % (pair_relations_file, tab_file, mci_file)
    sys.stderr.write('Converting format of pair relationship file via mxcload ...\n')
    os.system(mcxload_cmd)

    stats_file = outdir + 'GCF_details.txt'
    sf_handle = open(stats_file, 'w')
    sf_handle.write('\t'.join(['inflation parameter', 'GCF id', 'number of BGCs',
                               'samples with multiple BGCs in GCF', 'size of the SCC', 'mean number of OGs',
                               'stdev for number of OGs', 'min difference', 'max difference', 'annotations']) + '\n')

    if run_inflation_tests:
        mcl_inflation_params = [0.8, 1.4, 2, 4, 8]
        for i in mcl_inflation_params:
            mcl_out = outdir + 'mcl.' + str(i).replace('.', '_') + '.out'
            mcxdump_out = outdir + 'final_mcl.' + str(i).replace('.', '_') + '.out'
            mcl_cmd = 'mcl %s -I %f -o %s -te %d' % (mci_file, i, mcl_out, threads)
            mcxdump_cmd = 'mcxdump -icl %s -tabr %s -o %s' % (mcl_out, tab_file, mcxdump_out)

            sys.stderr.write('Running MCL with inflation parameter set to %f ...\n' % i)
            #os.system(mcl_cmd)

            sys.stderr.write('Dumping results in human-readable format ...\n')
            #os.system(mcxdump_cmd)

            clustered_bgcs = set([])
            with open(mcxdump_out) as omo:
                for j, gcf in enumerate(omo):
                    gcf = gcf.strip()
                    gcf_mems = gcf.split()
                    if len(gcf_mems) < 2: continue
                    diffs = set([])
                    samp_counts = defaultdict(int)
                    samp_ogs = defaultdict(set)
                    annots = set([])
                    for a, bgc1 in enumerate(gcf_mems):
                        samp_counts[bgc1.split('_')[0]] += 1
                        annots.add(bgc_annot[bgc1])
                        samp_ogs[bgc1.split('_')[0]] = samp_ogs[bgc1.split('_')[0]].union(bgc_cogs[bgc1])
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
                    gcf_stats = [i, 'gcf_' + str(j), len(gcf_mems), multi_same_sample, len(scc), statistics.mean(num_ogs), statistics.stdev(num_ogs), min(diffs), max(diffs), '; '.join(annots)]
                    sf_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
            singleton_bgcs = set([])
            for bgc in bgc_cogs:
                if not bgc in clustered_bgcs: singleton_bgcs.add(bgc)
            sf_handle.write('\t'.join([str(x) for x in ([i, 'singletons', len(singleton_bgcs)] + ['NA']*7)]) + '\n')
    else:
        mcl_out = outdir + 'mcl.' + str(i).replace('.', '_') + '.out'
        mcxdump_out = outdir + 'final_mcl.' + str(i).replace('.', '_') + '.out'
        mcl_cmd = 'mcl %s -I %f -o %s -te %d' % (mci_file, i, mcl_out, threads)
        mcxdump_cmd = 'mcxdump -icl %s -tabr %s -o %s' % (mcl_out, tab_file, mcxdump_out)

        sys.stderr.write('Running MCL with inflation parameter set to %f ...\n' % i)
        os.system(mcl_cmd)

        sys.stderr.write('Dumping results in human-readable format ...\n')
        os.system(mcxdump_cmd)

        gcf_listing_dir = outdir + 'GCF_listings/'
        if not os.path.isdir(gcf_listing_dir): os.system('mkdir %s' % gcf_listing_dir)
        with open(mcxdump_out) as omo:
            for j, gcf in enumerate(omo):
                gcf = gcf.strip()
                gcf_mems = gcf.split()
                if len(gcf_mems) < 2: continue
                outf_list = open(gcf_listing_dir + 'gcf_' + str(j+1) + '.txt', 'w')
                for bgc in gcf_mems:
                    sname = '_'.join(bgc.split('_')[:-1])
                    gbkpath = bgc_gbks[bgc]
                    outf_list.write('%s\t%s\n' % (gbkpath, sname))
                outf_list.close()

        clustered_bgcs = set([])
        with open(mcxdump_out) as omo:
            for j, gcf in enumerate(omo):
                gcf = gcf.strip()
                gcf_mems = gcf.split()
                if len(gcf_mems) < 2: continue
                diffs = set([])
                samp_counts = defaultdict(int)
                samp_ogs = defaultdict(set)
                annots = set([])
                for a, bgc1 in enumerate(gcf_mems):
                    samp_counts[bgc1.split('_')[0]] += 1
                    annots.add(bgc_annot[bgc1])
                    samp_ogs[bgc1.split('_')[0]] = samp_ogs[bgc1.split('_')[0]].union(bgc_cogs[bgc1])
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
                gcf_stats = [i, 'gcf_' + str(j), len(gcf_mems), multi_same_sample, len(scc), statistics.mean(num_ogs), statistics.stdev(num_ogs), min(diffs), max(diffs), '; '.join(annots)]
                sf_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
        singleton_bgcs = set([])
        for bgc in bgc_cogs:
            if not bgc in clustered_bgcs: singleton_bgcs.add(bgc)
        sf_handle.write('\t'.join([str(x) for x in ([i, 'singletons', len(singleton_bgcs)] + ['NA']*7)]) + '\n')
    sf_handle.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    """)

    parser.add_argument('-b', '--bgc_specs_file', help='BGC specifications file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-o', '--output_directory', help="Output directory.", required=True)
    parser.add_argument('-t', '--threads', type=int, help="Number of threads to use for MCL step.", required=False, default=1)
    parser.add_argument('-i', '--mcl_inflation', type=float, help="Inflation parameter to be used for MCL.", required=False, default=1.4)
    parser.add_argument('-r', '--run_inflation_tests', action='store_true', help="Run tests for selecting best inflation parameter for MCL analysis and exit.", default=False, required=False)
    args = parser.parse_args()

    bgclust(args.bgc_specs_file, args.orthofinder_matrix, args.output_directory, args.threads, args.mcl_inflation, args.run_inflation_tests)