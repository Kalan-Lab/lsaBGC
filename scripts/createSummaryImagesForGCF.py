#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict
from operator import itemgetter
from ete3 import Tree

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
    Program: createSummaryImagesForGCF.py
    Author: Rauf Salamzade
    Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_input_listing', help="Path to GCF listing file.", required=True)
    parser.add_argument('-p', '--gcf_popgene_result', help="Path to Population Genetics result directory.", required=True)
    parser.add_argument('-t', '--species_tree', help="Species phylogeny.", required=True)
    parser.add_argument('-o', '--output_dir', help="Path to result directory.", required=True)
    args = parser.parse_args()
    return args

def is_number(x):
    try:
        x = float(x)
        return True
    except:
        return False

def main():
    """
    Function to create a track file for visualizing BGC gene architecture across a phylogeny in the interactive tree
    of life (iTol)
    """

    myargs = create_parser()
    gcf_input_listing_file = myargs.gcf_input_listing
    gcf_popgene_result_dir = myargs.gcf_popgene_result
    species_tree_file = myargs.species_tree
    output_dir = myargs.output_dir

    try:
        assert(os.path.isfile(gcf_input_listing_file) and os.path.isdir(gcf_popgene_result_dir) and os.path.isfile(species_tree_file))
        output_dir = os.path.abspath(output_dir) + '/'
    except:
        raise RuntimeError()

    species_of_interest = set([])
    t = Tree(species_tree_file)
    for leaf in t:
        species_of_interest.add(str(leaf).strip('\n').lstrip('-'))

    popgen_cod_dir = gcf_popgene_result_dir + 'Codon_PopGen_Analyses/'

    consensus_similarity = defaultdict(lambda: defaultdict(list))
    for f in os.listdir(popgen_cod_dir):
        if not f.endswith('_sim_to_consensus.txt'): continue
        with open(popgen_cod_dir + f) as opcdf:
            for line in opcdf:
                line = line.strip()
                ls = line.split('\t')
                species = ls[1].split('_RS_')[0].split('_GB_')[0]
                consensus_similarity[ls[0]][species].append(float(ls[2]))

    popgen_stats_file = output_dir + 'Population_Genetics_Tracks.txt'
    consens_heatmap_file = output_dir + 'Consensus_Similarity_Heatmap.txt'
    species_count_file = output_dir + 'Species_GCF_Carriage_Track.txt'
    annotation_file = output_dir + 'Annotation_Track.txt'

    psf_handle = open(popgen_stats_file, 'w')
    chf_handle = open(consens_heatmap_file, 'w')
    scf_handle = open(species_count_file, 'w')
    af_handle = open(annotation_file, 'w')

    psf_handle.write('\t'.join(['hg', 'hg_order', 'hg_start', 'hg_end', 'samples_with_hg', 'tajimas_d', 'beta_rd']) + '\n')
    chf_handle.write('\t'.join(['hg', 'hg_order', 'hg_start', 'hg_end', 'hg_direction', 'label', 'median_consensus_difference']) + '\n')
    scf_handle.write('\t'.join(['label', 'isolates_with_gcf']) + '\n')
    af_handle.write('\t'.join(['hg', 'hg_order', 'hg_start', 'hg_end', 'manual_annotation', '#core', 'annotation', 'domains']) + '\n')

    data = []
    with open(gcf_popgene_result_dir + 'Ortholog_Group_Information.txt') as ogprf:
        for i, line in enumerate(ogprf):
            if i == 0: continue
            line = line.strip('\n')
            ls = line.split('\t')
            hg, annot, hg_order, hg_dir = ls[2:6]
            hg_med_len, core = ls[7:9]
            num_samples = ls[10]
            tajimas_d = ls[13]
            median_beta_rd = ls[16]
            domains = ls[-1]
            if int(num_samples) >= 3:
                data.append([hg, core, annot, domains, int(hg_order), hg_dir, float(hg_med_len), num_samples, tajimas_d, median_beta_rd])

    start = 1
    for i, hgi in enumerate(sorted(data, key=itemgetter(4))):
        end = start + hgi[6]
        psf_handle.write('\t'.join([str(x) for x in ([hgi[0], hgi[4], start, end] + hgi[-3:])]) + '\n')
        af_handle.write('\t'.join([hgi[0], hgi[4], start, end, 'XXXX', '#' + hgi[1:4]]) + '\n')
        for species in species_of_interest:
            hg_spec_consensus_differences = consensus_similarity[hgi[0]][species]
            if len(hg_spec_consensus_differences) > 0:
                median_consensus_diff = statistics.median(hg_spec_consensus_differences)
                chf_handle.write('\t'.join([str(x) for x in ([hgi[0], hgi[4], start, end, hgi[5], species, median_consensus_diff])]) + '\n')
        start += hgi[6] + 1

    species_iso_with_gcf = defaultdict(int)
    with open(gcf_input_listing_file) as ogilf:
        for line in ogilf:
            line = line.strip()
            ls = line.split('\t')
            species = ls[0].split('_RS_')[0].split('_GB_')[0]
            species_iso_with_gcf[species] += 1

    for species in species_of_interest:
        scf_handle.write('\t'.join([species, str(species_iso_with_gcf[species])]) + '\n')



main()