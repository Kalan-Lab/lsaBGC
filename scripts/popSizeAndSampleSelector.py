#!/usr/bin/env python

### Program: lsaBGC-AutoAnalyze.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2021, Kalan-Lab
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#	 list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#	 this list of conditions and the following disclaimer in the documentation
#	 and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#	 contributors may be used to endorse or promote products derived from
#	 this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import sys
from time import sleep
from operator import itemgetter
import argparse
from collections import defaultdict
from ete3 import Tree
from Bio import SeqIO
from lsaBGC.classes.Pan import Pan
from lsaBGC import util
import itertools

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_NJTREECONSTRUCTION = lsaBGC_main_directory + '/lsaBGC/Rscripts/createNJTree.R'
RSCRIPT_FOR_DEFINECLADES_FROM_PHYLO = lsaBGC_main_directory + '/lsaBGC/Rscripts/defineCladesFromPhylo.R'

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: popSizeAndSampleSelector.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-l', '--input_listing', help="Path to tab delimited file listing: (1) sample name (2) path to Prokka Genbank and (3) path to Prokka predicted proteome. This file is produced by lsaBGC-Process.py.", required=True, default=None)
    parser.add_argument('-o', '--output_directory', help="Parent output/workspace directory.", required=True)
    parser.add_argument('-k', '--sample_set', help="Sample set to keep in analysis. Should be file with one sample id per line.", required=False)
    parser.add_argument('-s', '--lineage_phylogeny', help="Path to species phylogeny. If not provided a MASH based neighborjoining tree will be constructed and used.", default=None, required=False)
    parser.add_argument('-lps', '--lower_num_populations', type=int, help='If population analysis specified, what is the lower number of populations to fit to . Use the script determinePopulationK.py to see how populations will look with k set to different values.', required=False, default=2)
    parser.add_argument('-ups', '--upper_num_populations', type=int, help='If population analysis specified, what is the number of populations to . Use the script determinePopulationK.py to see how populations will look with k set to different values.', required=False, default=20)
    parser.add_argument('-i', '--identity_cutoff', type=float, help='Identity to collapse samples at.', required=False, default=0.98)
    parser.add_argument('-c', '--cores', type=int, help="Total number of cores to use.", required=False, default=1)

    args = parser.parse_args()
    return args

def determineN50(assembly_fasta):
    """
    Code borrowed from:
    https://onestopdataanalysis.com/n50-genome/
    """
    scaffold_lengths = []
    with open(assembly_fasta) as oaf:
        for rec in SeqIO.parse(oaf, 'fasta'):
            scaffold_lengths.append(len(str(rec.seq)))
    tmp = []
    for tmp_number in set(scaffold_lengths):
        tmp += [tmp_number] * scaffold_lengths.count(tmp_number) * tmp_number
    tmp.sort()

    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
    return(median)

def singleLinkageClustering(similar_pairs):
    """
    Solution for single-linkage clustering taken from mimomu's repsonse in the stackoverflow page:
    https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements?lq=1
    """
    L = similar_pairs
    LL = set(itertools.chain.from_iterable(L))
    for each in LL:
        components = [x for x in L if each in x]
        for i in components:
            L.remove(i)
        L += [list(set(itertools.chain.from_iterable(components)))]
    return(L)

def lsaBGC_AutoAnalyze():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """
    myargs = create_parser()

    outdir = os.path.abspath(myargs.output_directory) + '/'
    input_listing_file = os.path.abspath(myargs.input_listing)

    try:
        assert (os.path.isfile(input_listing_file))
    except:
        raise RuntimeError(
            'Input listings file does not exist. Exiting now ...')

    if os.path.isdir(outdir):
        sys.stderr.write("Output directory exists. Continuing in 5 seconds ...\n ")
        sleep(5)
    else:
        os.system('mkdir %s' % outdir)

    """
    PARSE OPTIONAL INPUTS
    """

    sample_set_file = myargs.sample_set
    lineage_phylogeny_file = myargs.lineage_phylogeny
    lower_num_populations = myargs.lower_num_populations
    upper_num_populations = myargs.upper_num_populations
    identity_cutoff = myargs.identity_cutoff
    cores = myargs.cores

    """
    START WORKFLOW
    """

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for easier determination of results' provenance in the future.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [input_listing_file, outdir,  lineage_phylogeny_file, sample_set_file, identity_cutoff,
                        lower_num_populations, upper_num_populations, cores]
    parameter_names = ["Listing File of Prokka Annotation Files for Initial Set of Samples",
                       "Output Directory", "Phylogeny File in Newick Format", "Sample Retention Set", "Identity Cutoff",
                       "Lower Limit for Number of Populations", "Upper Limit for Number of Populations", "Cores"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Step 0: (Optional) Parse sample set retention specifications file, if provided by the user.
    sample_retention_set = util.getSampleRetentionSet(sample_set_file)

    Pan_Object = Pan(input_listing_file, logObject=logObject)
    logObject.info("Converting Genbanks from Expansion listing file into FASTA per sample.")
    gw_fasta_dir = outdir + 'Sample_Assemblies_in_FASTA/'
    if not os.path.isdir(gw_fasta_dir): os.system('mkdir %s' % gw_fasta_dir)
    gw_fasta_listing_file = outdir + 'Genome_FASTA_Listings.txt'
    Pan_Object.convertGenbanksIntoFastas(gw_fasta_dir, gw_fasta_listing_file)
    logObject.info("Successfully performed conversions.")

    logObject.info("Running MASH Analysis Between Genomes.")
    gw_pairwise_differences = util.calculateMashPairwiseDifferences(gw_fasta_listing_file, outdir, 'genome_wide', 10000,
                                                                    cores, logObject, prune_set=sample_retention_set)

    logObject.info("Ran MASH Analysis Between Genomes.")
    mash_matrix_file = outdir + 'MASH_Distance_Matrix.txt'
    mash_matrix_handle = open(mash_matrix_file, 'w')
    mash_matrix_handle.write('Sample/Sample\t' + '\t'.join([s for s in sorted(gw_pairwise_differences)]) + '\n')

    all_samples = set([])
    similar_pairs = []
    for s1 in sorted(gw_pairwise_differences):
        printlist = [s1]
        all_samples.add(s1)
        for s2 in sorted(gw_pairwise_differences):
            if (1.0 - gw_pairwise_differences[s1][s2]) >= identity_cutoff:
                pair_id = sorted([s1, s2])
                similar_pairs.append(pair_id)
            printlist.append(str(gw_pairwise_differences[s1][s2]))
        mash_matrix_handle.write('\t'.join(printlist) + '\n')
    mash_matrix_handle.close()

    if not lineage_phylogeny_file:
        # Run MASH Analysis Between Genomic Assemblies
        logObject.info("Using MASH estimated distances between genomes to infer neighbor-joining tree.")
        mash_nj_tree = outdir + 'MASH_NeighborJoining_Tree.nwk'

        # create neighbor-joining tree
        cmd = ['Rscript', RSCRIPT_FOR_NJTREECONSTRUCTION, mash_matrix_file, mash_nj_tree]
        try:
            util.run_cmd(cmd, logObject)
            assert (util.is_newick(mash_nj_tree))
            lineage_phylogeny_file = mash_nj_tree
        except Exception as e:
            logObject.error(
                "Had issues with creating neighbor joining tree and defining populations using treestructure.")
            raise RuntimeError(
                "Had issues with creating neighbor joining tree and defining populations using treestructure.")

    if sample_retention_set == None:
        sample_retention_set = all_samples

    redundant_samples = set([])
    for cluster in singleLinkageClustering(similar_pairs):
        max_n50_rep = [0.0, None]
        for sample in cluster:
            sample_assembly_fasta = gw_fasta_dir + sample + '.fasta'
            sample_assembly_n50 = determineN50(sample_assembly_fasta)
            if sample_assembly_n50 > max_n50_rep[0]:
                max_n50_rep = [sample_assembly_n50, sample]
        assert(max_n50_rep[1] != None)
        for sample in cluster:
            if sample != max_n50_rep[1]:
                redundant_samples.add(sample)

    sample_retention_set = sample_retention_set.difference(redundant_samples)

    # Pruning lineage phylogeny provided
    logObject.info("Pruning lineage phylogeny to retain only samples of interest.")
    update_lineage_phylogeny_file = outdir + 'Lineage_Phylogeny.Pruned.nwk'
    update_sample_list_file = outdir + 'Sample_Listing_Keep.txt'

    update_sample_list_handle = open(update_sample_list_file, 'w')
    for sample in sample_retention_set:
        update_sample_list_handle.write(sample + '\n')
    update_sample_list_handle.close()

    t = Tree(lineage_phylogeny_file)
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    t.prune(sample_retention_set)
    t.write(format=5, outfile=update_lineage_phylogeny_file)
    lineage_phylogeny_file = update_lineage_phylogeny_file
    logObject.info("Successfully refined lineage phylogeny for sample set of interest.")

    result_dir = outdir + 'Population_Mapping_Results/'
    if not os.path.isdir(result_dir): os.system('mkdir %s' % result_dir)
    for ps in range(lower_num_populations, upper_num_populations+1):
        ps_dir = result_dir + 'PopSize_' + str(ps) + '/'
        if not os.path.isdir(ps_dir): os.system('mkdir %s' % ps_dir)
        population_listing_file = ps_dir + 'Populations_Defined.txt'
        populations_on_nj_tree_pdf = ps_dir + 'Populations_on_Lineage_Tree.pdf'

        # Attempt population fitting
        cmd = ['Rscript', RSCRIPT_FOR_DEFINECLADES_FROM_PHYLO, lineage_phylogeny_file, str(ps),
               population_listing_file, populations_on_nj_tree_pdf]
        try:
            util.run_cmd(cmd, logObject)
        except Exception as e:
            logObject.error(
                "Had issues with creating neighbor joining tree and defining populations using cutree.")
            raise RuntimeError(
                "Had issues with creating neighbor joining tree and defining populations using cutree.")

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_AutoAnalyze()
