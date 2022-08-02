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
import traceback
from collections import defaultdict
from ete3 import Tree
from lsaBGC.classes.Pan import Pan
from lsaBGC import util
from Bio import SeqIO

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_NJTREECONSTRUCTION = lsaBGC_main_directory + '/lsaBGC/Rscripts/createNJTree.R'

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-AutoAnalyze.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

    Program to compute pairwise genome-wide similarities needed for the Beta-RD statistic by lsaBGC-Divergence.py and 
    lsaBGC-PopGene.py and also uses them to perform a quick dereplication (based on a dynamic algorithm using N50 to 
    keep track of best representatives). For more reliable dereplication we recommend using dRep by Dr. Matt Olm. 
    lsaBGC-PopGene.py and lsaBGC-Divergence.py each feature the "--sample_set" option to take in a list of sample names
    to use for downstream analyses.
    
    This program can also optionally construct a species tree (via neighbor joining) and perform 
    population classification of genomes/samples. 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-l', '--sample_annotation_listing',
                        help="Path to tab delimited file listing: (1) sample name (2) path to Prokka Genbank and (3) path to Prokka predicted proteome. This file is produced by lsaBGC-Process.py.",
                        required=True)
    parser.add_argument('-o', '--output_directory', help="Parent output/workspace directory.", required=True)
    parser.add_argument('-t', '--create_species_tree', action='store_true',
                        help="Use neighbor-joining to create a species tree off pairwise distances.",
                        default=False, required=False)
    parser.add_argument('-i', '--identity_cutoff', type=float, help='Cutoff for ANI/AAI for grouping samples together if population inference requested.', required=False, default=0.99)
    parser.add_argument('-s', '--shared_cutoff', type=float, help='Cutoff for shared genomic content for grouping samples together if population inference requested.', required=False, default=0.8)
    parser.add_argument('-m', '--method', help='Which method to use for estimating ANI or AAI? Options: "CompareM", "FastANI", "MASH". Default is FastANI.', required=False, default="FastANI")
    parser.add_argument('-c', '--cpus', type=int, help="Total number of cpus to use.", required=False, default=1)

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

def main():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """
    myargs = create_parser()

    sample_annotation_listing_file = os.path.abspath(myargs.sample_annotation_listing)
    outdir = os.path.abspath(myargs.output_directory) + '/'

    try:
        assert (os.path.isfile(sample_annotation_listing_file))
    except:
        raise RuntimeError('The input file for Sample Annotation Listings does not exist. Exiting now ...')

    if os.path.isdir(outdir):
        sys.stderr.write("Output directory exists. Continuing in 5 seconds ...\n ")
        sleep(5)
    else:
        os.system('mkdir %s' % outdir)

    """
    PARSE OPTIONAL INPUTS
    """

    create_species_tree = myargs.create_species_tree
    identity_cutoff = myargs.identity_cutoff
    shared_cutoff = myargs.shared_cutoff
    method = myargs.method
    cpus = myargs.cpus

    try:
        assert(method.upper() in set(['MASH', 'FASTANI', 'COMPAREM']))
    except:
        sys.stderr.write("Selected method for calculating ANI/AAI not one of the accepted options. Exiting now ...")

    """
    START WORKFLOW
    """

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for easier determination of results' provenance in the future.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [sample_annotation_listing_file, outdir, create_species_tree, identity_cutoff,
                        shared_cutoff, method, cpus]
    parameter_names = ["Input Listing File", "Output Directory", "Construct Species Tree?", "Identity Cutoff",
                       "Shared Cutoff", "Method for Inferring Genome-Wide ANI/AAI", "cpus"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    gw_fasta_listing_file = outdir + 'Genome_FASTA_Listings.txt'
    Pan_Object = Pan(sample_annotation_listing_file, logObject=logObject)
    logObject.info("Converting Genbanks from listing file into FASTA per sample.")
    gw_fasta_dir = outdir + 'Sample_Assemblies_in_FASTA/'
    util.setupReadyDirectory([gw_fasta_dir])
    Pan_Object.convertGenbanksIntoFastas(gw_fasta_dir, gw_fasta_listing_file)
    logObject.info("Successfully performed conversions.")

    logObject.info("Running %s Analysis to Assess Similarity Between Genomes." % method)
    gw_pairwise_identities = None
    gw_pairwise_shared = None
    if method.upper() == 'MASH':
        gw_pairwise_identities = util.calculateMashPairwiseDifferences(gw_fasta_listing_file, outdir, 'genome_wide',
                                                                       10000, cpus, logObject)
    elif method.upper() == 'FASTANI':
        fastani_results_file = outdir + 'fastANI_Results.txt'
        gw_pairwise_identities, gw_pairwise_shared = util.runFastANI(gw_fasta_listing_file, outdir,
                                                                     fastani_results_file, cpus, logObject)
    elif method.upper() == 'COMPAREM':
        comparem_results_dir = outdir + 'CompareM/'
        util.setupReadyDirectory([comparem_results_dir])
        gw_pairwise_identities, gw_pairwise_shared = util.runCompareM(gw_fasta_listing_file,
                                                                             comparem_results_dir, cpus, logObject)

    sample_assembly_n50s = {}
    with open(gw_fasta_listing_file) as ogf:
        for line in ogf:
            line = line.strip()
            s, sample_assembly_fasta = line.split('\t')
            sample_assembly_n50 = determineN50(sample_assembly_fasta)
            sample_assembly_n50s[s] = sample_assembly_n50

    logObject.info("Ran MASH/FastANI/CompareM Analysis Between Genomes.")
    mash_matrix_file = outdir + 'Distance_Matrix.txt'
    mash_matrix_handle = open(mash_matrix_file, 'w')
    mash_matrix_handle.write('Sample/Sample\t' + '\t'.join([s for s in sorted(gw_pairwise_identities)]) + '\n')

    genome_wise_est_file = outdir + 'Genome_Wide_Estimates.txt'
    genome_wise_est_handle = open(genome_wise_est_file, 'w')

    similar_samples = []
    similar_sample_set = set([])
    all_samples = set([])
    redundant_samples = set([])
    poor_n50_samples = set([])
    for i, s1 in enumerate(sorted(gw_pairwise_identities)):
        s1_n50 = sample_assembly_n50s[s1]
        if s1_n50 < 10000:
            poor_n50_samples.add(s1)
        printlist = [s1]
        all_samples.add(s1)
        for j, s2 in enumerate(sorted(gw_pairwise_identities)):
            if gw_pairwise_identities[s1][s2] >= identity_cutoff:
                if s1 != s2 and i < j:
                    coverage_met_flag = False
                    if method.upper() == 'FASTANI' or method.upper() == 'COMPAREM':
                        if gw_pairwise_shared[s1][s2] >= shared_cutoff and gw_pairwise_shared[s2][s1] >= shared_cutoff:
                            coverage_met_flag = True
                    elif method.upper() == 'MASH':
                        coverage_met_flag = True
                    if coverage_met_flag:
                        s2_n50 = sample_assembly_n50s[s2]
                        if s1_n50 >= s2_n50:
                            redundant_samples.add(s2)
                        else:
                            redundant_samples.add(s1)
            if s1 != s2:
                genome_wise_est_handle.write(s1 + '\t' + s2 + '\t' + str(gw_pairwise_identities[s1][s2]) + '\n')
            printlist.append(str(1.0-gw_pairwise_identities[s1][s2]))
        mash_matrix_handle.write('\t'.join(printlist) + '\n')
    mash_matrix_handle.close()
    genome_wise_est_handle.close()

    sample_retention_set = all_samples.difference(redundant_samples).difference(poor_n50_samples)
    update_sample_list_file = outdir + 'Sample_Listing_Keep.txt'
    update_sample_list_handle = open(update_sample_list_file, 'w')
    for sample in sample_retention_set:
        update_sample_list_handle.write(sample + '\n')
    update_sample_list_handle.close()

    if create_species_tree:
        logObject.info("Building ANI/AAI based tree using neighbor-joining tree.")
        mash_nj_tree = outdir + 'Neighbor_Joining_Tree.nwk'

        # create neighbor-joining tree
        cmd = ['Rscript', RSCRIPT_FOR_NJTREECONSTRUCTION, mash_matrix_file, mash_nj_tree]
        try:
            util.run_cmd(cmd, logObject)
            assert (util.is_newick(mash_nj_tree))
            lineage_phylogeny_file = mash_nj_tree
        except Exception as e:
            logObject.error("Had issues with creating neighbor joining tree and defining populations using treestructure.")
            raise RuntimeError("Had issues with creating neighbor joining tree and defining populations using treestructure.")

        # Pruning lineage phylogeny provided
        logObject.info("Pruning neighbor-joining tree to retain only samples after dereplication.")
        update_lineage_phylogeny_file = outdir + 'Neighbor_Joining_Tree_Dereplicated.nwk'


        t = Tree(lineage_phylogeny_file)
        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)
        t.prune(sample_retention_set)
        t.write(format=5, outfile=update_lineage_phylogeny_file)
        logObject.info("Successfully refined lineage phylogeny for sample set of interest.")

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    main()
