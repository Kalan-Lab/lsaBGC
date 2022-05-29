#!/usr/bin/env python

### Program: lsaBGC-Ready.py
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
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
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
import argparse
from Bio import SeqIO
from collections import defaultdict
from time import sleep
from lsaBGC import util
import subprocess
import traceback
import multiprocessing


def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Ready.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Program to convert existing antiSMASH (and optionally BiG-SCAPE) results and convert to input used by 
	the lsaBGC framework. Will run OrthoFinder2 on just proteins from antiSMASH BGCs.

	Note, to avoid issues with BiG-SCAPE clustering (if used instead lsaBGC-Cluster.py), please use distinct output 
	prefices for each sample so that BGC names do not overlap across samples (can happen if sample genomes were
	assembled by users and do not have unique identifiers). 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--genome_listing',
                        help='Tab-delimited, two column file where the first column is the sample/isolate/genome name and the second is the full path to the genome file (Genbank or FASTA)',
                        required=True)
    parser.add_argument('-a', '--antismash_listing',
                        help='Tab-delimited, two column file where the first column is the sample/isolate/genome name the second is the full path to an antiSMASH BGC prediction in Genbank format.',
                        required=True)
    parser.add_argument('-o', '--output_directory', help="Parent output/workspace directory.", required=True)
    parser.add_argument('-b', '--bigscape_results', help="Path to BiG-SCAPE results directory.", required=False,
                        default=None)
    parser.add_argument('-gg', '--genomes_as_genbanks', action='store_true',
                        help='Genomes used for initial antiSMASH analysis were in Genbank format with CDS features which have protein translations included.',
                        default=False, required=False)
    parser.add_argument('-a', '--annotate', action='store_true',
                        help='Perform annotation of BGC proteins using KOFam HMM profiles.')
    parser.add_argument('-i', '--input_listing',
                        help='Samples to retain in the analysis. Each sample should be listed on a single line.',
                        required=False, default=None)
    parser.add_argument('-c', '--cores', type=int,
                        help="Total number of cores/threads to use for running OrthoFinder2/prodigal.", required=False,
                        default=1)

    args = parser.parse_args()
    return args


def lsaBGC_Ready():
    """
	Void function which runs primary workflow for program.
	"""

    """
	PARSE REQUIRED INPUTS
	"""
    myargs = create_parser()

    outdir = os.path.abspath(myargs.output_directory) + '/'
    genome_listing_file = os.path.abspath(myargs.genome_listing)
    antismash_listing_file = os.path.abspath(myargs.antismash_listing)

    try:
        assert (os.path.isfile(genome_listing_file))
    except:
        raise RuntimeError('Issue with path to genome listing file.')

    try:
        assert (os.path.isfile(antismash_listing_file))
    except:
        raise RuntimeError('Issue with path to antiSMASH listing file.')

    if os.path.isdir(outdir):
        sys.stderr.write("Output directory exists. Continuing in 5 seconds ...\n ")
        sleep(5)
    else:
        os.system('mkdir %s' % outdir)

    """
	PARSE OPTIONAL INPUTS
	"""

    cores = myargs.cores
    genomes_as_genbanks = myargs.genomes_as_genbanks
    input_listing_file = myargs.input_listing
    bigscape_results_dir = myargs.bigscape_results
    annotate = myargs.annotate

    if input_listing_file != None:
        input_listing_file = os.path.abspath(input_listing_file)
        try:
            assert (os.path.isfile(input_listing_file))
        except:
            raise RuntimeError(
                "Input listing file to subset genomes from more comprehensive analysis provided but does not exist.")

    if bigscape_results_dir != None:
        try:
            bigscape_results_dir = os.path.abspath(bigscape_results_dir) + '/'
            assert (os.path.isdir(bigscape_results_dir))
        except:
            raise RuntimeError('Issue with BiG-SCAPE results directory.')

    """
	START WORKFLOW
	"""

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Process Genomes

    sample_genomes, format_prediction = util.parseSampleGenomes(genome_listing_file, input_listing_file, logObject)

    if format_prediction == 'mixed':
        logObject.error('Format of genomes provided is not consistently FASTA or Genbank, please check input.')
        raise RuntimeError('Format of genomes provided is not consistently FASTA or Genbank, please check input.')

    proteomes_directory = outdir + 'Predicted_Proteomes_Initial/'
    os.system('mkdir %s' % proteomes_directory)

    if not genomes_as_genbanks:
        ## genomes are provided as FASTAs
        prodigal_outdir = outdir + 'Prodigal_Gene_Calling/'
        genbanks_directory = outdir + 'Genomic_Genbanks_Initial/'
        os.system('mkdir %s %s' % (prodigal_outdir, genbanks_directory))
        util.processGenomes(sample_genomes, prodigal_outdir, proteomes_directory, genbanks_directory, logObject,
                            cores=cores, locus_tag_length=3)
    else:
        ## genomes are provided as Genbanks
        util.extractProteins(sample_genomes, proteomes_directory, logObject)

    # Step 1: Perform annotation if requested
    if annotate:


    # Step 2: Process antiSMASH Results and Add Annotations
    antismash_bgcs_directory = outdir + 'AntiSMASH_BGCs_Retagged/'
    if not os.path.isdir(antismash_bgcs_directory): os.system('mkdir %s' % antismash_bgcs_directory)
    sample_bgcs, bgc_to_sample = util.processAntiSMASHGenbanks(antismash_listing_file, antismash_bgcs_directory,
                                                               proteomes_directory, logObject)

    # Step 3: Extract BGC proteins from full predicted proteomes
    bgc_prot_directory = outdir + 'BGC_Proteins_per_Sample/'
    if not os.path.isdir(bgc_prot_directory): os.system('mkdir %s' % bgc_prot_directory)
    sample_bgc_proteins = util.extractProteinsFromAntiSMASHBGCs(sample_bgcs, bgc_prot_directory, logObject)

    # Step 4: Run OrthoFinder2 with Proteins from BGCs
    orthofinder_directory = outdir + 'OrthoFinder2_Results/'
    orthofinder_bgc_matrix_file = util.runOrthoFinder2(bgc_prot_directory, orthofinder_directory, logObject,
                                                       cores=cores)

    # Step 5: Determine thresholds for finding genome-wide paralogs
    blast_directory = outdir + 'BLASTing_of_Ortholog_Groups/'
    if not os.path.isdir(blast_directory): os.system('mkdir %s' % blast_directory)
    samp_hg_lts, lt_to_hg, paralogy_thresholds = util.determineParalogyThresholds(orthofinder_bgc_matrix_file,
                                                                                  bgc_prot_directory, blast_directory,
                                                                                  logObject, cores=cores)

    # Step 6: Identify Paralogs and Consolidate Differences between AntiSMASH and lsaBGC-Ready Gene Calling
    final_proteomes_directory = outdir + 'Predicted_Proteomes/'
    final_genbanks_directory = outdir + 'Genomic_Genbanks/'
    if not os.path.isdir(final_genbanks_directory): os.system('mkdir %s' % final_genbanks_directory)
    if not os.path.isdir(final_proteomes_directory): os.system('mkdir %s' % final_proteomes_directory)
    util.identifyParalogsAndCreateResultFiles(samp_hg_lts, lt_to_hg, sample_bgc_proteins, paralogy_thresholds,
                                              bgc_prot_directory,
                                              blast_directory, proteomes_directory, genbanks_directory,
                                              final_proteomes_directory,
                                              final_genbanks_directory, outdir, logObject, cores=cores)

    # Step 7: Process BiG-SCAPE Results and Create GCF Listings (if provided by user)
    if bigscape_results_dir != None:
        gcf_listings_directory = outdir + 'GCF_Listings/'
        if not os.path.isdir(gcf_listings_directory): os.system('mkdir %s' % gcf_listings_directory)
        util.createGCFListingsDirectory(sample_bgcs, bigscape_results_dir, gcf_listings_directory, logObject)

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)


if __name__ == '__main__':
    lsaBGC_Ready()