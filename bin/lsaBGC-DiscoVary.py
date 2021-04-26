# !/usr/bin/env python

### Program: lsaBGC-DiscoVary.py
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
from time import sleep
import argparse
from lsaBGC import util
from lsaBGC.classes.GCF import GCF

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-DiscoVary.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will construct a reference database of alleles for each homolog group within a GCF, afterwards
	it will map raw paired-end Illumina sequencing data to 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC specifications file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
    parser.add_argument('-p', '--paired_end_sequencing', help="Sequencing data specifications file. Tab delimited: 1st column contains metagenomic sample name, whereas 2nd and 3rd columns contain full paths to forward and reverse reads, respectively.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
    parser.add_argument('-a', '--codon_alignments', help="File listing the codon alignments for each homolog group in the GCF. Can be found as part of PopGene output.", required=True)
    parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)

    args = parser.parse_args()

    return args

def lsaBGC_DiscoVary():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """
    myargs = create_parser()

    gcf_listing_file = os.path.abspath(myargs.gcf_listing)
    paired_end_sequencing_file = os.path.abspath(myargs.paired_end_sequencing)
    orthofinder_matrix_file = os.path.abspath(myargs.orthofinder_matrix)
    codon_alignments_file = os.path.abspath(myargs.codon_alignments)
    outdir = os.path.abspath(myargs.output_directory) + '/'

    ### vet input files quickly
    try:
        assert (os.path.isfile(orthofinder_matrix_file))
        assert (os.path.isfile(gcf_listing_file))
        assert (os.path.isfile(paired_end_sequencing_file))
    except:
        raise RuntimeError('One or more of the input files provided, does not exist. Exiting now ...')

    if os.path.isdir(outdir):
        sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
        sleep(5)
    else:
        os.system('mkdir %s' % outdir)

    """
    PARSE OPTIONAL INPUTS
    """

    gcf_id = myargs.gcf_id
    cores = myargs.cores

    """
    START WORKFLOW
    """
    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, paired_end_sequencing_file, outdir, gcf_id, cores]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Paired-Sequencing Listing File",
                       "Output Directory", "GCF Identifier", "Cores"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Create GCF object
    GCF_Object = GCF(gcf_listing_file, gcf_id=gcf_id, logObject=logObject)

    # Step 1: Process GCF listings file
    logObject.info("Processing BGC Genbanks from GCF listing file.")
    GCF_Object.readInBGCGenbanks(comprehensive_parsing=True)
    logObject.info("Successfully parsed BGC Genbanks and associated with unique IDs.")

    # Step 2: Parse OrthoFinder Homolog vs Sample Matrix and associate each homolog group with a color
    logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
    gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy = util.parseOrthoFinderMatrix(orthofinder_matrix_file, GCF_Object.pan_genes)
    GCF_Object.inputHomologyInformation(gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy)
    logObject.info("Successfully parsed homolog matrix.")

    # Step 3: Create database of genes with surrounding flanks and, independently, cluster them into allele groups / haplotypes.
    logObject.info("Extracting and clustering GCF genes with their flanks.")

    genes_with_flanks_fasta = outdir + 'GCF_Genes.fasta'
    cd_hit_clusters_fasta_file = outdir + 'GCF_Genes_Clusters.fasta'
    cd_hit_nr_fasta_file = outdir + 'GCF_Genes_NR.fasta'
    bowtie2_db_prefix = outdir + 'GCF_Genes'

    GCF_Object.extractGeneWithFlanksAndCluster(genes_with_flanks_fasta, cd_hit_clusters_fasta_file, cd_hit_nr_fasta_file, bowtie2_db_prefix)

    logObject.info("Successfully extracted genes with flanks and clustered them into discrete haplotypes.")

    # Step 4: Align paired-end reads to database genes with surrounding flanks
    bowtie2_outdir = outdir + 'Bowtie2_Alignments/'
    if not os.path.isfile(bowtie2_outdir): os.system('mkdir %s' % bowtie2_outdir)
    logObject.info("Running Bowtie2 alignment of paired-end sequencing reads against database of GCF genes with surrounding flanking sequences.")
    util.runBowtie2Alignments(bowtie2_db_prefix, paired_end_sequencing_file, bowtie2_outdir, logObject, cores=cores)
    logObject.info("Bowtie2 alignments completed successfully!")

    # Step 5: Determine haplotypes found in samples and identify supported novelty SNVs
    results_outdir = outdir + 'SNV_Miner_and_Allele_Typing_Results/'
    if not os.path.isdir(results_outdir): os.system('mkdir %s' % results_outdir)
    logObject.info("Beginning typing of homolog group alleles and mining of novel SNVs.")
    GCF_Object.runSNVMining(paired_end_sequencing_file, cd_hit_nr_fasta_file, bowtie2_outdir, results_outdir, cores=cores)
    logObject.info("Successfully typed alleles and mined for novel SNVs.")

    # Step 6: Construct summary matrices
    logObject.info("Consolidating allele typing results into matrix formats.")
    GCF_Object.createSummaryMatricesForMetaNovelty(paired_end_sequencing_file, results_outdir, outdir)
    logObject.info("Successfully constructed matrices of allele typings.")

    # Step 7: Create Novelty Report
    logObject.info("Generating report of novel SNVs found across paire-end sequencing reads.")
    GCF_Object.generateNoveltyReport(codon_alignments_file, results_outdir, outdir)
    logObject.info("Successfully generated novelty SNV report.")

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_DiscoVary()
