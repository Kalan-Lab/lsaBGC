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
from lsaBGC import processing
from lsaBGC.classes.GCF import GCF

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-DiscoVary.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will construct a reference database of alleles for each homolog group within a GCF, afterwards
	it will map raw paired-end Illumina sequencing data to each non-redundant instance of a homolog group and 
	report: (i) support for different alleles of homolog groups being present in the reads and (ii) identify any
	novel SNVs.	
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC specifications file. Tab delimited: 1st column contains path to BGC Genbanks and 2nd column contains sample name.', required=True)
    parser.add_argument('-s', '--sequencing_readsets_listing', help="Sequencing data specifications file. Tab delimited: 1st column contains metagenomic sample name, whereas 2nd, 3rd, and so on columns contain full paths to sequencing reads.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-l', '--input_listing', type=str, help="Tab delimited text file for samples with three columns: (1) sample name\n(2) path to whole-genome generated Genbank file (*.gbk), and (3)path to whole-genome generated\npredicted-proteome file (*.faa).", required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
    parser.add_argument('-a', '--codon_alignments', help="File listing the codon alignments for each homolog group in the GCF.\nCan be found as part of PopGene output.", required=True)
    parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is antiSMASH.', default='antiSMASH', required=False)
    parser.add_argument('-awl', '--ambiguity_window_length', help='Length around beginning, end and ambiguous sites in codon alignments to avoid calling SNVs on due to more challenging alignment for reads. Default is 50.', required=False, default=50)
    parser.add_argument('-ch', '--core_homologs', nargs="+", help="List of homolog group identifiers comprising the core of the BGC/GCF.", required=False, default=[])
    parser.add_argument('-ap', '--allow_phasing', action='store_true', help="Allow phasing with DESMAN. Requires manual installation of\nDESMAN (not through conda) and $PATH to be updated.", required=False, default=False)
    parser.add_argument('-c', '--cpus', type=int, help="The number of cpus to use.", required=False, default=1)

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
    input_listing_file = os.path.abspath(myargs.input_listing)
    paired_end_sequencing_file = os.path.abspath(myargs.sequencing_readsets_listing)
    orthofinder_matrix_file = os.path.abspath(myargs.orthofinder_matrix)
    codon_alignments_file = os.path.abspath(myargs.codon_alignments)
    outdir = os.path.abspath(myargs.output_directory) + '/'

    ### vet input files quickly
    try:
        assert (os.path.isfile(orthofinder_matrix_file))
        assert (os.path.isfile(gcf_listing_file))
        assert (os.path.isfile(input_listing_file))
        assert (os.path.isfile(paired_end_sequencing_file))
        assert (os.path.isfile(codon_alignments_file))
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

    bgc_prediction_software = myargs.bgc_prediction_software.upper()
    gcf_id = myargs.gcf_id
    allow_phasing = myargs.allow_phasing
    core_hg_set = myargs.core_homologs
    ambiguity_window_length = myargs.ambiguity_window_length
    cpus = myargs.cpus

    try:
        assert (bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
    except:
        raise RuntimeError('BGC prediction software option is not a valid option.')

    """
    START WORKFLOW
    """
    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)
    version_string = util.parseVersionFromSetupPy()
    logObject.info('Running lsaBGC version %s' % version_string)

    # Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future records.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, paired_end_sequencing_file, outdir, gcf_id,
                        bgc_prediction_software, ambiguity_window_length, len(core_hg_set) > 0, allow_phasing, cpus]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Paired-Sequencing Listing File",
                       "Output Directory", "GCF Identifier", "BGC Prediction Software", "Ambiguity Window Length",
                       "Core Homologs Manually Provided", "Allow DESMAN Phasing?", "cpus"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Create GCF object
    GCF_Object = GCF(gcf_listing_file, gcf_id=gcf_id, logObject=logObject)

    # Step 1: Process GCF listings file
    logObject.info("Processing BGC Genbanks from GCF listing file.")
    GCF_Object.readInBGCGenbanks(comprehensive_parsing=True, prediction_method=bgc_prediction_software)
    logObject.info("Successfully parsed BGC Genbanks and associated with unique IDs.")

    # Step 2: Parse OrthoFinder Homolog vs Sample Matrix and associate each homolog group with a color
    logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
    gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy = util.parseOrthoFinderMatrix(orthofinder_matrix_file, GCF_Object.pan_genes)
    GCF_Object.inputHomologyInformation(gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy)
    GCF_Object.identifyKeyHomologGroups()
    if len(core_hg_set) > 0:
        GCF_Object.core_homologs = set(core_hg_set)
    logObject.info("Successfully parsed homolog matrix.")

    # Step 3: Process annotation files related to input sample sets
    logObject.info("Parsing annotation files provided for set of samples to incorporate into analysis.")
    input_sample_prokka_data = processing.readInAnnotationFilesForExpandedSampleSet(input_listing_file, logObject=logObject)
    logObject.info("Successfully parsed new sample annotation files.")

    # Step 4: Build HMMs for homolog groups observed in representative BGCs for GCF
    ### This function is originally intended for the expansion functionality in lsaBGC, but here we
    ### use it to just get a better sense of how paralogous different genes in the GCF might be.
    logObject.info("Determining non-unique positions along codon multiple sequence alignments.")
    hg_nonunique_positions = util.determineNonUniqueRegionsAlongCodonAlignment(outdir, input_sample_prokka_data, codon_alignments_file, cpus=cpus, logObject=logObject)
    logObject.info("Marked non-unique positions along codon MSAs!")

    # Step 5: Create database of genes with surrounding flanks and, independently, cluster them into allele groups / haplotypes.
    logObject.info("Extracting and clustering GCF genes with their flanks.")
    genes_representative_fasta = outdir + 'GCF_Genes_Representatives.fasta'
    genes_fasta = outdir + 'GCF_Genes.fasta'
    bowtie2_db_prefix = outdir + 'GCF_Genes_Representatives'
    GCF_Object.extractGenesAndCluster(genes_representative_fasta, genes_fasta, codon_alignments_file, bowtie2_db_prefix)
    logObject.info("Successfully extracted genes and clustered them into discrete alleles.")

    # Step 6: Align paired-end reads to database genes with surrounding flanks
    bowtie2_outdir = outdir + 'Bowtie2_Alignments/'
    if not os.path.isfile(bowtie2_outdir): os.system('mkdir %s' % bowtie2_outdir)
    logObject.info("Running Bowtie2 alignment of paired-end sequencing reads against database of GCF genes with surrounding flanking sequences.")
    util.runBowtie2Alignments(bowtie2_db_prefix, paired_end_sequencing_file, bowtie2_outdir, logObject, cpus=cpus)
    logObject.info("Bowtie2 alignments completed successfully!")

    # Step 7: Parse bowtie2 alignments found per sample and identify support for SNVs
    results_outdir = outdir + 'Alignment_Parsing/'
    if not os.path.isdir(results_outdir): os.system('mkdir %s' % results_outdir)
    logObject.info("Beginning typing of homolog group alleles and mining of novel SNVs.")
    GCF_Object.runSNVMining(paired_end_sequencing_file, genes_representative_fasta, codon_alignments_file, bowtie2_outdir, results_outdir, cpus=cpus)
    logObject.info("Successfully typed alleles and mined for novel SNVs.")

    # Step 8: Decide on GCF presence, determine consensus/haplotypes for homolog groups, and generate novelty report
    # and create matrix of most closely resembling reference alleles
    phased_alleles_outdir = outdir + 'Phased_Homolog_Group_Sequences/'
    if not os.path.isdir(phased_alleles_outdir): os.system('mkdir %s' % phased_alleles_outdir)
    logObject.info("Phasing or determining consensus allele and reporting of novel SNVs.")
    GCF_Object.phaseAndSummarize(paired_end_sequencing_file, codon_alignments_file, results_outdir, phased_alleles_outdir, outdir, hg_nonunique_positions, cpus=cpus, allow_phasing=allow_phasing, ambiguity_window_length=ambiguity_window_length)
    logObject.info("Successfully constructed matrices of allele typings.")

    # Step 9: Filter low coverage gene instances and construct gene-phylogenies
    comp_hg_phylo_outdir = outdir + 'Comprehensive_Homolog_Group_Phylogenies/'
    if not os.path.isdir(comp_hg_phylo_outdir): os.system('mkdir %s' % comp_hg_phylo_outdir)
    logObject.info("Filtering low coverage gene instances and construct gene-phylogenies.")
    GCF_Object.generateGenePhylogenies(codon_alignments_file, phased_alleles_outdir, comp_hg_phylo_outdir, hg_nonunique_positions)
    logObject.info("Successfully generated gene-specific phylogenies.")

    # Step 10: Determine similarity in BGC content between pairs of samples.
    logObject.info("Determining pairwise differences in BGC content between samples.")
    GCF_Object.calculatePairwiseDifferences(paired_end_sequencing_file, results_outdir, outdir)
    logObject.info("Successfully calculated pairwise differences between samples.")

    # Write checkpoint file for lsaBGC-AutoAnalyze.py
    checkpoint_file = outdir + 'CHECKPOINT.txt'
    checkpoint_handle = open(checkpoint_file, 'w')
    checkpoint_handle.write('lsaBGC-DiscoVary completed successfully!')
    checkpoint_handle.close()

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_DiscoVary()
