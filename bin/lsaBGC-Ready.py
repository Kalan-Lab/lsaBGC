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

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Ready.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Program to convert existing BGC predictions, e.g. from antiSMASH, DeepBGC, and GECCO, (and optionally BiG-SCAPE) 
	results and convert to input used by the lsaBGC suite (make it "ready" for lsaBGC analysis). Will run OrthoFinder2 
	on just proteins from antiSMASH BGCs. If BiG-SCAPE results are not provided, users have the option to run 
	lsaBGC-Cluster instead which implements algorithms designed for clustering complete instances of BGCs from 
	completed/finished genomic assemblies.

    There are scripts from creating input listing files from AntiSMASH result directories and directories with 
    genomes (recursivelyIdentifyBGCGenbanks.py & listAllGenomesInDirectory.py) if you don't want 
    to write them yourself. The listing files as inputs instead of directories help ensure that sample mapping 
    between genomes (FASTA or Genbank) and BGC Genbanks and should be manually investigated to ensure proper
    linking.
    
	ALGORITHMIC OVERVIEWS/CONSIDERATIONS:
	***************************************************************************************************************** 
    -*-  OrthoFinder2 modes:
            * Genome_Wide: Run OrthoFinder2 as intended with all primary sample full genome-wide proteomes.
              [LOW-THROUGHPUT (<200 Genomes)]
            * BGC_Only: OrthoFinder2 is run across samples/genomes accounting for only BGC embedded proteins. 
              Genome-wide paralogs for orthogroups are subsequently identified by using orthogroup specific cutoffs
              based on the percent identity and coverage thresholds determined for each orthogroup (the minimum 
              perc. id and coverage observed within BGC proteins belonging to the same orthgroup). 
              [DEFAULT; MEDIUM-THROUGHPUT (>200 but <500 genomes)]
            * COMING SOON: palo - scalable genome-wide orthology determination.
	   
    -*-  To avoid issues with processing BiG-SCAPE results (if used instead lsaBGC-Cluster.py), please use distinct 
         output prefices for each sample when running antiSMASH so that BGC names do not overlap across samples 
         (can happen if sample genomes were assembled by users and do not have unique identifiers). If issues persist
         please consider using lsaBGC-Cluster.py, we use similar methods to BiG-SCAPE, though the algorithms are 
         mostly designed for complete BGCs in mind for lsaBGC-Cluster.py, while BiG-SCAPE has some nice settings 
         to handle fragmented BGCs. lsaBGC-Expansion/AutoExpansion are specifically designed for detecting fragmented
         GCF instances in draft assemblies and can be run on the initial "primary" genome set as well.
	***************************************************************************************************************** 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--genome_listing',
                        help='Tab-delimited, two column file for primary samples (ideally with high-quality or complete genomes)\nwhere the first column is the sample/isolate/genome name and the second is the\nfull path to the genome file (Genbank or FASTA).\nCheck note above about available scripts to automatically create this.',
                        required=True)
    parser.add_argument('-d', '--additional_genome_listing', help='Tab-delimited, two column file for samples with additional/draft\ngenomes (same format as for the "--genome_listing" argument). The genomes/BGCs of these\nsamples won\'t be used in ortholog-grouping of proteins and clustering of BGCs, but will simply have gene\ncalling run for them. This will enable more sensitive/expanded detection of GCF instances later\nusing lsaBGC-Expansion/AutoExpansion.\nCheck note above about available scripts to automatically create this.',
                        required=False, default=None)
    parser.add_argument('-l', '--bgc_genbank_listing',
                        help='Tab-delimited, two column file listing BGC predictions results for primary samples\n(those from the "--genome_listing" argument), where the first column is the sample name and the second\nis the full path to BGC prediction in Genbank format.',
                        required=True)
    parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is antiSMASH.', default='antiSMASH', required=False)
    parser.add_argument('-b', '--bigscape_results', help='Path to BiG-SCAPE results directory of antiSMASH/DeepBGC/GECCO results predicted in primary\ngenomes.Please make sure the sample names match what is provided for "--genome_listings".', required=False,
                        default=None)
    parser.add_argument('-o', '--output_directory', help='Parent output/workspace directory.', required=True)
    parser.add_argument('-m', '--orthofinder_mode', help='Method for running OrthoFinder2. (Options: Genome_Wide, BGC_Only). Default is BGC_Only.', required=False, default='BGC_Only')
    parser.add_argument('-a', '--annotate', action='store_true',
                        help='Perform annotation of BGC proteins using KOfam HMM profiles.', required=False, default=False)
    parser.add_argument('-lc', '--lsabgc_cluster', action='store_true', help='Run lsaBGC-Cluster with default parameters. Note, we recommend running lsaBGC-Cluster manually\nand exploring parameters through its ability to generate a user-report for setting clustering parameters.', required=False, default=False)
    parser.add_argument('-le', '--lsabgc_expansion', action='store_true', help='Run lsaBGC-AutoExpansion with default parameters. Assumes either "--bigscape_results" or\n"--lsabgc_cluster" is specified.', default=False, required=False)
    parser.add_argument('-c', '--cores', type=int,
                        help="Total number of cores/threads to use for running OrthoFinder2/prodigal.", required=False,
                        default=1)
    parser.add_argument('-k', '--keep_intermediates', action='store_true', help='Keep intermediate directories / files which are likely not useful for downstream analyses.', required=False, default=False)
    parser.add_argument('-spe', '--skip_primary_expansion', action='store_true', help='Skip expansion on primary genomes as well.', required=False, default=False)
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
    bgc_genbank_listing_file = os.path.abspath(myargs.bgc_genbank_listing)

    try:
        assert (os.path.isfile(genome_listing_file))
    except:
        raise RuntimeError('Issue with reading genome listing file for samples with complete genomic assemblies.')

    try:
        assert (os.path.isfile(bgc_genbank_listing_file))
    except:
        raise RuntimeError('Issue with path to BGC predictions Genbanks listing file.')

    if os.path.isdir(outdir):
        sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
        sleep(5)
    else:
        os.system('mkdir %s' % outdir)

    int_outdir = outdir + 'Intermediate_Files/'
    fin_outdir = outdir + 'Final_Results/'
    util.setupReadyDirectory([int_outdir, fin_outdir])

    """ 
	PARSE OPTIONAL INPUTS
	"""

    additional_genome_listing_file = myargs.additional_genome_listing
    bgc_prediction_software = myargs.bgc_prediction_software.upper()
    orthofinder_mode = myargs.orthofinder_mode
    cores = myargs.cores
    bigscape_results_dir = myargs.bigscape_results
    annotate = myargs.annotate
    run_lsabgc_cluster = myargs.lsabgc_cluster
    run_lsabgc_expansion = myargs.lsabgc_expansion
    keep_intermediates = myargs.keep_intermediates
    skip_primary_expansion = myargs.skip_primary_expansion

    try:
        assert(bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
    except:
        raise RuntimeError('BGC prediction software option is not a valid option.')

    if additional_genome_listing_file != None:
        try:
            additional_genome_listing_file = os.path.abspath(additional_genome_listing_file)
            assert (os.path.isfile(additional_genome_listing_file))
        except:
            raise RuntimeError('Issue with reading genome listing file for samples with additional genomic assemblies.')

    if bigscape_results_dir != None:
        try:
            bigscape_results_dir = os.path.abspath(bigscape_results_dir) + '/'
            assert (os.path.isdir(bigscape_results_dir))
        except:
            raise RuntimeError('Issue with BiG-SCAPE results directory.')

    kofam_hmm_file = None
    kofam_pro_list = None
    if annotate:
        try:
            kofam_db_location = lsaBGC_main_directory + 'db/kofam_location_paths.txt'
            assert(os.path.isfile(kofam_db_location))
            with open(kofam_db_location) as okdlf:
                for line in okdlf:
                    kofam_pro_list, kofam_hmm_file = line.strip().split('\t')
                    assert(os.path.isfile(kofam_hmm_file) and os.path.isfile(kofam_pro_list))
        except:
            raise RuntimeError("It appears KOfam database was not setup or setup successfully. Please run/rerun the script setup_annotation_dbs.py and report to Github issues if issue persists.")

    """
	START WORKFLOW
	"""

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)
    logObject.info("Saving parameters for future records.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [genome_listing_file, bgc_genbank_listing_file, outdir, additional_genome_listing_file,
                        bgc_prediction_software, orthofinder_mode, bigscape_results_dir, annotate, run_lsabgc_cluster,
                        run_lsabgc_expansion, keep_intermediates, cores]
    parameter_names = ["Primary Genome Listing File", "BGC Predictions Genbanks Listing File", "Output Directory",
                       "Additional Genome Listing File", "BGC Prediction Software", "OrthoFinder Mode",
                       "BiG-SCAPE Results Directory", "Perform KOfam Annotation?", "Run lsaBGC-Cluster Analysis?",
                       "Run lsaBGC-AutoExpansion Analysis?", "Keep Intermediate Files/Directories?", "Number of Cores"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Step 1: Process Primary Genomes
    sample_genomes, format_prediction = util.parseSampleGenomes(genome_listing_file, logObject)

    if format_prediction == 'mixed':
        logObject.error('Format of genomes provided is not consistently FASTA or Genbank, please check input.')
        raise RuntimeError('Format of genomes provided is not consistently FASTA or Genbank, please check input.')

    if len(sample_genomes) > 300:
        logObject.warning('Over 300 primary genomes specified! Consider de-replicating further as this will produce a lot of small files and require certain environmental conditions to be set to allow parallel reading of files. Waiting 5 seconds before proceeding.')
        sleep(5)
    if len(sample_genomes) > 500:
        logObject.error('Currently the number of primary genomes allowed is 500, which will produce 250,000 files during all-vs-all alignment by OrthoFinder2. Please dereplicate the primary number of genomes further to avoid excesive file creation and runtime.')
        raise RuntimeError('Currently the number of primary genomes allowed is 500, which will produce 250,000 files during all-vs-all alignment by OrthoFinder2. Please dereplicate the primary number of genomes further to avoid excesive file creation and runtime.')

    # Check ulimit settings!

    try:
        num_genomes = len(sample_genomes)
        num_files = num_genomes*num_genomes
        uH = subprocess.check_output('ulimit -Hn', shell=True)
        uS = subprocess.check_output('ulimit -Sn', shell=True)
        uH = int(uH.decode('utf-8').strip())
        uS = int(uS.decode('utf-8').strip())
        if num_files > uS:
            logObject.error("Too many files will be produced and need to be read at once by OrthoFinder2. Luckily, you can resolve this quite easily via: ulimit -n 250000 . After running the command rerun lsaBGC-Ready.py and you should not get stuck here again.")
            raise RuntimeError("Too many files will be produced and need to be read at once by OrthoFinder2. Luckily, you can resolve this quite easily via: ulimit -n 250000 . After running the command rerun lsaBGC-Ready.py and you should not get stuck here again.")
        if num_files > uH:
            logObject.error("Too many files will be produced and need to be read at once by OrthoFinder2. Your system requires root privs to change this, which I do not recommend. See the following OrthoFinder2 Github issue for more details: https://github.com/davidemms/OrthoFinder/issues/384")
            raise RuntimeError("Too many files will be produced and need to be read at once by OrthoFinder2. Your system requires root privs to change this, which I do not recommend. See the following OrthoFinder2 Github issue for more details: https://github.com/davidemms/OrthoFinder/issues/384")
        if num_files < uS and num_files < uH:
            logObject.info("Great news! You have correctly set ulimit settings and we believe OrthoFinder2 should be able to run smoothly!")
    except:
        logObject.error("Difficulties validating ulimit settings are properly set to allow for successful OrthoFinder2 run.")
        raise RuntimeError("Difficulties validating ulimit settings are properly set to allow for successful OrthoFinder2 run.")

    proteomes_directory = outdir + 'Predicted_Proteomes_Initial/'
    util.setupReadyDirectory([proteomes_directory])

    if format_prediction == 'fasta':
        ## genomes are provided as FASTAs
        prodigal_outdir = outdir + 'Prodigal_Gene_Calling/'
        genbanks_directory = outdir + 'Genomic_Genbanks_Initial/'
        util.setupReadyDirectory([prodigal_outdir, genbanks_directory])

        # Note, locus tags of length 3 are used within lsaBGC to mark samples whose BGCs were integral in defining GCFs.
        util.processGenomes(sample_genomes, prodigal_outdir, proteomes_directory, genbanks_directory, logObject,
                            cores=cores, locus_tag_length=3)
        sample_genomes = util.updateSampleGenomesWithGenbanks(genbanks_directory)
    else:
        ## genomes are provided as Genbanks with CDS features
        util.extractProteins(sample_genomes, proteomes_directory, logObject)

    # Step 2: Process Additional Genomes
    additional_sample_annotation_listing_file = int_outdir + 'Additional_Sample_Annotation_Files.txt'
    if additional_genome_listing_file != None:
        additional_sample_genomes, additional_format_prediction = util.parseSampleGenomes(additional_genome_listing_file, logObject)
        if additional_format_prediction == 'mixed':
            logObject.error('Format of additional genomes provided is not consistently FASTA or Genbank, please check input.')
            raise RuntimeError('Format of additional genomes provided is not consistently FASTA or Genbank, please check input.')

        additional_proteomes_directory = outdir + 'Predicted_Proteomes_Additional/'
        util.setupReadyDirectory([additional_proteomes_directory])
        if additional_format_prediction == 'fasta':
            additional_prodigal_outdir = outdir + 'Prodigal_Gene_Calling_Additional/'
            additional_genbanks_directory = outdir + 'Genomic_Genbanks_Additional/'
            util.setupReadyDirectory([additional_prodigal_outdir, additional_genbanks_directory])

            # Note, locus tags of length 4 are used within lsaBGC to mark samples with additional genomes where we ultimately
            # find them via lsaBGC-Expansion.
            util.processGenomes(additional_sample_genomes, additional_prodigal_outdir, additional_proteomes_directory, additional_genbanks_directory,
                                logObject, cores=cores, locus_tag_length=4)
        else:
            # genomes are provided as Genbanks with CDS features
            util.extractProteins(additional_sample_genomes, additional_proteomes_directory, logObject)

        additional_sample_annotation_listing_handle = open(additional_sample_annotation_listing_file, 'w')
        for f in os.listdir(additional_proteomes_directory):
            sample = f.split('.faa')[0]
            if additional_format_prediction == 'fasta':
                additional_sample_annotation_listing_handle.write(sample + '\t' + additional_genbanks_directory + sample + '.gbk' + '\t' + additional_proteomes_directory + f + '\n')
            else:
                additional_sample_annotation_listing_handle.write(sample + '\t' + additional_sample_genomes[sample] + '\t' + additional_proteomes_directory +f + '\n')
        additional_sample_annotation_listing_handle.close()

    # Step 3: Process BGC Genbank Results and Add Annotations
    bgcs_directory = outdir + 'BGCs_Retagged/'
    util.setupReadyDirectory([bgcs_directory])

    if bgc_prediction_software == 'DEEPBGC':
        deepbgc_split_directory = outdir + 'Split_DeepBGC_Genbanks/'
        util.setupReadyDirectory([deepbgc_split_directory])
        bgc_genbank_listing_file = util.splitDeepBGCGenbank(bgc_genbank_listing_file, deepbgc_split_directory, outdir, logObject)
    sample_bgcs, bgc_to_sample = util.processBGCGenbanks(bgc_genbank_listing_file, bgc_prediction_software,
                                                         sample_genomes, bgcs_directory, proteomes_directory, logObject)

    # Step 4: Extract BGC proteins from full predicted proteomes
    bgc_prot_directory = outdir + 'BGC_Proteins_per_Sample/'
    util.setupReadyDirectory([bgc_prot_directory])

    sample_bgc_proteins = util.extractProteinsFromBGCs(sample_bgcs, bgc_prot_directory, logObject)

    # Step 5: Perform KOfam annotation if requested and update BGCs (including references to them)
    protein_annotations = None
    if annotate:
        ko_annot_directory = outdir + 'KOfam_Annotations/'
        if not os.path.isdir(ko_annot_directory): os.system('mkdir %s' % ko_annot_directory)
        protein_annotations = util.performKOFamAnnotation(sample_bgc_proteins, bgc_prot_directory,
                                                                       ko_annot_directory, kofam_hmm_file,
                                                                       kofam_pro_list, logObject, cores=cores)
        bgcs_directory_updated = outdir + 'BGCs_Retagged_and_Annotated/'
        util.setupReadyDirectory([bgcs_directory_updated])

        sample_bgcs_update, bgc_to_sample_update, sample_bgc_proteins_update = util.updateBGCGenbanksToIncludeAnnotations(protein_annotations, bgc_to_sample, sample_bgc_proteins,
                                                                                                    bgcs_directory,
                                                                                                    bgcs_directory_updated,
                                                                                                    logObject)
        bgcs_directory = bgcs_directory_updated
        sample_bgcs = sample_bgcs_update
        bgc_to_sample = bgc_to_sample_update
        sample_bgc_proteins = sample_bgc_proteins_update

    orthofinder_directory = outdir + 'OrthoFinder2_Results/'

    final_proteomes_directory = outdir + 'Predicted_Proteomes/'
    final_genbanks_directory = outdir + 'Genomic_Genbanks/'
    util.setupReadyDirectory([final_proteomes_directory, final_genbanks_directory])
    primary_sample_annotation_listing_file = int_outdir + 'Primary_Sample_Annotation_Files.txt'
    primary_bgc_listing_file = int_outdir + 'Primary_Sample_BGC_Files.txt'
    primary_orthofinder_matrix_file = int_outdir + 'Orthogroups.tsv'

    if orthofinder_mode.upper() == 'GENOME_WIDE':
        # Step 6 - GW: Incorporate proteins only found in BGC Genbanks into genome-wide predicted proteomes + genbanks
        util.incorporateBGCProteinsIntoProteomesAndGenbanks(sample_bgc_proteins, sample_genomes, protein_annotations,
                                                            bgc_prot_directory, proteomes_directory, final_proteomes_directory,
                                                            final_genbanks_directory, primary_sample_annotation_listing_file,
                                                            primary_bgc_listing_file, logObject)

        # Step 7 - GW: Run OrthoFinder2 with genome-wide predicted proteomes
        orthofinder_bgc_matrix_file = util.runOrthoFinder2(final_proteomes_directory, orthofinder_directory, logObject, cores=cores)
        os.system('mv %s %s' % (orthofinder_bgc_matrix_file, primary_orthofinder_matrix_file))
    elif orthofinder_mode.upper() == 'BGC_ONLY':
        # Step 6 - BO: Run OrthoFinder2 with Proteins from BGCs
        orthofinder_bgc_matrix_file = util.runOrthoFinder2(bgc_prot_directory, orthofinder_directory, logObject, cores=cores)

        # Step 7 - BO: Determine thresholds for finding genome-wide paralogs
        blast_directory = outdir + 'BLASTing_of_Ortholog_Groups/'
        util.setupReadyDirectory([blast_directory])
        samp_hg_lts, lt_to_hg, paralogy_thresholds = util.determineParalogyThresholds(orthofinder_bgc_matrix_file,
                                                                                      bgc_prot_directory, blast_directory,
                                                                                      logObject, cores=cores)
        # Step 8 - BO: Identify Paralogs and Consolidate Differences between BGC prediction software and lsaBGC-Ready Gene Calling
        util.identifyParalogsAndCreateResultFiles(samp_hg_lts, lt_to_hg, sample_bgc_proteins, paralogy_thresholds,
                                                  protein_annotations, bgc_prot_directory, blast_directory,
                                                  proteomes_directory, sample_genomes, final_proteomes_directory,
                                                  final_genbanks_directory, primary_sample_annotation_listing_file,
                                                  primary_bgc_listing_file, primary_orthofinder_matrix_file, logObject,
                                                  cores=cores)

    prim_samps_with_bgcs = set([])
    additional_lines_to_append = []
    with open(primary_sample_annotation_listing_file) as opsalf:
        for line in opsalf:
            line = line.strip()
            s,gbk,faa = line.split('\t')
            prim_samps_with_bgcs.add(s)
            if additional_genome_listing_file != None and not skip_primary_expansion:
                additional_lines_to_append.append(line)

    # handle some primary genomes not having any BGC predictions!
    primary_sample_annotation_listing_handle = open(primary_sample_annotation_listing_file, 'a+')
    for f in os.listdir(proteomes_directory):
        s = f.split('.faa')[0]
        if s in prim_samps_with_bgcs: continue
        faa = proteomes_directory + f
        fin_faa = final_proteomes_directory + f
        gbk = genbanks_directory + s + '.gbk'
        fin_gbk = final_genbanks_directory + s + '.gbk'
        os.system('cp %s %s' % (faa, fin_faa))
        os.system('cp %s %s' % (gbk, fin_gbk))
        primary_sample_annotation_listing_handle.write(s + '\t' + fin_gbk + '\t' + fin_faa + '\n')
        if additional_genome_listing_file != None and not skip_primary_expansion:
            additional_lines_to_append.append(s + '\t' + fin_gbk + '\t' + fin_faa)

    primary_sample_annotation_listing_handle.close()

    if additional_genome_listing_file != None and not skip_primary_expansion:
        additional_sample_annotation_listing_handle = open(additional_sample_annotation_listing_file, 'a+')
        additional_sample_annotation_listing_handle.write('\n'.join(additional_lines_to_append) + '\n')
        additional_sample_annotation_listing_handle.close()

    # Step 9: Process BiG-SCAPE Results and Create GCF Listings (if provided by user) or Run lsaBGC-Cluster if requested.
    gcf_listings_directory = None
    if bigscape_results_dir != None:
        bigscape_reformat_directory = outdir + 'BiG_SCAPE_Results_Reformatted/'
        gcf_listings_directory = bigscape_reformat_directory + 'GCF_Listings/'
        if not os.path.isdir(bigscape_reformat_directory) or not os.path.isdir(gcf_listings_directory):
            util.setupReadyDirectory([bigscape_reformat_directory, gcf_listings_directory])
        util.createGCFListingsDirectory(sample_bgcs, bgc_to_sample, bigscape_results_dir, gcf_listings_directory,
                                        logObject)
    elif run_lsabgc_cluster:
        lsabgc_cluster_results_dir = outdir + 'lsaBGC_Cluster_Results/'
        lsabgc_cluster_cmd = ['lsaBGC-Cluster.py', '-b', primary_bgc_listing_file, '-m', primary_orthofinder_matrix_file,
                              '-c', str(cores), '-o', lsabgc_cluster_results_dir, '-r', '0.7', '-i', '4.0', '-j', '20.0']
        try:
            subprocess.call(' '.join(lsabgc_cluster_cmd), shell=True, stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                            executable='/bin/bash')
            if logObject:
                logObject.info('Successfully ran: %s' % ' '.join(lsabgc_cluster_cmd))
            gcf_listings_directory = lsabgc_cluster_results_dir + 'GCF_Listings/'
        except:
            if logObject:
                logObject.error('Had an issue running: %s' % ' '.join(lsabgc_cluster_cmd))
                logObject.error(traceback.format_exc())
            raise RuntimeError('Had an issue running: %s' % ' '.join(lsabgc_cluster_cmd))

    # Step 10: Run lsaBGC-AutoExpansion if requested
    if run_lsabgc_expansion and ((bigscape_results_dir != None and os.path.isdir(bigscape_results_dir)) or
                                 (run_lsabgc_cluster and os.path.isdir(lsabgc_cluster_results_dir))) and \
                                additional_genome_listing_file != None and os.path.isdir(gcf_listings_directory):
        lsabgc_expansion_results_dir = outdir + 'lsaBGC_AutoExpansion_Results/'
        lsabgc_expansion_cmd = ['lsaBGC-AutoExpansion.py', '-g', gcf_listings_directory, '-m',
                               int_outdir + 'Orthogroups.tsv', '-l', int_outdir + 'Primary_Sample_Annotation_Files.txt',
                               '-e', int_outdir + 'Additional_Sample_Annotation_Files.txt', '-q', '-c', str(cores),
                               '-o', lsabgc_expansion_results_dir]
        try:
            subprocess.call(' '.join(lsabgc_expansion_cmd), shell=True, stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                            executable='/bin/bash')
            if logObject:
                logObject.info('Successfully ran: %s' % ' '.join(lsabgc_expansion_cmd))
        except:
            if logObject:
                logObject.error('Had an issue running: %s' % ' '.join(lsabgc_expansion_cmd))
                logObject.error(traceback.format_exc())
            raise RuntimeError('Had an issue running: %s' % ' '.join(lsabgc_expansion_cmd))

    # Step 11: Create Final Results Directory
    if not keep_intermediates:
        util.selectFinalResultsAndCleanUp(outdir, fin_outdir, logObject)

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_Ready()