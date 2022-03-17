# !/usr/bin/env python

### Program: lsaBGC-HMMExpansion.py
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
from lsaBGC import processing, util
from lsaBGC.classes.GCF import GCF

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Expansion.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program takes in a list of representative BGCs belonging to a single GCF (for instance those identified in 
	high quality genomic assemblies) and then searches for new instances of the GCF in an expansion set of draft 
	assemblies. 
	
	The program builds profile HMMs for homolog groups identified by OrthoFinder from the representative BGCs, 
	searches for their presence in the comprehensive set of BGCs (from draft genomes), applys a classical HMM to 
	identify neighborhoods of genes on draft assembly scaffolds which match the GCF profile, and performs conditional
	checks to filter out false-positive findings.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name while the 2nd column is the path to an AntiSMASH BGC in Genbank format.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-l', '--initial_listing', type=str, help="Tab delimited text file for samples with three columns: (1) sample name (2) Prokka generated Genbank file (*.gbk), and (3) Prokka generated predicted-proteome file (*.faa). Please remove troublesome characters in the sample name.", required=True)
    parser.add_argument('-e', '--expansion_listing', type=str, help="Tab delimited text file for samples in the expansion set with three columns: (1) sample name (2) Prokka generated Genbank file (*.gbk), and (3) Prokka generated predicted-proteome file (*.faa). Please remove troublesome characters in the sample name.", required=True)
    parser.add_argument('-o', '--output_directory', help="Path to output directory.", required=True)
    parser.add_argument('-ms', '--min_segment_size', type=float, help="The minimum number of homolog groups needed to report discrete segments of the GCF. Default is 5.", required=False, default=5)
    parser.add_argument('-mcs', '--min_segment_core_size', type=float, help="The minimum number of core (to the known instances of the GCF) homololog groups needed to report discrete segments of the GCF. Default is 3.", required=False, default=3)
    parser.add_argument('-sct', '--syntenic_correlation_threshold', type=float, help="The minimum syntenic correlation needed to at least one known GCF instance to report segment.", required=False, default=0.8)
    parser.add_argument('-tg', '--transition_from_gcf_to_gcf', type=float, help="GCF to GCF state transition probability for HMM. Should be between 0.0 and 1.0. Default is 0.9.", required=False, default=0.9)
    parser.add_argument('-tb', '--transition_from_bg_to_bg', type=float, help="Background to background state transition probability for HMM. Should be between 0.0 and 1.0. Default is 0.9.", required=False, default=0.9)
    parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)
    parser.add_argument('-q', '--quick_mode', action='store_true', help="Run in quick-mode. Instead of running HMMScan for each homolog group, a consensus sequence is emitted and Diamond is used for searching instead. Method inspired by Melnyk et al. 2019", required=False, default=False)
    parser.add_argument('-no', '--no_orthogroup_matrix', action='store_true', help="Avoid writing the updated OrthoFinder matrix at the end.", required=False, default=False)
    parser.add_argument('-z', '--pickle_expansion_annotation_data', help="Pickle file with serialization of annotation data in the expansion listing file.", required=False, default=None)

    args = parser.parse_args()

    return args

def lsaBGC_Expansion():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """
    myargs = create_parser()

    gcf_listing_file = os.path.abspath(myargs.gcf_listing)
    orthofinder_matrix_file = os.path.abspath(myargs.orthofinder_matrix)
    initial_listing_file = os.path.abspath(myargs.initial_listing)
    expansion_listing_file = os.path.abspath(myargs.expansion_listing)
    outdir = os.path.abspath(myargs.output_directory) + '/'

    ### vet input files quickly
    try:
        assert (os.path.isfile(orthofinder_matrix_file))
        assert (os.path.isfile(gcf_listing_file))
        assert (os.path.isfile(initial_listing_file))
        assert (os.path.isfile(expansion_listing_file))
    except:
        raise RuntimeError('One or more of the input files provided, does not exist. Exiting now ...')

    if os.path.isdir(outdir):
        sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n")
        sleep(5)
    else:
        os.system('mkdir %s' % outdir)

    """
    PARSE OPTIONAL INPUTS
    """

    gcf_id = myargs.gcf_id
    cores = myargs.cores
    min_segment_size = myargs.min_segment_size
    min_segment_core_size = myargs.min_segment_core_size
    syntenic_correlation_threshold = myargs.syntenic_correlation_threshold
    transition_from_gcf_to_gcf = myargs.transition_from_gcf_to_gcf
    transition_from_bg_to_bg = myargs.transition_from_bg_to_bg
    quick_mode = myargs.quick_mode
    no_orthogroup_matrix = myargs.no_orthogroup_matrix
    pickle_expansion_annotation_data_file = myargs.pickle_expansion_annotation_data

    """
    START WORKFLOW
    """
    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, initial_listing_file, expansion_listing_file, outdir,
                        gcf_id, cores, min_segment_size, min_segment_core_size, syntenic_correlation_threshold,
                        transition_from_gcf_to_gcf, transition_from_bg_to_bg, quick_mode, no_orthogroup_matrix,
                        pickle_expansion_annotation_data_file]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File",
                       "Listing File of Prokka Annotation Files for Initial Set of Samples",
                       "Listing File of Prokka Annotation Files for Expansion/Additional Set of Samples",
                       "Output Directory", "GCF Identifier", "Cores", "Minimum Size of Segments",
                       "Minimum Core Size of Segments", "Syntenic Correlation Threshold",
                       "HMM Transition Probability from GCF to GCF",
                       "HMM Transition Probability from Background to Background", "Run Expansion in Quick Mode?",
                       "Skip rewriting Expanded OrthoGroup CSV File?",
                       "Pickle File with Annotation Data in Expansion Listing for Quick Loading"]
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
    GCF_Object.identifyKeyHomologGroups()
    logObject.info("Successfully parsed homolog matrix.")

    # Step 3: Process annotation files related to input and expanded sample sets
    logObject.info("Parsing annotation file provided in expansion listing file for larger set of samples to incorporate into analysis.")
    initial_sample_prokka_data = processing.readInAnnotationFilesForExpandedSampleSet(initial_listing_file, logObject=logObject)
    expanded_sample_prokka_data = processing.readInAnnotationFilesForExpandedSampleSet(expansion_listing_file, logObject=logObject)
    logObject.info("Successfully parsed new sample annotation files.")

    # Step 4: Build HMMs for homolog groups observed in representative BGCs for GCF
    logObject.info("Building profile HMMs of homolog groups observed in representative BGCs for GCF.")
    GCF_Object.constructHMMProfiles(outdir, initial_sample_prokka_data, cores=cores, quick_mode=quick_mode)
    logObject.info("HMM profiles constructed and concatenated successfully!")

    # Step 5: Search HMM profiles in proteomes from comprehensive set of BGCs
    logObject.info("Searching for homolog group HMMs in proteins extracted from comprehensive list of BGCs.")
    GCF_Object.runHMMScan(outdir, expanded_sample_prokka_data, cores=cores, quick_mode=quick_mode,
                          annotation_pickle_file=pickle_expansion_annotation_data_file)
    logObject.info("Successfully found new instances of GCF in new sample set.")

    # Step 6: Determine whether samples' assemblies feature GCF of interest
    logObject.info("Searching for homolog group HMMs in proteins extracted from comprehensive list of BGCs.")
    GCF_Object.identifyGCFInstances(outdir, expanded_sample_prokka_data, orthofinder_matrix_file,
                                    min_size=min_segment_size, min_core_size=min_segment_core_size,
                                    gcf_to_gcf_transition_prob=transition_from_gcf_to_gcf,
                                    background_to_background_transition_prob=transition_from_bg_to_bg,
                                    syntenic_correlation_threshold=syntenic_correlation_threshold,
                                    no_orthogroup_matrix=no_orthogroup_matrix, cores=cores)
    logObject.info("Successfully found new instances of GCF in new sample set.")

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_Expansion()