# !/usr/bin/env python

### Program: lsaBGC-Refiner.py
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
	Program: lsaBGC-Refiner.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will take in a list of homologous (ideally orthologous) BGC genbanks belonging to a single GCF and 
	whittle them down to include only annotations/features in between user specified homolog groups. It is particularly
	useful for curation of GCFs which featuere distinct BGCs aggregated together due to close physical proximity as 
	described in: https://msystems.asm.org/content/6/2/e00057-21/article-info 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name while the 2nd column is the path to a BGC prediction in Genbank format.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder homolog by sample matrix.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-o', '--output_directory', help="Output directory.", required=True)
    parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is antiSMASH.', default='antiSMASH', required=False)
    parser.add_argument('-b1', '--first_boundary_homolog', help="Identifier for the first homolog group to be used as boundary for pruning BGCs..", required=True)
    parser.add_argument('-b2', '--second_boundary_homolog', help="Identifier for the second homolog group to be used as boundary for pruning BGCs.", required=True)
    args = parser.parse_args()
    return args

def lsaBGC_Refiner():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """
    myargs = create_parser()

    gcf_listing_file = os.path.abspath(myargs.gcf_listing)
    orthofinder_matrix_file = os.path.abspath(myargs.orthofinder_matrix)
    outdir = os.path.abspath(myargs.output_directory) + '/'

    ### vet input files quickly
    try:
        assert (os.path.isfile(orthofinder_matrix_file))
        assert (os.path.isfile(gcf_listing_file))
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
    bgc_prediction_software = myargs.bgc_prediction_software.upper()
    first_boundary_homolog = myargs.first_boundary_homolog
    second_boundary_homolog = myargs.second_boundary_homolog


    try:
        assert(bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
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

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future records.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, outdir, gcf_id, bgc_prediction_software,
                        first_boundary_homolog, second_boundary_homolog]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Output Directory", "GCF Identifier",
                       "BGC Prediction Software", "First Boundary Homolog", "Second Boundary Homolog"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Create GCF object
    GCF_Object = GCF(gcf_listing_file, gcf_id=gcf_id, logObject=logObject)

    # Step 1: Process GCF listings file
    logObject.info("Processing BGC Genbanks from GCF listing file.")
    GCF_Object.readInBGCGenbanks(comprehensive_parsing=True, prediction_method=bgc_prediction_software)
    logObject.info("Successfully parsed BGC Genbanks and associated with unique IDs.")

    # Step 2: Parse OrthoFinder Homolog vs Sample Matrix
    logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
    gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy = util.parseOrthoFinderMatrix(orthofinder_matrix_file, GCF_Object.pan_genes)
    GCF_Object.inputHomologyInformation(gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy)
    logObject.info("Successfully parsed homolog matrix.")

    # Check whether boundary homolog groups are associated with GCF
    try:
        assert(first_boundary_homolog in hg_genes.keys())
        assert(second_boundary_homolog in hg_genes.keys())
    except Exception as e:
        logObject.error("Unable to determine one or both boundary homolog groups in set of homolog groups associated with GCF!")
        raise RuntimeError("Unable to determine one or both boundary homolog groups in set of homolog groups associated with GCF!")

    # Step 3: Refine BGCs based on user specifications.
    logObject.info("Beginning refinement of BGCs.")
    new_gcf_listing_file = outdir + gcf_id + '.txt'
    GCF_Object.refineBGCGenbanks(new_gcf_listing_file, outdir, first_boundary_homolog, second_boundary_homolog)
    logObject.info("Refinement of BGC Genbanks complete!")

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_Refiner()