# !/usr/bin/env python

### Program: lsaBGC-Divergence.py
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
	Program: lsaBGC-Divergence.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will calculate Beta-RD, the ratio of the estimated ANI between orthologous BGCs from two samples to the
	estimated genome-wide ANI, for all pairs of samples featuring a BGC belonging to a focal GCF of interest.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC specifications file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
    parser.add_argument('-a', '--assembly_listing', help="Sequencing data specifications file. Tab delimited: 1st column contains metagenomic sample name, whereas 2nd and 3rd columns contain full paths to forward and reverse reads, respectively.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
    parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)
    parser.add_argument('-s', '--sketch_size', type=int, help="The sketch size, number of kmers to use in fingerprinting", required=False, default=10000)

    args = parser.parse_args()

    return args

def lsaBGC_Divergence():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """
    myargs = create_parser()

    gcf_listing_file = os.path.abspath(myargs.gcf_listing)
    assembly_listing_file = os.path.abspath(myargs.assembly_listing_file)
    outdir = os.path.abspath(myargs.output_directory) + '/'

    ### vet input files quickly
    try:
        assert (os.path.isfile(assembly_listing_file))
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
    cores = myargs.cores
    sketch_size = myargs.sketch_size

    """
    START WORKFLOW
    """

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, assembly_listing_file, outdir, sketch_size, gcf_id, cores]
    parameter_names = ["GCF Listing File", "Assembly Listing File", "Output Directory", "MASH Sketch Size", "GCF Identifier", "Cores"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Create GCF object
    GCF_Object = GCF(gcf_listing_file, gcf_id=gcf_id, logObject=logObject)

    # Step 1: Extract Genbank Sequences into FASTA
    logObject.info("Converting BGC Genbanks from GCF listing file into FASTA per sample.")
    gcf_fasta_listing_file = outdir + 'GCF_Listings.fasta'
    GCF_Object.convertGenbanksIntoFastas(gcf_listing_file, gcf_fasta_listing_file)
    logObject.info("Successfully performed conversion and partitioning by sample.")

    # Step 2: Run MASH Analysis Between BGCs in GCF
    logObject.info("Running MASH Analysis Between BGCs in GCF.")
    bgc_pairwise_differences = util.calculateMashPairwiseDifferences(gcf_fasta_listing_file, outdir, 'gcf', sketch_size, cores, logObject)
    logObject.info("Ran MASH Analysis Between BGCs in GCF.")

    # Step 3: Run MASH Analysis Between Genomic Assemblies
    logObject.info("Running MASH Analysis Between Genomes.")
    gw_pairwise_differences = util.calculateMashPairwiseDifferences(assembly_listing_file, outdir, 'genome_wide', sketch_size, cores, logObject)
    logObject.info("Ran MASH Analysis Between Genomes.")

    # Step 4: Calculate and report Beta-RD statistic for all pairs of samples/isolates
    logObject.info("Beginning generation of report.")
    samples_gcf = gcf_pairwise_differences.keys()
    samples_gw = gw_pairwise_differences.keys()
    samples_intersect = samples_gw.intersection(samples_gcf)
    try:
        final_report = outdir + 'relative_divergence_report.txt'
        final_report_handle = open(final_report, 'w')
        final_report_handle.write('gcf_id\tsample_1\tsample_2\tbeta_rd\n')
        for i, s1 in enumerate(samples_intersect):
            for j, s2 in enumerate(samples_intersect):
                if i < j:
                    gcf_dist = bgc_pairwise_differences[s1][s2]
                    gcf_sim = 1.0 - gcf_dist
                    gw_dist = gw_pairwise_differences[s1][s2]
                    gw_sim = 1.0 - gw_dist

                    if gw_sim != 0.0:
                        beta_rd = gcf_sim/gw_sim
                        final_report_handle.write('%s\t%s\t%s\t%f\n' % (gcf_id, s1, s2, beta_rd))
                    else:
                        logObject.warning('Samples %s and %s had an estimated ANI of 0.0 and thus not reported!' % (s1, s2))
        final_report_handle.close()
        logObject.info("Successfully computed Beta-RD statistic between pairs of samples to measure BGC similarity relative to Genome-Wide similarity.")
    except:
        error_message = "Had issues attempting to calculate Beta-RD stat between pairs of samples to measure BGC similarity relative to Genome-Wide similarity."
        logObject.error(error_message)
        raise RuntimeError(error_message)

    # Close logging object and exit
    lsaBGC.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_Divergence()