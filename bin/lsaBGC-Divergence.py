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
import traceback
from time import sleep
import argparse
from collections import defaultdict
from lsaBGC import util
from lsaBGC.classes.GCF import GCF
from lsaBGC.classes.Pan import Pan

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-Divergence.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will calculate Beta-RD, the ratio of the estimated amino acid distances between orthologous BGCs from 
	two samples to the expected differences based on core protein alignments performed by requesting GToTree analysis in 
	lsaBGC-Ready, for all pairs of samples featuring a BGC belonging to a focal GCF of interest.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC specifications file. Tab delimited: 1st column contains path to BGC Genbank and 2nd column contains sample name.', required=True)
    parser.add_argument('-l', '--input_listing', help="Path to tab delimited file listing: (1) sample name (2) path to Prokka Genbank and (3) path to Prokka predicted proteome. This file is produced by lsaBGC-Process.py.", required=True)
    parser.add_argument('-a', '--codon_alignments', help="File listing the codon alignments for each homolog group in the GCF. Can be found as part of PopGene output.", required=True)
    parser.add_argument('-w', '--expected_similarities', help="Path to file listing expected similarities between genomes/samples. This is\ncomputed most easily by running lsaBGC-Ready.py with '-t' specified, which will estimate\nsample to sample similarities based on alignment used to create species phylogeny.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
    parser.add_argument('-k', '--sample_set', help="Sample set to keep in analysis. Should be file with one sample id per line.", required=False)
    parser.add_argument('-n', '--use_codon', help="Expected sample to sample similarities are reflective of DNA distances instead of protein distances (e.g. if FastANI or MASH were used in computeGenomeWideDistances.py).", required=False, default=False)
    parser.add_argument('-c', '--cpus', type=int, help="The number of cpus to use.", required=False, default=1)
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
    input_listing_file = os.path.abspath(myargs.input_listing)
    codon_alignments_file = os.path.abspath(myargs.codon_alignments)
    expected_distances = myargs.expected_similarities
    outdir = os.path.abspath(myargs.output_directory) + '/'

    ### vet input files quickly
    try:
        assert (os.path.isfile(input_listing_file))
        assert (os.path.isfile(gcf_listing_file))
        assert (os.path.isfile(codon_alignments_file))
        assert (os.path.isfile(expected_distances))

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

    use_codon_flag = myargs.use_codon
    sample_set_file = myargs.sample_set
    gcf_id = myargs.gcf_id
    cpus = myargs.cpus

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
    parameter_values = [gcf_listing_file, input_listing_file, codon_alignments_file, outdir,
                        gcf_id, sample_set_file, expected_distances, use_codon_flag, cpus]
    parameter_names = ["GCF Listing File", "Input Listing File of Prokka Annotation Files for All Samples",
                       "File Listing the Location of Codon Alignments for Each Homolog Group",
                       "Output Directory", "GCF Identifier", "Retention Sample Set",
                       "File with Expected Amino Acid Differences Between Genomes/Samples", "Use Codon Distances",
                       "cpus"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Step 0: (Optional) Parse sample set retention specifications file, if provided by the user.
    sample_retention_set = util.getSampleRetentionSet(sample_set_file)

    # Step 1: Parse expected differences results as the genome-wide similarity estimates
    gw_pairwise_similarities = defaultdict(lambda: defaultdict(float))
    try:
        with open(expected_distances) as of:
            for line in of:
                line = line.strip()
                ls = line.split('\t')
                s1, s2, sim = ls[:3]
                sim = float(sim)
                gw_pairwise_similarities[s1][s2] = sim
    except Exception as e:
        error_message = 'Had issues reading expected distances file: %s' % expected_distances
        logObject.error(error_message)
        raise RuntimeError(error_message)

    # Step 2: Determine Sequence Similarity from Codon Alignments
    logObject.info("Determining similarities in BGC content and sequence space between pairs of samples.")
    bgc_pairwise_similarities = util.determineBGCSequenceSimilarityFromCodonAlignments(codon_alignments_file, cpus=cpus, use_translation=(not use_codon_flag))
    logObject.info("Finished determining BGC specific similarity between pairs of samples.")

    # Step 3: Calculate and report Beta-RD statistic for all pairs of samples/isolates
    logObject.info("Beginning generation of report.")
    samples_bgc = set(bgc_pairwise_similarities.keys())
    try:
        final_report = outdir + 'Relative_Divergence_Report.txt'
        final_report_handle = open(final_report, 'w')
        final_report_handle.write('gcf_id\tsample_1\tsample_2\tbeta_rd\tgw_seq_sim\tgcf_seq_sim\tgcf_content_sim\n')
        for i, s1 in enumerate(samples_bgc):
            for j, s2 in enumerate(samples_bgc):
                if sample_retention_set and (not s1 in sample_retention_set or not s2 in sample_retention_set): continue
                if i >= j: continue
                gcf_seq_sim = bgc_pairwise_similarities[s1][s2][0]
                gcf_con_sim = bgc_pairwise_similarities[s1][s2][1]
                gw_seq_sim = gw_pairwise_similarities[s1][s2]

                if gw_seq_sim != 0.0:
                    if gcf_seq_sim != 'NA':
                        beta_rd = min([gcf_seq_sim/gw_seq_sim, 2.0])
                        final_report_handle.write('%s\t%s\t%s\t%f\t%f\t%f\t%f\n' % (gcf_id, s1, s2, beta_rd, gw_seq_sim, gcf_seq_sim, gcf_con_sim))
                    else:
                        logObject.warning('Samples %s and %s had no GCF-related genes/positions in common and are thus not reported!' % (s1, s2))
                else:
                    beta_rd = 2.0
                    final_report_handle.write('%s\t%s\t%s\t%f\t%f\t%f\t%f\n' % (gcf_id, s1, s2, beta_rd, gw_seq_sim, gcf_seq_sim, gcf_con_sim))
                    logObject.warning('Samples %s and %s had an estimated ANI of 0.0, thus default beta-rd of 2.0 is assigned!' % (s1, s2))

        final_report_handle.close()
        logObject.info("Successfully computed Beta-RD statistic between pairs of samples to measure BGC similarity relative to Genome-Wide similarity.")
    except Exception as e:
        error_message = "Had issues attempting to calculate Beta-RD stat between pairs of samples to measure BGC similarity relative to Genome-Wide similarity."
        logObject.error(error_message)
        logObject.error(traceback.format_exc())
        raise RuntimeError(error_message)

    # Write checkpoint file for lsaBGC-AutoAnalyze.py
    checkpoint_file = outdir + 'CHECKPOINT.txt'
    checkpoint_handle = open(checkpoint_file, 'w')
    checkpoint_handle.write('lsaBGC-Divergence completed successfully!')
    checkpoint_handle.close()

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_Divergence()