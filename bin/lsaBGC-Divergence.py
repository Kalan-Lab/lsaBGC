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

	This program will calculate Beta-RD, the ratio of the estimated ANI between orthologous BGCs from two samples to the
	estimated genome-wide ANI, for all pairs of samples featuring a BGC belonging to a focal GCF of interest.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC specifications file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
    parser.add_argument('-l', '--input_listing', help="Path to tab delimited file listing: (1) sample name (2) path to Prokka Genbank and (3) path to Prokka predicted proteome. This file is produced by lsaBGC-Process.py.", required=False, default=None)
    parser.add_argument('-a', '--codon_alignments', help="File listing the codon alignments for each homolog group in the GCF. Can be found as part of PopGene output.", required=True)
    parser.add_argument('-pr', '--precomputed_ani_result', help="Path to MASH or FastANI results.", required=False)
    parser.add_argument('-pi', '--precomputed_ani_input', help="Path to two-column, tab-delimited file showing FASTA used as input for MASH/FastANI in second column with sample name in the first column.", required=False)
    parser.add_argument('-f', '--fastani', action='store_true', help="Use FastANI instead of MASH for estimating genome-wide ANI between pairs of samples. Takes much longer but is more accurate.", required=False)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
    parser.add_argument('-k', '--sample_set', help="Sample set to keep in analysis. Should be file with one sample id per line.", required=False)
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
    input_listing_file = os.path.abspath(myargs.input_listing)
    codon_alignments_file = os.path.abspath(myargs.codon_alignments)
    outdir = os.path.abspath(myargs.output_directory) + '/'

    ### vet input files quickly
    try:
        assert (os.path.isfile(input_listing_file))
        assert (os.path.isfile(gcf_listing_file))
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

    sample_set_file = myargs.sample_set
    gcf_id = myargs.gcf_id
    cores = myargs.cores
    sketch_size = myargs.sketch_size
    precomputed_ani_result_file = myargs.precomputed_ani_result
    precomputed_ani_input_file = myargs.precomputed_ani_input
    used_fastani = myargs.fastani

    """
    START WORKFLOW
    """

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, input_listing_file, codon_alignments_file, outdir, sketch_size, used_fastani,
                        gcf_id, sample_set_file, precomputed_ani_result_file, precomputed_ani_input_file, cores]
    parameter_names = ["GCF Listing File", "Input Listing File of Prokka Annotation Files for All Samples",
                       "File Listing the Location of Codon Alignments for Each Homolog Group",
                       "Output Directory", "MASH Sketch Size", 'Used FastANI for ANI Estimation?',
                       "GCF Identifier", "Retention Sample Set", "Precomputed MASH/FastANI Results File",
                       "Precomputed MASH/FastANI Input File", "Cores"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Create GCF object
    GCF_Object = GCF(gcf_listing_file, gcf_id=gcf_id, logObject=logObject)

    # Step 0: (Optional) Parse sample set retention specifications file, if provided by the user.
    sample_retention_set = util.getSampleRetentionSet(sample_set_file)

    # Step 1: Extract Genbank Sequences into FASTA and Run MASH Analysis Between Genomic Assemblies
    gw_pairwise_differences = None
    if not os.path.isfile(precomputed_ani_result_file) or not os.path.isfile(precomputed_ani_input_file):
        logObject.info("Converting BGC Genbanks from GCF listing file into FASTA per sample.")
        gcf_fasta_listing_file = outdir + 'GCF_Listings.fasta'
        gcf_fasta_dir = outdir + 'Sample_GCF_FASTAs/'
        if not os.path.isdir(gcf_fasta_dir): os.system('mkdir %s' % gcf_fasta_dir)
        GCF_Object.convertGenbanksIntoFastas(gcf_fasta_dir, gcf_fasta_listing_file)
        logObject.info("Successfully performed conversion and partitioning by sample.")

        # Extract Genbank Sequences into FASTA
        Pan_Object = Pan(input_listing_file, logObject=logObject)
        logObject.info("Converting Genbanks from Expansion listing file into FASTA per sample.")
        gw_fasta_dir = outdir + 'Sample_Expansion_FASTAs/'
        if not os.path.isdir(gw_fasta_dir): os.system('mkdir %s' % gw_fasta_dir)
        gw_fasta_listing_file = outdir + 'Genome_FASTA_Listings.txt'
        Pan_Object.convertGenbanksIntoFastas(gw_fasta_dir, gw_fasta_listing_file)
        logObject.info("Successfully performed conversions.")

        logObject.info("Running MASH Analysis Between Genomes.")
        gw_pairwise_differences = util.calculateMashPairwiseDifferences(gw_fasta_listing_file, outdir, 'genome_wide', sketch_size, cores, logObject)
        logObject.info("Ran MASH Analysis Between Genomes.")
    else:
        fasta_to_name = {}
        gw_pairwise_similarities = defaultdict(lambda: defaultdict(float))
        try:
            with open(precomputed_mash_input_file) as oflf:
                for line in oflf:
                    line = line.strip()
                    ls = line.split('\t')
                    fasta_to_name[ls[1]] = ls[0]
        except:
            error_message = "Had issues reading the FASTA listing file %s" % precomputed_mash_input_file
            logObject.error(error_message)
            raise RuntimeError(error_message)
        try:
            with open(precomputed_mash_result_file) as of:
                for line in of:
                    line = line.strip()
                    ls = line.split('\t')
                    f1, f2, dist = ls[:3]
                    dist = float(dist)
                    n1 = fasta_to_name[f1]
                    n2 = fasta_to_name[f2]
                    gw_pairwise_differences[n1][n2] = dist
        except Exception as e:
            error_message = 'Had issues reading the output of MASH dist analysis in file: %s' % precomputed_mash_result_file
            logObject.error(error_message)
            raise RuntimeError(error_message)

    # Step 2: Determine Sequence Similarity from Codon Alignments
    logObject.info("Determining similarities in BGC content and sequence space between pairs of samples.")
    bgc_pairwise_similarities = util.determineBGCSequenceSimilarityFromCodonAlignments(codon_alignments_file)
    logObject.info("Finished determining BGC specific similarity between pairs of samples.")

    # Step 3: Calculate and report Beta-RD statistic for all pairs of samples/isolates
    logObject.info("Beginning generation of report.")
    samples_bgc = set(bgc_pairwise_similarities.keys())
    samples_gw = set(gw_pairwise_differences.keys())
    samples_intersect = samples_gw.intersection(samples_bgc)
    try:
        final_report = outdir + 'Relative_Divergence_Report.txt'
        final_report_handle = open(final_report, 'w')
        final_report_handle.write('gcf_id\tsample_1\tsample_2\tbeta_rd\tgw_seq_sim\tgcf_seq_sim\tgcf_content_sim\n')
        for i, s1 in enumerate(samples_intersect):
            for j, s2 in enumerate(samples_intersect):
                if sample_retention_set and (not s1 in sample_retention_set or not s2 in sample_retention_set): continue
                if i >= j: continue
                gcf_seq_sim = bgc_pairwise_similarities[s1][s2][0]
                gcf_con_sim = bgc_pairwise_similarities[s1][s2][1]
                #gcf_sim = 1.0 - gcf_dist
                gw_dist = gw_pairwise_differences[s1][s2]
                gw_seq_sim = 1.0 - gw_dist

                if gw_seq_sim != 0.0:
                    if gcf_seq_sim != 'NA':
                        beta_rd = gcf_seq_sim/gw_seq_sim
                        final_report_handle.write('%s\t%s\t%s\t%f\t%f\t%f\t%f\n' % (gcf_id, s1, s2, beta_rd, gw_seq_sim, gcf_seq_sim, gcf_con_sim))
                    else:
                        logObject.warning('Samples %s and %s had no genes/positions in common and are thus not reported!' % (s1, s2))
                else:
                    logObject.warning('Samples %s and %s had an estimated ANI of 0.0 and are thus not reported!' % (s1, s2))
        final_report_handle.close()
        logObject.info("Successfully computed Beta-RD statistic between pairs of samples to measure BGC similarity relative to Genome-Wide similarity.")
    except Exception as e:
        error_message = "Had issues attempting to calculate Beta-RD stat between pairs of samples to measure BGC similarity relative to Genome-Wide similarity."
        logObject.error(error_message)
        logObject.error(traceback.format_exc())
        raise RuntimeError(error_message)

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_Divergence()