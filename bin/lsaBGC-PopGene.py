# !/usr/bin/env python

### Program: lsaBGC-PopGene.py
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
from collections import defaultdict
from lsaBGC import util
from lsaBGC.classes.GCF import GCF

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-PopGene.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program investigates conservation and population genetic related statistics for each homolog group 
	observed in BGCs belonging to a single GCF.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name while the 2nd column is the path to a BGC prediction in Genbank format.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-o', '--output_directory', help="Path to output directory.", required=True)
    parser.add_argument('-k', '--sample_set', help="Sample set to keep in analysis. Should be file with one sample id per line.", required=False)
    parser.add_argument('-s', '--species_phylogeny', help="The species phylogeny in Newick format. Specifies that dN/dS should be calculated by comparing extant homolog-group sequences to ancestral state reconstruction. Does not currently work!!!", required=False, default=None)
    parser.add_argument('-u', '--population_classification', help='Popualation classifications for each sample. Tab delemited: 1st column lists sample name while the 2nd column is an identifier for the population the sample belongs to.', required=False, default=None)
    parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO).\nDefault is antiSMASH.', default='antiSMASH', required=False)
    parser.add_argument('-c', '--cpus', type=int, help="The number of cpus to use.", required=False, default=1)
    parser.add_argument('-d', '--regular_mafft', action='store_true', help="Run mafft --linsi and not the MAGUS divide-and-conquer approach which allows for scalability and more efficient computing.", default=False, required=False)
    parser.add_argument('-e', '--each_pop', action='store_true', help='Run analyses individually for each population as well.', required=False, default=False)
    parser.add_argument('-t', '--filter_for_outliers', action='store_true', help='Filter instances of homolog groups which deviate too much from the median gene length observed for the initial set of proteins.', required=False, default=False)
    parser.add_argument('-w', '--expected_similarities', help="Path to file listing expected similarities between genomes/samples. This is\ncomputed most easily by running lsaBGC-Ready.py with '-t' specified, which will estimate\nsample to sample similarities based on alignment used to create species phylogeny.", required=False, default=None)
    args = parser.parse_args()

    return args

def lsaBGC_PopGene():
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

    species_phylogeny = myargs.species_phylogeny
    sample_set_file = myargs.sample_set
    gcf_id = myargs.gcf_id
    bgc_prediction_software = myargs.bgc_prediction_software.upper()
    cpus = myargs.cpus
    population_classification_file = myargs.population_classification
    run_for_each_pop = myargs.each_pop
    filter_for_outliers = myargs.filter_for_outliers
    expected_distances = myargs.expected_similarities
    regular_mafft = myargs.regular_mafft

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
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, outdir, gcf_id, bgc_prediction_software,
                        population_classification_file, sample_set_file, species_phylogeny, run_for_each_pop,
                        filter_for_outliers, expected_distances, regular_mafft, cpus]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Output Directory", "GCF Identifier",
                       "BGC Prediction Software", "Populations Specification/Listing File", "Sample Retention Set",
                       "Species Phylogeny Newick File", "Run Analysis for Each Population",
                       "Filter for Outlier Homolog Group Instances", "File with Expected Amino Acid Differences Between Genomes/Samples",
                       "AAI Similarity Instead of ANI", "Use Regular MAFFT - not MAGUS?", "cpus"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Create GCF object
    GCF_Object = GCF(gcf_listing_file, gcf_id=gcf_id, logObject=logObject)

    # Step 0: (Optional) Parse sample set retention specifications file, if provided by the user.
    sample_retention_set = util.getSampleRetentionSet(sample_set_file)

    # Step 0B: (Optional) Parse sample to sample genome-wide relationships
    gw_pairwise_similarities = None
    if expected_distances and os.path.isfile(expected_distances):
        try:
            gw_pairwise_similarities = defaultdict(lambda: defaultdict(float))
            with open(expected_distances) as of:
                for line in of:
                    line = line.strip()
                    ls = line.split('\t')
                    s1, s2, sim = ls[:3]
                    sim = float(sim)
                    gw_pairwise_similarities[s1][s2] = sim
        except Exception as e:
            error_message = 'Had issues reading the output of expected distances results in file: %s' % expected_distances
            logObject.error(error_message)
            raise RuntimeError(error_message)

    # Step 1: Process GCF listings file
    logObject.info("Processing BGC Genbanks from GCF listing file.")
    GCF_Object.readInBGCGenbanks(comprehensive_parsing=True, prune_set=sample_retention_set,
                                 prediction_method=bgc_prediction_software)
    logObject.info("Successfully parsed BGC Genbanks and associated with unique IDs.")

    # Step 2: Parse OrthoFinder Homolog vs Sample Matrix
    logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
    gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy = util.parseOrthoFinderMatrix(orthofinder_matrix_file, GCF_Object.pan_genes)
    GCF_Object.inputHomologyInformation(gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy)
    GCF_Object.identifyKeyHomologGroups()
    logObject.info("Successfully parsed homolog matrix.")

    # Step 3: Calculate homolog order index (which can be used to roughly predict order of homologs within BGCs)
    GCF_Object.determineHgOrderIndex()

    # Step 4: (Optional) Parse population and sample inclusion specifications file, if provided by user.
    if population_classification_file:
        logObject.info("User provided information on populations, parsing this information.")
        GCF_Object.readInPopulationsSpecification(population_classification_file, prune_set=sample_retention_set)

    # Step 5: Create codon alignments if not provided a directory with them (e.g. one produced by lsaBGC-See.py)
    logObject.info("User requested construction of phylogeny from SCCs in BGC! Beginning phylogeny construction.")
    logObject.info("Beginning process of creating protein alignments for each homolog group using mafft, then translating these to codon alignments using PAL2NAL.")
    GCF_Object.constructCodonAlignments(outdir, only_scc=False, cpus=cpus, list_alignments=True, filter_outliers=False)
    if filter_for_outliers:
        GCF_Object.constructCodonAlignments(outdir, only_scc=False, cpus=cpus, list_alignments=True, filter_outliers=True, use_magus=(not regular_mafft))
    logObject.info("All codon alignments for SCC homologs now successfully achieved!")

    # Step 6: Analyze codon alignments and parse population genetics and conservation stats
    logObject.info("Beginning population genetics analyses of each codon alignment.")
    populations = [None]
    population_analysis_on = False
    if population_classification_file:
        populations = populations + list(sorted(set(GCF_Object.sample_population.values())))
        population_analysis_on = True
    if run_for_each_pop:
        for pop in populations:
            GCF_Object.runPopulationGeneticsAnalysis(outdir, cpus=cpus, population=pop, filter_outliers=False, population_analysis_on=population_analysis_on, gw_pairwise_similarities=gw_pairwise_similarities, use_translation=True, species_phylogeny=species_phylogeny)
            if filter_for_outliers:
                GCF_Object.runPopulationGeneticsAnalysis(outdir, cpus=cpus, population=pop, filter_outliers=True, population_analysis_on=population_analysis_on, gw_pairwise_similarities=gw_pairwise_similarities, use_translation=True, species_phylogeny=species_phylogeny)
    else:
        GCF_Object.runPopulationGeneticsAnalysis(outdir, cpus=cpus, population=None, filter_outliers=False, population_analysis_on=population_analysis_on, gw_pairwise_similarities=gw_pairwise_similarities, use_translation=True, species_phylogeny=species_phylogeny)
        if filter_for_outliers:
            GCF_Object.runPopulationGeneticsAnalysis(outdir, cpus=cpus, population=None, filter_outliers=True, population_analysis_on=population_analysis_on, gw_pairwise_similarities=gw_pairwise_similarities, use_translation=True, species_phylogeny=species_phylogeny)
    logObject.info("Successfully ran population genetics and evolutionary analyses of each codon alignment.")

    # Write checkpoint file for lsaBGC-AutoAnalyze.py
    checkpoint_file = outdir + 'CHECKPOINT.txt'
    checkpoint_handle = open(checkpoint_file, 'w')
    checkpoint_handle.write('lsaBGC-PopGene completed successfully!')
    checkpoint_handle.close()

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_PopGene()