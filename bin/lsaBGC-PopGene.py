# !/usr/bin/env python

### Program: lsaBGC-PopGene.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

import os
import sys
from time import sleep
import argparse
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

    parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name while the 2nd column is the path to an AntiSMASH BGC in Genbank format.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-o', '--output_directory', help="Path to output directory.", required=True)
    parser.add_argument('-p', '--population_classification', help='Popualation classifications for each sample. Tab delemited: 1st column lists sample name while the 2nd column is an identifier for the population the sample belongs to.', required=False, default=None)
    parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)
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

    gcf_id = myargs.gcf_id
    cores = myargs.cores
    population_classification_file = myargs.population_classification

    """
    START WORKFLOW
    """
    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, outdir, gcf_id, population_classification_file, cores]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Output Directory", "GCF Identifier",
                       "Populations Specification/Listing File", "Cores"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Create GCF object
    GCF_Object = GCF(gcf_listing_file, gcf_id=gcf_id, logObject=logObject)

    # Step 1: Process GCF listings file
    logObject.info("Processing BGC Genbanks from GCF listing file.")
    GCF_Object.readInBGCGenbanks(comprehensive_parsing=True)
    logObject.info("Successfully parsed BGC Genbanks and associated with unique IDs.")

    # Step 2: Parse OrthoFinder Homolog vs Sample Matrix
    logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
    gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy = util.parseOrthoFinderMatrix(orthofinder_matrix_file, GCF_Object.pan_genes)
    GCF_Object.inputHomologyInformation(gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy)
    logObject.info("Successfully parsed homolog matrix.")

    # Step 3: Calculate homolog order index (which can be used to roughly predict order of homologs within BGCs)
    GCF_Object.determineHgOrderIndex()

    # Step 4: (Optional) Parse population specifications file, if provided by user
    sample_population = None
    if population_classification_file:
        logObject.info("User provided information on populations, parsing this information.")
        GCF_Object.readInPopulationsSpecification(population_classification_file)

    # Step 5: Create codon alignments if not provided a directory with them (e.g. one produced by lsaBGC-See.py)
    logObject.info("User requested construction of phylogeny from SCCs in BGC! Beginning phylogeny construction.")
    logObject.info("Beginning process of creating protein alignments for each homolog group using mafft, then translating these to codon alignments using PAL2NAL.")
    GCF_Object.constructCodonAlignments(outdir, only_scc=False, cores=cores, list_alignments=True)
    logObject.info("All codon alignments for SCC homologs now successfully achieved!")

    # Step 6: Analyze codon alignments and parse population genetics and conservation stats
    logObject.info("Beginning population genetics analysis of each codon alignment.")
    GCF_Object.runPopulationGeneticsAnalysis(outdir, cores=cores)
    logObject.info("")

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_PopGene()