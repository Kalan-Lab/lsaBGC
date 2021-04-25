# !/usr/bin/env python

### Program: lsaBGC-HMMExpansion.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

import os
import sys
from time import sleep
import argparse
from lsaBGC import processing, util
from lsaBGC.classes.GCF import GCF

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-HMMExpansion.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program takes in a list of representative BGCs belonging to a single GCF from samples which have been analyzed 
	via OrthoFinder and also a more comprehensive lists of BGCs (genbanks + proteomes) from barebones AntiSMASH 
	analysis on the comprehensive set of available genomes for the lineage/taxa of focus. The program builds profile 
	HMMs for homolog groups identified by OrthoFinder from the representative BGCs, searches for their presence in 
	the comprehensive set of BGCs (from draft genomes), and expands the OrthoFinder homolog vs. sample matrix to incorporate
	new samples and append BGCs from the draft genomes into the list of BGCs belonging to the GCF in question. 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name while the 2nd column is the path to an AntiSMASH BGC in Genbank format.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-a', '--assembly_listing', type=str, help="Tab delimited text file. First column is the sample name and the second is the path to its assembly in FASTA format. Please remove troublesome characters in the sample name.", required=True)
    parser.add_argument('-o', '--output_directory', help="Path to output directory.", required=True)
    parser.add_argument('-l', '--lineage', type=str, help="The lineage under investigation.", required=False,
                        default="Lineage")
    parser.add_argument('-cp', '--conda_path', type=str,
                        help="Path to anaconda/miniconda installation directory itself.", required=True)
    parser.add_argument('-pe', '--prokka_env_path', type=str, help="Path to conda environment for Prokka.",
                        required=True)
    parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)
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
    assembly_listing_file = os.path.abspath(myargs.assembly_listing)
    outdir = os.path.abspath(myargs.output_directory) + '/'
    conda_path = myargs.conda_path
    prokka_env_path = os.path.abspath(myargs.prokka_env_path) + '/'

    prokka_load_code = '. %s/etc/profile.d/conda.sh && conda activate %s &&' % (conda_path, prokka_env_path)

    ### vet input files quickly
    try:
        assert (os.path.isfile(orthofinder_matrix_file))
        assert (os.path.isfile(gcf_listing_file))
        assert (os.path.isfile(assembly_listing_file))
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
    lineage = myargs.lineage

    """
    START WORKFLOW
    """
    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, assembly_listing_file, outdir, gcf_id, lineage, cores]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Assembly Listing File", "Output Directory", "GCF Identifier", "Lineage", "Cores"]
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

    # Step 3: Build HMMs for homolog groups observed in representative BGCs for GCF
    logObject.info("Building profile HMMs of homolog groups observed in representative BGCs for GCF.")
    GCF_Object.constructHMMProfiles(outdir, cores=cores)
    logObject.info("HMM profiles constructed and concatenated successfully!")

    # Step 4: Process assemblies
    logObject.info("Parse sample assemblies from listing file.")
    sample_assemblies = processing.readInAssemblyListing(assembly_listing_file, logObject)
    logObject.info("Successfully parsed sample assemblies.")

    # Step 5: Run Prokka for gene calling
    locus_tag_length = 4
    prokka_outdir = outdir + 'Prokka_Results/'
    prokka_proteomes_dir = prokka_outdir + 'Prokka_Proteomes/'
    prokka_genbanks_dir = prokka_outdir + 'Prokka_Genbanks/'
    """
    try:
        if not os.path.isdir(prokka_outdir): os.system('mkdir %s' % prokka_outdir)
        if not os.path.isdir(prokka_proteomes_dir): os.system('mkdir %s' % prokka_proteomes_dir)
        if not os.path.isdir(prokka_genbanks_dir): os.system('mkdir %s' % prokka_genbanks_dir)
    except:
        logObject.error("Can't create Prokka results directories. Exiting now ...")
        raise RuntimeError("Can't create Prokka results directories. Exiting now ...")

    logObject.info("Running/setting-up Prokka for all samples!")
    processing.runProkka(sample_assemblies, prokka_outdir, prokka_proteomes_dir, prokka_genbanks_dir, prokka_load_code,
                         lineage, cores, locus_tag_length, logObject, skip_annotation_flag=True)
    logObject.info("Successfully ran/set-up Prokka.")
    """

    # Step 6: Search HMMs in proteomes from comprehensive set of BGCs
    logObject.info("Searching for homolog group HMMs in proteins extracted from comprehensive list of BGCs.")
    GCF_Object.runHMMScanAndAssignBGCsToGCF(outdir, prokka_genbanks_dir, prokka_proteomes_dir, orthofinder_matrix_file, cores=cores)

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_Expansion()