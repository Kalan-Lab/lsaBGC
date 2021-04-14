# !/usr/bin/env python

### Program: lsaBGC-Refiner.py
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
	Program: lsaBGC-Refiner.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will take in a list of homologous (ideally orthologous) BGC genbanks belonging to a single GCF and 
	whittles them down to include only annotations/features in between user specified homolog groups. It is particularly
	useful for curation of GCFs which featuere distinct BGCs aggregated together due to close physical proximity as 
	described in: https://msystems.asm.org/content/6/2/e00057-21/article-info 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name while the 2nd column is the path to an AntiSMASH BGC in Genbank format.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder homolog by sample matrix.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-o', '--output_directory', help="Output directory.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-b1', '--first_boundary_homolog', help="Identifier for the first homolog group to be used as boundary for pruning BGCs..", required=True)
    parser.add_argument('-b2', '--second_boundary_homolog', help="Identifier for the second homolog group to be used as boundary for pruning BGCs.", required=True)
    parser.add_argument('-c', '--cores', type=int, help="Number of cores to use for MCL step.", required=False, default=1)
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
    cores = myargs.cores
    first_boundary_homolog = myargs.first_boundary_homolog
    second_boundary_homolog = myargs.second_boundary_homolog

    """
    START WORKFLOW
    """
    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, outdir, gcf_id, first_boundary_homolog, second_boundary_homolog, cores]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Output Directory", "GCF Identifier",
                       "First Boundary Homolog", "Second Boundary Homolog", "Cores"]
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

    # Check whether boundary homolog groups are associated with GCF
    try:
        assert(first_boundary_homolog in hg_genes.keys())
        assert(second_boundary_homolog in hg_genes.keys())
    except Exception as e:
        logObject.error("Unable to determine one or both boundary homolog groups in set of homolog groups associated with GCF!")
        raise RuntimeError("Unable to determine one or both boundary homolog groups in set of homolog groups associated with GCF!")

    # Step 3: Refine BGCs
    logObject.info("Beginning refinement of BGCs!")
    new_gcf_listing_file = outdir + gcf_id + '.txt'
    GCF_Object.refineBGCGenbanks(new_gcf_listing_file, first_boundary_homolog, second_boundary_homolog)
    logObject.info("iTol track written and automatic plot via gggenes/ggtree (R) rendered!")

    # Step 5: (Optional) Create phylogeny from single-copy-core homologs from BGCs across samples (single copy in samples, not BGCs)
    if create_core_gcf_phylogeny:
        logObject.info("User requested construction of phylogeny from SCCs in BGC! Beginning phylogeny construction.")
        if codon_alignments_dir == None:
            logObject.info("Codon alignments were not provided, so beginning process of creating protein alignments for each homolog group using mafft, then translating these to codon alignments using PAL2NAL.")
            GCF_Object.constructCodonAlignments(outdir, only_scc=True, cores=cores)
            logObject.info("All codon alignments for SCC homologs now successfully achieved!")
        else:
            GCF_Object.codo_alg_dir = codon_alignments_dir
            logObject.info("Codon alignments were provided by user: %s.\nMoving forward to phylogeny construction with FastTree2." % codon_alignments_dir)

        # Step 6: Create phylogeny using FastTree2 after creating concatenated BGC alignment and processing to remove
        # sites with high rates of missing data.
        logObject.info("Creating phylogeny using FastTree2 after creating concatenated BGC alignment and processing to remove sites with high rates of missing data!")

        GCF_Object.constructGCFPhylogeny(outdir + 'BGC_SCCs_Concatenated.fasta', outdir + 'BGC_SCCs_Concatenated.nwk')
        GCF_Object.modifyPhylogenyForSamplesWithMultipleBGCs(outdir + 'BGC_SCCs_Concatenated.nwk', outdir + 'BGC_SCCs_Concatenated.edited.nwk')

        GCF_Object.visualizeGCFViaR(outdir + 'BGCs_Visualization.gggenes.txt',
                                    outdir + 'BGCs_Visualization.heatmap.txt',
                                    outdir + 'species_phylogeny.edited.nwk',
                                    outdir + 'BGC_Visualization.BGC_phylogeny.pdf')

    # Close logging object and exit
    util.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_Refiner()