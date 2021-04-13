# !/usr/bin/env python

### Program: lsaBGC-See.py
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
	Program: lsaBGC-See.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will create automatic visuals depicting genes across a species or BGC-specific phylogeny as well as
	iTol tracks visualizing BGCs from a single GCF across a species tree. Alternatively, if a species tree is not 
	available, it will also create a phylogeny based on single copy core genes of the GCF.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name while the 2nd column is the path to an AntiSMASH BGC in Genbank format.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-o', '--output_directory', help="Output directory.", required=True)
    parser.add_argument('-i', '--gcf_id', help="GCF identifier.", required=False, default='GCF_X')
    parser.add_argument('-s', '--species_phylogeny', help="The species phylogeny in Newick format.", required=False, default=None)
    parser.add_argument('-c', '--cores', type=int, help="Number of cores to use for MCL step.", required=False, default=1)
    parser.add_argument('-p', '--create_core_gcf_phylogeny', action='store_true', help="Create phylogeny from core COGs.", required=False, default=False)
    parser.add_argument('-a', '--codon_alignments_dir', help="Path to directory with codon alignments. Will redo if not provided.", required=False, default=None)
    args = parser.parse_args()
    return args

def lsaBGC_See():
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
    species_phylogeny = myargs.species_phylogeny
    cores = myargs.cores
    create_core_gcf_phylogeny = myargs.create_core_gcf_phylogeny
    codon_alignments_dir = myargs.codon_alignments_dir

    """
    START WORKFLOW
    """

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = util.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, outdir, gcf_id, species_phylogeny, cores,
                        create_core_gcf_phylogeny, codon_alignments_dir]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Output Directory", "GCF Identifier",
                       "Species Phylogeny Newick File", "Cores", "Create GCF Phylogeny?", "Codon Alignments Directory"]
    util.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Create GCF object
    GCF_Object = GCF(gcf_listing_file, gcf_id=gcf_id, logObject=logObject)

    # Step 1: Process GCF listings file
    logObject.info("Processing BGC Genbanks from GCF listing file.")
    GCF_Object.readInBGCGenbanks(comprehensive_parsing=True)
    logObject.info("Successfully parsed BGC Genbanks and associated with unique IDs.")

    # Step 2: If species phylogeny was provided, edit it to feature duplicate leaves for isolates which have multiple
    # BGCs in the GCF.
    if species_phylogeny:
        logObject.info("Altering species phylogeny to reflect multiple BGCs per sample/isolate.")
        GCF_Object.modifyPhylogenyForSamplesWithMultipleBGCs(species_phylogeny, outdir + 'species_phylogeny.edited.nwk')
        logObject.info("Successfully edited species phylogeny.")

    # Step 3: Parse OrthoFinder Homolog vs Sample Matrix and associate each homolog group with a color
    logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
    gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy = util.parseOrthoFinderMatrix(orthofinder_matrix_file, GCF_Object.pan_genes)
    GCF_Object.inputHomologyInformation(gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy)
    GCF_Object.assignColorsToHGs(GCF_Object.gene_to_hg, GCF_Object.bgc_genes)
    logObject.info("Successfully parsed homolog matrix.")

    # Step 4: Create iTol and gggenes (R) tracks for visualizing BGCs of GCF across a phylogeny.
    logObject.info("Create iTol tracks for viewing BGCs of GCF across phylogeny. Note, should be used to annotate edited species phylogeny or BGC SCC phylogeny as some samples could have multiple BGCs!")
    GCF_Object.createItolBGCSeeTrack(outdir + 'BGCs_Visualization.iTol.txt')
    GCF_Object.visualizeGCFViaR(outdir + 'BGCs_Visualization.gggenes.txt',
                                outdir + 'BGCs_Visualization.heatmap.txt',
                                outdir + 'species_phylogeny.edited.nwk',
                                outdir + 'BGC_Visualization.species_phylogeny.pdf')
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
    lsaBGC_See()