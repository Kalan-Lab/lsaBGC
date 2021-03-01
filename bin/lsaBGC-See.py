# !/usr/bin/env python

### Program: lsaBGC-See.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

import os
import sys
from time import sleep
import argparse
from lsaBGC import lsaBGC

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: lsaBGC-See.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will create iTol tracks visualizing BGCs from a single GCF across a species tree. Alternatively, if a 
	species tree is not available, it will also create a phylogeny based on single copy core genes of the GCF.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name while the 2nd column is the path to an AntiSMASH BGC in Genbank format.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-o', '--output_directory', help="Output directory.", required=True)
    parser.add_argument('-l', '--dataset_label', help="Dataset label for iTol track", required=False, default='lsaBGC-See')
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

    dataset_label = myargs.dataset_label
    species_phylogeny = myargs.species_phylogeny
    cores = myargs.cores
    create_core_gcf_phylogeny = myargs.create_core_gcf_phylogeny
    codon_alignments_dir = myargs.codon_alignments_dir

    """
    START WORKFLOW
    """

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = lsaBGC.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, outdir, dataset_label, species_phylogeny, cores,
                        create_core_gcf_phylogeny, codon_alignments_dir]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Output Directory", "Dataset Label",
                       "Species Phylogeny Newick File", "Cores", "Create GCF Phylogeny?", "Codon Alignments Directory"]
    lsaBGC.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Step 1: Process GCF listings file
    logObject.info("Processing BGC Genbanks from GCF listing file.")
    bbgc_gbk, bgc_genes, comp_gene_info, all_genes, bgc_sample, sample_bgcs = lsaBGC.readInBGCGenbanksPerGCF(gcf_listing_file, logObject)
    logObject.info("Successfully parsed BGC Genbanks and associated with unique IDs.")

    # Step 2: If species phylogeny was provided, edit it to feature duplicate leaves for isolates which have multiple
    # BGCs in the GCF.
    if species_phylogeny:
        logObject.info("Altering species phylogeny to reflect multiple BGCs per sample/isolate.")
        lsaBGC.modifyPhylogenyForSamplesWithMultipleBGCs(species_phylogeny, sample_bgcs, outdir + 'species_phylogeny.edited.nwk', logObject)
        logObject.info("Successfully edited species phylogeny.")

    # Step 3: Parse OrthoFinder Homolog vs Sample Matrix and associate each homolog group with a color
    logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
    gene_to_cog, cog_genes, cog_median_gene_counts = lsaBGC.parseOrthoFinderMatrix(orthofinder_matrix_file, all_genes)
    cog_to_color = lsaBGC.assignColorsToCOGs(gene_to_cog.values())
    logObject.info("Successfully parsed homolog matrix.")

    # Step 4: Create iTol tracks for viewing BGCs of GCF across phylogeny.
    logObject.info("Create iTol tracks for viewing BGCs of GCF across phylogeny. Note, should be used to annotate edited species phylogeny or BGC SCC phylogeny as some samples could have multiple BGCs!")
    lsaBGC.createItolBGCSeeTrack(bgc_genes, gene_to_cog, cog_to_color, comp_gene_info, dataset_label, outdir, logObject)
    logObject.info("iTol track written!")

    # Step 5: (Optional) Create phylogeny from single-copy-core homologs from BGCs across samples (single copy in samples, not BGCs)
    if create_core_gcf_phylogeny:
        logObject.info("User requested construction of phylogeny from SCCs in BGC! Beginning phylogeny construction.")
        if codon_alignments_dir == None:
            logObject.info("Codon alignments were not provided, so beginning process of creating protein alignments for each homolog group using mafft, then translating these to codon alignments using PAL2NAL.")
            codon_alignments_dir = lsaBGC.constructCodonAlignments(bgc_sample, cog_genes, comp_gene_info, outdir, cores, logObject)
            logObject.info("All codon alignments for SCC homologs now successfully achieved!")
        else:
            logObject.info("Codon alignments were provided by user. Moving forward to phylogeny construction with FastTree2.")

        # Step 6: Create phylogeny using FastTree2 after creating concatenated BGC alignment and processing to remove
        # sites with high rates of missing data.
        logObject.info("Creating phylogeny using FastTree2 after creating concatenated BGC alignment and processing to remove sites with high rates of missing data!")
        lsaBGC.constructBGCPhylogeny(codon_alignments_dir, outdir + 'BGC_SCCs_Concatenated', logObject)
        lsaBGC.modifyPhylogenyForSamplesWithMultipleBGCs(species_phylogeny, sample_bgcs, outdir + 'BGC_SCCs_Concatenated.edited.nwk', logObject)
        logObject.info("Phylogeny created successfully!")

    # Close logging object and exit
    lsaBGC.closeLoggerObject(logObject)
    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_See()