# !/usr/bin/env python

### Program: lsaBGC-PopGene.py
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
	Program: lsaBGC-PopGene.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will create iTol tracks visualizing BGCs from a single GCF across a species tree. Alternatively, if a 
	species tree is not available, it will also create a phylogeny based on single copy core genes of the GCF.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC listings file for a gcf. Tab delimited: 1st column lists sample name while the 2nd column is the path to an AntiSMASH BGC in Genbank format.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-o', '--output_directory', help="Path to output directory.", required=True)
    parser.add_argument('-a', '--codon_alignments_dir', help="Path to directory with codon alignments. Will redo if not provided.", required=False, default=None)
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
    species_phylogeny_file = os.path.abspath(myargs.species_phylogeny)
    outdir = os.path.abspath(myargs.output_directory) + '/'

    ### vet input files quickly
    try:
        assert (os.path.isfile(orthofinder_matrix_file))
        assert (os.path.isfile(gcf_listing_file))
        assert (os.path.isfile(species_phylogeny_file))
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

    cores = myargs.cores
    codon_alignments_dir = myargs.codon_alignments_dir
    population_classification_file = myargs.population_classification

    """
    START WORKFLOW
    """

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = lsaBGC.createLoggerObject(log_file)

    # Step 1: Process GCF listings file
    logObject.info("Processing BGC Genbanks from GCF listing file.")
    bbgc_gbk, bgc_genes, comp_gene_info, all_genes, bgc_sample, sample_bgcs = lsaBGC.readInBGCGenbanksPerGCF(gcf_listing_file, logObject)
    logObject.info("Successfully parsed BGC Genbanks and associated with unique IDs.")

    # Step 2: Parse OrthoFinder Homolog vs Sample Matrix
    logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
    gene_to_cog, cog_genes, cog_median_gene_counts = lsaBGC.parseOrthoFinderMatrix(orthofinder_matrix_file, all_genes)
    logObject.info("Successfully parsed homolog matrix.")

    # Step 3: Create codon alignments if not provided a directory with them (e.g. one produced by lsaBGC-See.py)
    logObject.info("User requested construction of phylogeny from SCCs in BGC! Beginning phylogeny construction.")
    if codon_alignments_dir == None:
        logObject.info("Codon alignments were not provided, so beginning process of creating protein alignments for each homolog group using mafft, then translating these to codon alignments using PAL2NAL.")
        codon_alignments_dir = lsaBGC.constructCodonAlignments(bgc_sample, cog_genes, comp_gene_info, outdir, cores, logObject)
        logObject.info("All codon alignments for SCC homologs now successfully achieved!")
    else:
        logObject.info("Codon alignments were provided by user. Moving forward to phylogeny construction with FastTree2.")

    # Step 4:
    plots_dir = outdir + 'Codon_MSA_Plots/'
    if not os.path.isdir(plots_dir): os.system('mkdir %s' % plots_dir)
    cog_stats = {}
    for f in os.listdir(codon_alignments_dir):
        cog = f.split('.msa.fna')[0]
        codon_alignment_fasta = codon_alignments_dir + f
        parsed_stats = lsaBGC.parseCodonAlignmentStats(cog, codon_alignment_fasta, plots_dir)
        cog_stats[cog] = parsed_stats

    # Step 6: Create phylogeny using FastTree2 after creating concatenated BGC alignment and processing to remove
    # sites with high rates of missing data.
    logObject.info("Creating phylogeny using FastTree2 after creating concatenated BGC alignment and processing to remove sites with high rates of missing data!")

final_output_handle = open(outdir + 'Ortholog_Group_Information.txt', 'w')
final_output_handle.write('\t'.join(
    ['cog', 'cog_order_index', 'cog_median_copy_count', 'annotation', 'is_core', 'median_gene_length', 'bgcs_with_cog',
     'samples_with_cog', 'Tajimas_D', 'core_codons', 'total_variable_codons', 'nonsynonymous_codons',
     'synonymous_codons', 'dn_ds', 'total_effective_codons', 'nonsynonymous_effective_codons',
     'synonymous_effective_codons', 'effective_dn_ds', 'cog_phylogenetic_breadth', 'differentially_variable_domains',
     'all_domains', 'coinfinder_coevolved_genes']) + '\n')
for cog, order_rank in sorted(cog_order_scores.items(), key=itemgetter(1)):
    final_output_handle.write('\t'.join([str(x) for x in [cog, order_rank, cog_median_counts[cog]] + cog_stats[cog] + [
        '; '.join(coevolved_cogs[cog])]]) + '\n')
final_output_handle.close()

if __name__ == '__main__':
    lsaBGC_PopGene()
