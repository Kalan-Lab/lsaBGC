# !/usr/bin/env python

### Program: lsaBGC-MetaNovelty.py
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
	Program: lsaBGC-MetaNovelty.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	This program will construct a reference database of alleles for each homolog group within a GCF, afterwards
	it will map raw paired-end Illumina sequencing data to 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--gcf_listing', help='BGC specifications file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
    parser.add_argument('-i', '--paired_end_sequencing', help="Sequencing data specifications file. Tab delimited: 1st column contains metagenomic sample name, whereas 2nd and 3rd columns contain full paths to forward and reverse reads, respectively.", required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
    parser.add_argument('-s', '--popgene_stats', help="Resulting output from lsaBGC-PopGene for focal GCF.")
    parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)
    args = parser.parse_args()

    return args

def lsaBGC_MetaNovelty():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """
    myargs = create_parser()

    gcf_listing_file = os.path.abspath(myargs.gcf_listing)
    paired_end_sequencing_file = os.path.abspath(myargs.paired_end_sequencing)
    orthofinder_matrix_file = os.path.abspath(myargs.orthofinder_matrix)
    outdir = os.path.abspath(myargs.output_directory) + '/'

    ### vet input files quickly
    try:
        assert (os.path.isfile(orthofinder_matrix_file))
        assert (os.path.isfile(gcf_listing_file))
        assert (os.path.isfile(paired_end_sequencing_file))
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

    """
    START WORKFLOW
    """

    # create logging object
    log_file = outdir + 'Progress.log'
    logObject = lsaBGC.createLoggerObject(log_file)

    # Step 0: Log input arguments and update reference and query FASTA files.
    logObject.info("Saving parameters for future provedance.")
    parameters_file = outdir + 'Parameter_Inputs.txt'
    parameter_values = [gcf_listing_file, orthofinder_matrix_file, paired_end_sequencing_file, outdir, cores]
    parameter_names = ["GCF Listing File", "OrthoFinder Orthogroups.csv File", "Paired-Sequencing Listing File",
                       "Output Directory", "Cores"]
    lsaBGC.logParametersToFile(parameters_file, parameter_names, parameter_values)
    logObject.info("Done saving parameters!")

    # Step 1: Process GCF listings file
    logObject.info("Processing BGC Genbanks from GCF listing file.")
    bgc_gbk, bgc_genes, comp_gene_info, all_genes, bgc_sample, sample_bgcs = lsaBGC.readInBGCGenbanksPerGCF(gcf_listing_file, logObject)
    logObject.info("Successfully parsed BGC Genbanks and associated with unique IDs.")

    # Step 2: Parse OrthoFinder Homolog vs Sample Matrix and associate each homolog group with a color
    logObject.info("Starting to parse OrthoFinder homolog vs sample information.")
    gene_to_cog, cog_genes, cog_median_gene_counts = lsaBGC.parseOrthoFinderMatrix(orthofinder_matrix_file, all_genes)
    logObject.info("Successfully parsed homolog matrix.")

    # Step 3: Create database of genes with surrounding flanks and, independently, cluster them into allele groups / haplotypes.
    logObject.info("Extracting and clustering GCF genes with their flanks.")
    bowtie2_reference, instances_to_haplotypes = lsaBGC.extractGeneWithFlanksAndCluster(gene_to_cog,)
    logObject.info("Successfully extracted genes with flanks and clustered them into discrete haplotypes.")

    # Step 4: Align paired-end reads to database genes with surrounding flanks
    bowtie2_dir = outdir + 'Bowtie2_Alignments/'
    if not os.path.isfile(bowtie2_dir): os.system('mkdir %s' % bowtie2_dir)

    bowtie2_cores = cores
    bowtie2_pool_size = 1
    if cores >= 4:
        bowtie2_cores = 4
        bowtie2_pool_size = int(cores / 4)

    bowtie2_args = []
    process_args = []
    with open(sequence_specs_file) as ossf:
        for line in ossf:
            line = line.strip()
            sample, frw_read, rev_read = line.split('\t')
            bowtie2_args.append([sample, frw_read, rev_read, bgc_ref_concat_fasta.split('.fasta')[0], bowtie2_dir, bowtie2_cores])
            process_args.append([sample, bowtie2_dir + sample + '.filtered.sorted.bam', bgc_ref_concat_fasta, cog_gene_to_rep, cog_read_dir])


    # Step 5: Determine haplotypes found in samples and identify supported novelty SNVs
    results_dir = outdir + 'Parsed_Results/'

    # Step 6: Construct summary matrices

def bgmeta(bgc_specs_file, sequence_specs_file, orthofinder_matrix, outdir, cores):
    bgc_ref_concat_fasta = outdir + 'BGCs_Collapsed.fasta'
    bowtie2_dir = outdir + 'Bowtie2_Alignments_BGCs/'
    cog_read_dir = outdir + 'Read_Partitioning_Among_COG_Alleles/'
    #os.system('mkdir %s %s' % (bowtie2_dir, cog_read_dir))

    bowtie2_cores = cores
    bowtie2_pool_size = 1
    if cores >= 4:
        bowtie2_cores = 4
        bowtie2_pool_size = int(cores / 4)

    bowtie2_args = []
    process_args = []
    with open(sequence_specs_file) as ossf:
        for line in ossf:
            line = line.strip()
            sample, frw_read, rev_read = line.split('\t')
            bowtie2_args.append([sample, frw_read, rev_read, bgc_ref_concat_fasta.split('.fasta')[0], bowtie2_dir, bowtie2_cores])
            process_args.append([sample, bowtie2_dir + sample + '.filtered.sorted.bam', bgc_ref_concat_fasta, cog_gene_to_rep, cog_read_dir])

    #p = multiprocessing.Pool(bowtie2_pool_size)
    #p.map(bowtie2Alignment, bowtie2_args)
    #p.close()

    p = multiprocessing.Pool(cores)
    p.map(assessBestAllelicMatches, process_args)
    p.close()



    sys.exit(0)

if __name__ == '__main__':
    lsaBGC_MetaNovelty()