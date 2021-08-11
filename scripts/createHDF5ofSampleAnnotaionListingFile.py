# !/usr/bin/env python

import os
import sys
import h5py
import traceback
from time import sleep
import argparse
import multiprocessing
from collections import defaultdict
from lsaBGC import util
from lsaBGC import processing

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: createHDF5ofSampleAnnotationListingFile.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-l', '--annotation_listing', help='Path to tab delimited text file for samples with three columns: (1) sample name (2) Prokka generated Genbank file (*.gbk), and (3) Prokka generated predicted-proteome file (*.faa). Please remove troublesome characters in the sample name.', required=True)
    parser.add_argument('-o', '--output_file', help="Output file.", required=True)
    parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)
    args = parser.parse_args()

    return args

def main():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """
    myargs = create_parser()

    annotation_listing_file = os.path.abspath(myargs.annotation_listing)
    output_file = os.path.abspath(myargs.output_file)

    ### vet input files quickly
    try:
        assert (os.path.isfile(annotation_listing_file))
    except:
        raise RuntimeError('Input file does not exist. Exiting now ...')

    if os.path.isfile(output_file):
        sys.stderr.write("Output file exists. Overwriting in 5 seconds ...\n ")
        sleep(5)

    """
    PARSE OPTIONAL INPUTS
    """

    cores = myargs.cores

    """
    START WORKFLOW
    """

    sample_prokka_data = processing.readInAnnotationFilesForExpandedSampleSet(annotation_listing_file)

    with multiprocessing.Manager() as manager:
        sample_gbk_info = manager.dict()
        genbanks = []
        for sample in sample_prokka_data:
            sample_genbank = sample_prokka_data[sample]['genbank']
            genbanks.append([sample, sample_genbank, sample_gbk_info])

        with manager.Pool(cores) as pool:
            pool.map(util.parseGenbankAndFindBoundaryGenes, genbanks)

        with h5py.File(output_file, 'w') as of:
            paths_grp = of.create_group("paths")
            gbk_info_grp = of.create_group("gbk_info")
            for sample in sample_gbk_info:
                sample_gbk_path = sample_prokka_data[sample]['genbank']
                sample_prot_path = sample_prokka_data[sample]['predicted_proteome']
                sample_paths_grp = paths_grp.create_group(sample)
                sample_paths_grp["genbank"] = sample_gbk_path
                sample_paths_grp["predicted_proteome"] = sample_prot_path

                gene_to_scaff, scaff_genes, bound_genes, gito, goti = sample_gbk_info[sample]
                sample_gbk_info_grp = gbk_info_grp.create_group(sample)
                sample_gbk_info_grp["gene_to_scaff"] = dict(gene_to_scaff)
                sample_gbk_info_grp["scaff_genes"] = dict(scaff_genes)
                sample_gbk_info_grp["bound_genes"] = set(bound_genes)
                sample_gbk_info_grp["gito"] = dict(gito)
                sample_gbk_info_grp["goti"] = dict(goti)

    # Exit program
    sys.exit(0)

if __name__ == '__main__':
    main()