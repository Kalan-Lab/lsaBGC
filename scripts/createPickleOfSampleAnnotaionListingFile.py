# !/usr/bin/env python

import os
import sys
import _pickle as cPickle
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
	Program: createPickleOfSampleAnnotationListingFile.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-l', '--annotation_listing',
                        help='Path to tab delimited text file for samples with three columns: (1) sample name (2) Prokka generated Genbank file (*.gbk), and (3) Prokka generated predicted-proteome file (*.faa). Please remove troublesome characters in the sample name.',
                        required=True)
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

        pickle_data = {'paths': dict(sample_prokka_data), 'gbk_info': dict(sample_gbk_info)}
        of = open(output_file, 'wb')
        cPickle.dump(pickle_data, of)
        of.close()

    # Exit program
    sys.exit(0)


if __name__ == '__main__':
    main()
