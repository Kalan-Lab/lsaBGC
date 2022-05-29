import os
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict
from time import sleep
from lsaBGC import util
import subprocess
import traceback
import multiprocessing

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-3])

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: setup_annotation_dbs.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
        
    Downloads KOfam profile HMMs and annotation information. By using this script and downloading the 
    data you must agree to their LICENSE.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-p', '--download_path',
                        help='Path to where database is installed.',
                        required=False, default=lsaBGC_main_directory + 'db/')

    args = parser.parse_args()
    return args


def setup_annot_dbs():
    myargs = create_parser()

    download_path = os.path.abspath(myargs.download_path) + '/'

    try:
        assert(os.path.isdir(download_path))
    except:
        sys.stderr('Error: Provided directory for downloading annotation files does not exist!')

    # download files in requested directory
    try:
        os.chdir(download_path)
        os.system('wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz')
        os.system('wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz')
        os.system('gunzip *.gz')
        ko_annot_info_file = download_path + 'ko_list'
        ko_phmm_file = download_path + 'profile.hmm'
        assert(os.path.isfile(ko_annot_info_file))
        assert(os.path.isdir(download_path + 'profiles/'))

        if os.path.isfile(ko_phmm_file):
            os.system('rm -f %s' % ko_phmm_file)
        for f in os.listdir(download_path + 'profiles/'):
            os.system('cat %s >> %s' % (download_path + 'profiles/' + f, ko_phmm_file))
        assert(os.path.isfile(ko_phmm_file))
        listing_file = lsaBGC_main_directory + 'db/kofam_location_paths.txt'
        listing_handle = open(listing_file, 'w')
        listing_handle.write(ko_annot_info_file + '\t' + ko_phmm_file + '\n')
        linding_handle.close()
    except:
        sys.stderr('Error: issues with downloading or seting up annotation database files! Please post to Github issues if unable to figure out!')

    sys.exit(0)

setup_annot_dbs()