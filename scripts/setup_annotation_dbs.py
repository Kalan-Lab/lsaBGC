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

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

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
        sys.stderr.write('Error: Provided directory for downloading annotation files does not exist!\n')

    # download files in requested directory
    try:
        os.chdir(download_path)

        ko_annot_info_file = download_path + 'ko_list'
        ko_phmm_file = download_path + 'profile.hmm'

        if not os.path.isfile(ko_annot_info_file) or not os.path.isfile(ko_phmm_file):
            os.system('rm -f %s' % ko_annot_info_file)
            # Download KOfam HMMs
            print('Setting up KO database!')
            os.system('wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz')
            os.system('wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz')
            os.system('gunzip ko_list.gz')
            os.system('tar -zxvf profiles.tar.gz')

            assert(os.path.isfile(ko_annot_info_file))
            assert(os.path.isdir(download_path + 'profiles/'))

            if os.path.isfile(ko_phmm_file):
                os.system('rm -f %s' % ko_phmm_file)
            for f in os.listdir(download_path + 'profiles/'):
                os.system('cat %s >> %s' % (download_path + 'profiles/' + f, ko_phmm_file))

            assert(os.path.isfile(ko_phmm_file))
            listing_file = lsaBGC_main_directory + 'db/database_location_paths.txt'
            listing_handle = open(listing_file, 'w')
            listing_handle.write('ko\t' + ko_annot_info_file + '\t' + ko_phmm_file + '\n')
            listing_handle.close()
            os.system('rm -rf %s %s' % (download_path + 'profiles/', download_path + 'profiles.tar.gz'))

        pgap_info_file = download_path + 'hmm_PGAP.tsv'
        pgap_hmm_file = download_path + 'PGAP.hmm'

        if not os.path.isfile(pgap_info_file) or not os.path.isfile(pgap_hmm_file):
            # Download PGAP HMMs (Includes TIGR)
            print('Setting up PGAP database!')
            os.system('wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz')
            os.system('wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv')
            os.system('tar -zxvf hmm_PGAP.HMM.tgz')

            assert(os.path.isfile(pgap_info_file))
            assert(os.path.isdir(download_path + 'hmm_PGAP/'))

            for f in os.listdir(download_path + 'hmm_PGAP/'):
                os.system('cat %s >> %s' % (download_path + 'hmm_PGAP/' + f, pgap_hmm_file))
            assert(os.path.isfile(ko_phmm_file))
            listing_file = lsaBGC_main_directory + 'db/database_location_paths.txt'
            listing_handle = open(listing_file, 'a+')
            listing_handle.write('pgap\t' + pgap_info_file + '\t' + pgap_hmm_file + '\n')
            listing_handle.close()
            os.system('rm -rf %s %s' % (download_path + 'hmm_PGAP/', download_path + 'hmm_PGAP.HMM.tgz'))

    except:
        sys.stderr.write('Error: issues with downloading or seting up annotation database files! Please post to Github issues if unable to figure out!\n')

    sys.exit(0)

setup_annot_dbs()