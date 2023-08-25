import os
import sys
import argparse
from Bio import SeqIO

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
    parser.add_argument('-nk', '--no_ko', action='store_true',
                        help='Do not download KO database.',
                        required=False, default=False)
    parser.add_argument('-dsh', '--download_scc_hmms', action='store_true',
                        help='Download SCC HMMs from GToTree for Actinobacteria, Bacteria, and Universal models. Added for Docker.',
                        required=False, default=False)

    args = parser.parse_args()
    return args

def setup_annot_dbs():
    #TODO Create substructe in db subdirectory to allow force deletion of all previous database download attempts.

    myargs = create_parser()
    download_path = os.path.abspath(myargs.download_path) + '/'
    no_ko = myargs.no_ko
    download_scc_hmms = myargs.download_scc_hmms

    try:
        assert(os.path.isdir(download_path))
    except:
        sys.stderr.write('Error: Provided directory for downloading annotation files does not exist!\n')

    listing_file = lsaBGC_main_directory + 'db/database_location_paths.txt'

    # download files in requested directory
    try:
        os.chdir(download_path)

        ko_annot_info_file = download_path + 'ko_list'
        ko_phmm_file = download_path + 'profile.hmm'

        if (not os.path.isfile(ko_annot_info_file) or not os.path.isfile(ko_phmm_file)) and (not no_ko):
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
            listing_handle = open(listing_file, 'w')
            listing_handle.write('ko\t' + ko_annot_info_file + '\t' + ko_phmm_file + '\n')
            listing_handle.close()
            os.system('rm -rf %s %s' % (download_path + 'profiles/', download_path + 'profiles.tar.gz'))

        pgap_info_file = download_path + 'hmm_PGAP.tsv'
        pgap_phmm_file = download_path + 'PGAP.hmm'

        if not os.path.isfile(pgap_info_file) or not os.path.isfile(pgap_hmm_file):
            # Download PGAP HMMs (Includes TIGR)
            print('Setting up PGAP database!')
            os.system('wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz')
            os.system('wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv')

            extract_hmm_dir = 'hmm_PGAP.HMM/'
            os.mkdir(extract_hmm_dir)
            os.system(' '.join(['tar', '-zxf', 'hmm_PGAP.HMM.tgz', '-C', 'hmm_PGAP.HMM/']))
            assert (os.path.isfile(pgap_info_file))
            assert (os.path.isdir(download_path + 'hmm_PGAP.HMM/'))
            for folder, subs, files in os.walk(extract_hmm_dir):
                for filename in files:
                    if filename.endswith('.HMM') or filename.endswith('.hmm'):
                        hmm_file_path = os.path.abspath(folder + '/' + filename)
                        os.system(' '.join(['cat', hmm_file_path, '>>', pgap_phmm_file]))
            assert(os.path.isfile(pgap_phmm_file))
            listing_handle = open(listing_file, 'a+')
            listing_handle.write('pgap\t' + pgap_info_file + '\t' + pgap_phmm_file + '\n')
            listing_handle.close()
            os.system('rm -rf %s %s' % (download_path + 'hmm_PGAP.HMM/', download_path + 'hmm_PGAP.HMM.tgz'))

        mibig_faa_file = download_path + 'mibig_prot_seqs_3.1.fasta'
        mibig_dmnd_file = download_path + 'mibig.dmnd'
        mibig_info_file = download_path + 'mibig_info.txt'

        if not os.path.isfile(mibig_info_file) or not os.path.isfile(mibig_dmnd_file):
            # Download MIBiG
            print('Setting up MIBiGv3.1 database!')
            os.system('wget https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_3.1.fasta')
            assert(os.path.isfile(mibig_faa_file))
            os.system(' '.join(['diamond', 'makedb', '--in', mibig_faa_file, '-d', mibig_dmnd_file]))
            assert(os.path.isfile(mibig_dmnd_file))
            mdf_handle = open(mibig_info_file, 'w')
            with open(mibig_faa_file) as omf:
                for rec in SeqIO.parse(omf, 'fasta'):
                    mdf_handle.write(rec.id + '\t' + rec.description + '\n')
            mdf_handle.close()
            assert(os.path.isfile(mibig_info_file))

            listing_handle = open(listing_file, 'a+')
            listing_handle.write('mibig\t' + mibig_info_file + '\t' + mibig_dmnd_file + '\n')
            listing_handle.close()

        actino_hmm_file = download_path + 'Actinobacteria.hmm'
        bacteria_hmm_file = download_path + 'Bacteria.hmm'
        universal_hmm_file = download_path + 'Universal_Hug_et_al.hmm'

        if download_scc_hmms and not (os.path.isfile(actino_hmm_file) and os.path.isfile(bacteria_hmm_file) and os.path.isfile(universal_hmm_file)):
            print('Setting up GToTree SCC HMMs!')

            if not os.path.isfile(actino_hmm_file):
                os.system('wget https://zenodo.org/record/7860735/files/Actinobacteria.hmm?download=1 -O ' + actino_hmm_file)
                assert (os.path.isfile(actino_hmm_file))

            if not os.path.isfile(bacteria_hmm_file):
                os.system('wget https://zenodo.org/record/7860735/files/Bacteria.hmm?download=1 -O ' + bacteria_hmm_file)
                assert (os.path.isfile(bacteria_hmm_file))

            if not os.path.isfile(universal_hmm_file):
                os.system('wget https://zenodo.org/record/7860735/files/Universal_Hug_et_al.hmm?download=1 -O ' + universal_hmm_file)
                assert (os.path.isfile(universal_hmm_file))

    except:
        sys.stderr.write('Error: issues with downloading or seting up annotation database files! Please post to Github issues if unable to figure out!\n')

    sys.exit(0)

setup_annot_dbs()
