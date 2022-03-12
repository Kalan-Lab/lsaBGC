import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from collections import defaultdict
from scipy import stats, spatial
from operator import itemgetter
import pysam

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: scrutinizeNovelSNVSupport.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--novel_snvs', help="Final Potentially_Novel_SNVs report.", required=True, default=None)
    parser.add_argument('-o', '--output', help="Path to output file - filtered novel SNVs report.", required=True, default=None)
    parser.add_argument('-b', '--comparator_bam_listings', help="File where each line has three columns: (1) sample name, and (2) path to a BAM file.", required=False, default=None)
    args = parser.parse_args()
    return args

def main():
    myargs = create_parser()

    novel_snvs_file = os.path.abspath(myargs.novel_snvs)
    comparator_bam_listings_file = myargs.comparator_bam_listings
    output_file = myargs.output

    try:
        assert(os.path.isfile(novel_snvs_file))
        assert(os.path.isfile(comparator_bam_listings_file))
    except:
        raise RuntimeError('One or more input files do not exist. Exiting now ...')

    # parse and check for read mapping
    reads_perfectly_mapping = defaultdict(set)
    with open(comparator_bam_listings_file) as ocblf:
        for i, line in enumerate(ocblf):
            line = line.strip('\n')
            sample, alignment_bam_file = line.split('\t')

            try:
                bam_file_handle = pysam.AlignmentFile(alignment_bam_file, "rb")
                for read_alignment in bam_file_handle.fetch():
                    read_name = read_alignment.query_name
                    perfect_match = True
                    for b in read_alignment.get_aligned_pairs(with_seq=True):
                        if b[0] == None or b[1] == None:
                            perfect_match = False
                        else:
                            if b[2].islower():
                                perfect_match = False
                    if perfect_match:
                        reads_perfectly_mapping[sample].add(read_name)
            except:
                raise RuntimeError('Difficulty reading in reference BAM alignment file %s on line %d.' % (alignment_bam_file, i))

    output_handle = open(output_file, 'w')
    with open(novel_snvs_file) as onsf:
        for i, line in enumerate(onsf):
            line = line.strip('\n')
            if i == 0:
                output_handle.write(line + '\n')
                continue
            ls = line.split('\t')
            sample = ls[1]
            reads = set([x.strip() for x in ls[-1].split(',')])

            ma_found_flag = False
            for r in reads:
                if r in reads_perfectly_mapping[sample]:
                    ma_found_flag = True
            if ma_found_flag:
                output_handle.write('\t'.join(ls[:-1]) + '\t' + ', '.join(retained_reads) + '\n')
    output_handle.close()

    # Exit program
    sys.exit(0)

if __name__ == '__main__':
    main()
