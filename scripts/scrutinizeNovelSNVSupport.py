import os
import sys
from Bio import SeqIO
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
    parser.add_argument('-f', '--focal_alignment_dir', help="Path to DiscoVary Result output subdirectory Alignment_Results/ for focal GCF of interest.", required=True, default=None)
    parser.add_argument('-r', '--reference_genes', help="Path to GCF_Genes_Representatives.fasta file in focal GCF's DiscoVary Results Directory")
    parser.add_argument('-o', '--output', help="Path to output file - filtered novel SNVs report.", required=True, default=None)
    parser.add_argument('-b', '--comparator_bam_listings', help="File where each line has three columns: (1) sample name, (2) path to a BAM file, and (3) path to a reference FASTA.", required=False, default=None)
    parser.add_argument('-k', '--kraken_listings', help="File where each line has two columns: (1) sample name and (2) path to Kraken2 results file.", required=False, default=None)
    parser.add_argument('-t', '--kraken_taxa', nargs="+", help="List of taxa id in Kraken database to retain. Should be surrounded by quotes to avoid issues.", required=False, default=[])
    parser.add_argument('-m', '--minimum_depth', type=int, help="Minimum number of reads needed for reporting a novel SNV.", required=False, default=5)
    parser.add_argument('-s', '--stringent', action='store_true', help="Perform more stringent assessment, also tosses out reads which align equally-well by score to a comparator BAM (requires more careful thinking about comparator set of BAMs).", required=False, default=False)
    args = parser.parse_args()
    return args

def main():
    myargs = create_parser()

    novel_snvs_file = os.path.abspath(myargs.novel_snvs)
    focal_alignment_dir = os.path.abspath(myargs.focal_alignment_dir) + '/'
    reference_genes_file = os.path.abspath(myargs.reference_genes)
    comparator_bam_listings_file = myargs.comparator_bam_listings
    kraken_listings_file = myargs.kraken_listings
    kraken_taxa = myargs.kraken_taxa
    output_file = myargs.output

    minimum_depth = myargs.minimum_depth
    stringent = myargs.stringent

    output = os.path.abspath(myargs.output)

    try:
        assert(os.path.isfile(novel_snvs_file))
        assert(os.path.isdir(focal_alignment_dir))
        assert(os.path.isfile(reference_genes_file))
        try:
            assert(os.path.isfile(comparator_bam_listings_file) or os.path.isfile(kraken_listings_file))
        except:
            raise RuntimeError("Neither comparator bam listing or kraken result listings provided - thus can't filter.")
    except:
        raise RuntimeError('One or more input files do not exist. Exiting now ...')

    # parse novel_snvs_file for snv supporting reads
    snv_supporting_reads = defaultdict(set)
    with open(novel_snvs_file) as onsf:
        for i, line in enumerate(onsf):
            if i == 0: continue
            line = line.strip('\n')
            ls = line.split('\t')
            sample = ls[1]
            reads = set([x.strip() for x in ls[-1].split(',')])
            ref_genome = ls[-7]
            ref_gene = ls[-6]
            ref_id = ref_genome + '|' + ref_gene
            snv_supporting_reads[sample] = snv_supporting_reads[sample].union(reads)

    ref_gene_lens = {}
    with open(reference_genes_file) as orgf:
        for rec in SeqIO.parse(orgf, 'fasta'):
            ref_gene_lens[rec.id] = len(str(rec.seq))

    focal_read_alignment_scores = defaultdict(list)
    for f in os.listdir(focal_alignment_dir):
        if not f.endswith('_topaligns.sorted.bam'): continue
        sample = f.split('_topaligns.sorted.bam')[0]
        bam_file_handle = pysam.AlignmentFile(focal_alignment_dir + f, "rb")
        for ref_gene in ref_gene_lens:
            reference = '|'.join(ref_gene.split('|')[-2:])
            for read_alignment in bam_file_handle.fetch(ref_gene, 0, ref_gene_lens[ref_gene]):
                read_name = read_alignment.query_name
                read_ascore = float(read_alignment.tags[0][1])
                if read_name in snv_supporting_reads[sample]:
                    focal_read_alignment_scores[read_name].append([sample, read_ascore])

    top_read_reflexive_alignment_scores = defaultdict(lambda: defaultdict(float))
    for read in focal_read_alignment_scores:
        for i, read_info in enumerate(sorted(focal_read_alignment_scores[read], key=itemgetter(1), reverse=True)):
            if i == 0:
                top_read_reflexive_alignment_scores[read_info[0]][read] = read_info[1]

    reads_with_conflicting_support = defaultdict(set)
    if os.path.isfile(comparator_bam_listings_file):
        # parse and check for read mapping
        with open(comparator_bam_listings_file) as ocblf:
            for i, line in enumerate(ocblf):
                line = line.strip('\n')
                sample, alignment_bam_file, reference_fasta_file = line.split('\t')

                ref_seq_lens = defaultdict(int)
                try:
                    with open(reference_fasta_file) as orff:
                        for rec in SeqIO.parse(orff, 'fasta'):
                            ref_seq_lens[rec.id] = len(str(rec.seq))
                except:
                    raise RuntimeError('Difficulty reading in reference FASTA file %s on line %d.' % (reference_fasta_file, i))

                bam_file_handle = None
                try:
                    bam_file_handle = pysam.AlignmentFile(alignment_bam_file, "rb")
                except:
                    raise RuntimeError('Difficulty reading in reference BAM alignment file %s on line %d.' % (alignment_bam_file, i))

                for ref in ref_seq_lens:
                    for read_alignment in bam_file_handle.fetch(ref, 0, ref_seq_lens[ref]):
                        read_name = read_alignment.query_name
                        read_ascore = float(read_alignment.tags[0][1])
                        if read_name in top_read_reflexive_alignment_scores.keys() and \
                            (read_ascore > top_read_reflexive_alignment_scores[sample][read_name] or
                             (read_ascore == top_read_reflexive_alignment_scores[sample][read_name] and stringent)):
                            reads_with_conflicting_support[sample].add(read_name)

    if os.path.isfile(kraken_listings_file) and len(kraken_taxa) == 0:
        kraken_taxa = set(kraken_taxa)
        # parse and check for read mapping
        try:
            with open(kraken_listings_file) as oklf:
                for line in oklf:
                    line = line.strip('\n')
                    sample, kraken_file = line.split('\t')
                    with open(kraken_file) as okf:
                        for l in okf:
                            l = l.strip('\n').split('\t')
                            if not l[1] + '/1' in snv_supporting_reads[sample] and not ls[1] + '/2' in snv_supporting_reads[sample]:
                                continue
                            if not l[2] in kraken_taxa:
                                if l[1] + '/1' in snv_supporting_reads:
                                    reads_with_conflicting_support[sample].add(l[1] + '/1')
                                elif ls[1] + '/2' in snv_supporting_reads:
                                    reads_with_conflicting_support[sample].add(l[1] + '/2')
        except:
            raise RuntimeError('Difficulty reading Kraken2 file %s.' % kraken_file)

    output_handle = open(output_file, 'w')
    with open(novel_snvs_file) as onsf:
        for i, line in enumerate(onsf):
            line = line.strip('\n')
            if i == 0:
                output_handle.write(line + '\n')
                continue
            ls = line.split('\t')
            reads = set([x.strip() for x in ls[-1].split(',')])
            retained_reads = [r for r in reads if not r in reads_with_conflicting_support]
            if len(retained_reads) >= minimum_depth:
                output_handle.write('\t'.join(ls[:-1]) + '\t' + ', '.join(retained_reads) + '\n')
    output_handle.close()

    # Exit program
    sys.exit(0)

if __name__ == '__main__':
    main()
