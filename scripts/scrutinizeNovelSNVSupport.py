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
	
	Program to further scrutunize support for novel SNVs found by lsaBGC-DiscoVary through using results from
	additional alignment (e.g. to whole-genome database) mappings and Kraken2 taxonomic read classification. 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--novel_snvs', help="Final Potentially_Novel_SNVs report.", required=True, default=None)
    parser.add_argument('-f', '--focal_alignment_dir', help="Path to DiscoVary Result output subdirectory Alignment_Results/ for focal GCF of interest.", required=True, default=None)
    parser.add_argument('-r', '--reference_genes', help="Path to GCF_Genes_Representatives.fasta file in focal GCF's DiscoVary Results Directory")
    parser.add_argument('-o', '--output', help="Path to output file - filtered novel SNVs report.", required=True, default=None)
    parser.add_argument('-b', '--comparator_bam_listings', help="File where each line has three columns: (1) sample name, and (2) path to a BAM file.", required=False, default=None)
    parser.add_argument('-k', '--kraken_listings', help="File where each line has two columns: (1) sample name and (2) path to Kraken2 results file.", required=False, default=None)
    parser.add_argument('-t', '--kraken_taxa', nargs="+", help="List of taxa id in Kraken database to retain. Should be surrounded by quotes to avoid issues.", required=False, default=[])
    parser.add_argument('-l', '--max_read_length', help="Max read length.", default=150, required=False)
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
    max_read_length = myargs.max_read_length

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
            #ref_id = ref_genome + '|' + ref_gene
            snv_supporting_reads[sample] = snv_supporting_reads[sample].union(reads)

    ref_gene_lens = {}
    ref_gene_seqs = {}
    with open(reference_genes_file) as orgf:
        for rec in SeqIO.parse(orgf, 'fasta'):
            ref_gene_lens[rec.id] = len(str(rec.seq))
            ref_gene_seqs[rec.id] = str(rec.seq)

    focal_read_alignment_scores = defaultdict(list)
    for f in os.listdir(focal_alignment_dir):
        if not f.endswith('_topaligns.sorted.bam'): continue
        sample = f.split('_topaligns.sorted.bam')[0]
        bam_file_handle = pysam.AlignmentFile(focal_alignment_dir + f, "rb")
        for ref_gene in ref_gene_lens:
            for read_alignment in bam_file_handle.fetch(ref_gene, 0, ref_gene_lens[ref_gene]):
                read_name = read_alignment.query_name
                read_ascore = float(read_alignment.tags[0][1])
                if read_name in snv_supporting_reads[sample]:

                    ref_positions = set([])
                    for b in read_alignment.get_aligned_pairs(with_seq=False):
                        if b[1] == None or b[0] == None: continue
                        ref_pos = b[1] + 1
                        ref_positions.add(ref_pos)

                    min_ref_pos = min(ref_positions)
                    max_ref_pos = max(ref_positions)
                    on_reference_edge = False
                    ref_seq = None
                    if min_ref_pos == 1 or max_ref_pos == ref_gene_lens[ref_gene]:
                            on_reference_edge = True
                            if max_ref_pos == ref_gene_lens[ref_gene]:
                                ref_seq = ref_gene_seqs[ref_gene][min_ref_pos-1:]
                            else:
                                ref_seq = ref_gene_seqs[ref_gene][min_ref_pos-1:max_ref_pos]
                    focal_read_alignment_scores[read_name].append([sample, read_ascore, on_reference_edge, ref_seq])

    top_read_reflexive_alignment_scores = defaultdict(lambda: defaultdict(float))
    for read in focal_read_alignment_scores:
        for i, read_info in enumerate(sorted(focal_read_alignment_scores[read], key=itemgetter(1), reverse=True)):
            if i == 0:
                top_read_reflexive_alignment_scores[read_info[0]][read] = [read_info[1], read_info[2], read_info[3]]

    reads_with_conflicting_support = defaultdict(set)
    if os.path.isfile(comparator_bam_listings_file):
        # parse and check for read mapping
        with open(comparator_bam_listings_file) as ocblf:
            for i, line in enumerate(ocblf):
                line = line.strip('\n')
                sample, alignment_bam_file = line.split('\t')

                try:
                    bam_file_handle = pysam.AlignmentFile(alignment_bam_file, "rb")

                    top_alignments = defaultdict(list)
                    for read_alignment in bam_file_handle.fetch():
                        read_name = read_alignment.query_name
                        read_ascore = float(read_alignment.tags[0][1])
                        if sample in top_read_reflexive_alignment_scores.keys() and \
                            read_name in top_read_reflexive_alignment_scores[sample] and \
                            (read_ascore > top_read_reflexive_alignment_scores[sample][read_name][0] or
                             (read_ascore == top_read_reflexive_alignment_scores[sample][read_name][0] and stringent)):

                            if top_read_reflexive_alignment_scores[sample][read_name][1]:
                                ref_seq = ""
                                for b in read_alignment.get_aligned_pairs(with_seq=True):
                                    if b[2]: ref_seq += b[2].upper()
                                top_alignments[read_name].append([read_ascore, ref_seq])
                            else:
                                reads_with_conflicting_support[sample].add(read_name)

                    for r in top_alignments:
                        for i, ta in enumerate(sorted(top_alignments[r], key=itemgetter(0), reverse=True)):
                            if i == 0:
                                bgc_ref_seq = top_read_reflexive_alignment_scores[sample][r][2].upper()
                                bgc_ref_seq_rc = str(Seq(bgc_ref_seq).reverse_complement())
                                if not (bgc_ref_seq in ta[1] or bgc_ref_seq_rc in ta[1]):
                                    reads_with_conflicting_support[sample].add(read_name)
                except:
                    raise RuntimeError('Difficulty reading in reference BAM alignment file %s on line %d.' % (alignment_bam_file, i))

    if kraken_listings_file and os.path.isfile(kraken_listings_file) and len(kraken_taxa) > 0:
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
                                if l[1] + '/1' in snv_supporting_reads[sample]:
                                    reads_with_conflicting_support[sample].add(l[1] + '/1')
                                elif ls[1] + '/2' in snv_supporting_reads[sample]:
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
            sample = ls[1]
            reads = set([x.strip() for x in ls[-1].split(',')])
            retained_reads = [r for r in reads if not r in reads_with_conflicting_support[sample]]
            if len(retained_reads) >= minimum_depth:
                output_handle.write('\t'.join(ls[:-2]) + '\t' + str(len(retained_reads)) + '\t' + ', '.join(retained_reads) + '\n')
    output_handle.close()

    # Exit program
    sys.exit(0)

if __name__ == '__main__':
    main()
