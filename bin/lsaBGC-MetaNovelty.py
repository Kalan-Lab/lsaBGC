"""
Author: Rauf Salamzade
Kalan Lab - UW Madison - MMI
01/08/2021
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from operator import itemgetter
import argparse
import multiprocessing
import statistics
import pysam
import copy

# define global variables
comp_gene_info = {}
bgc_cog_genes = {}
CORES = 1

def read_pair_generator(bam, region_string=None, start=None, stop=None):
    """
    Function taken from: https://www.biostars.org/p/306041/
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region_string, start=start, stop=stop):
        if not read.is_proper_pair or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def defineKnownAlleleGroups(input):
    cog, nucl_seq_dir, nucl_cdhit_dir = input
    cog_nucl_fasta = nucl_seq_dir + '/' + cog + '.fasta'
    cog_cdhit_result = nucl_cdhit_dir + '/' + cog + '.fasta'

    cog_nucl_handle = open(cog_nucl_fasta, 'w')
    for gene in bgc_cog_genes[cog]:
        cog_nucl_handle.write('>' + comp_gene_info[gene]['bgc_name'] + '|' + gene + '\n' + comp_gene_info[gene]['nucl_seq'] + '\n')
    cog_nucl_handle.close()
    os.system("cd-hit-est -i %s -o %s -G 1 -g 1 -d 0 -n 10 -M 2000 -c 0.95 -aL 0.95 -aS 0.95 -T 1" % (cog_nucl_fasta, cog_cdhit_result))

def bowtie2Alignment(input_args):
    sample, frw_read, rev_read, bgc_cdhit, bowtie2_dir, cores = input_args
    try:
        sam_file = bowtie2_dir + sample + '.sam'
        bam_file = bowtie2_dir + sample + '.bam'
        bam_file_sorted = bowtie2_dir + sample + '.sorted.bam'
        bam_file_filtered = bowtie2_dir + sample + '.filtered.bam'
        bam_file_filtered_sorted = bowtie2_dir + sample + '.filtered.sorted.bam'

        os.system("bowtie2 --very-sensitive --no-mixed --no-discordant --no-unal -a -x %s -1 %s -2 %s -S %s -p %d" % (bgc_cdhit, frw_read, rev_read, sam_file, cores))
        os.system("samtools view -h -Sb %s > %s" % (sam_file, bam_file))
        os.system("samtools sort -@ %d %s -o %s" % (cores, bam_file, bam_file_sorted))
        os.system("samtools index %s" % bam_file_sorted)

        bam_handle = pysam.AlignmentFile(bam_file_sorted, 'rb')
        filt_bam_handle = pysam.AlignmentFile(bam_file_filtered, "wb", template=bam_handle)

        for read in bam_handle.fetch():
            if not read.is_proper_pair or read.is_supplementary: continue
            filt_bam_handle.write(read)

        filt_bam_handle.close()
        bam_handle.close()
        os.system("samtools sort -@ %d %s -o %s" % (cores, bam_file_filtered, bam_file_filtered_sorted))
        os.system("samtools index %s" % bam_file_filtered_sorted)
        os.system("rm -f %s %s %s %s %s" % (sam_file, bam_file, bam_file_sorted, bam_file_filtered, bam_file_sorted + '.bai'))
    except:
        os.system('rm -f %s/%s*' % (bowtie2_dir, sample))

def assessBestAllelicMatches(input_args):
    try:
        sample, bam_alignment, ref_fasta, cog_gene_to_rep, res_dir = input_args



        cog_rep_genes = defaultdict(set)
        for g, r in cog_gene_to_rep.items():
            cog_rep_genes[r].add(g)

        if not os.path.isfile(bam_alignment): return
        result_file = res_dir + sample + '.txt'
        outf = open(result_file, 'w')
        outf.write('\t'.join(['# cog', 'allele_representative', 'reads', 'reads_with_novelty', 'reads_uniquely_mapping', 'reads_uniquely_mapping_with_novelty']) + '\n')

        bam_handle = pysam.AlignmentFile(bam_alignment, 'rb')

        topaligns_file = res_dir + sample + '_topaligns.bam'
        topaligns_file_sorted = res_dir + sample + '_topaligns.sorted.bam'
        topaligns_handle = pysam.AlignmentFile(topaligns_file, "wb", template=bam_handle)

        unialigns_file = res_dir + sample + '_unialigns.bam'
        unialigns_file_sorted = res_dir + sample + '_unialigns.sorted.bam'
        unialigns_handle = pysam.AlignmentFile(unialigns_file, "wb", template=bam_handle)

        for cog, cog_genes in bgc_cog_genes.items():
            read_ascores_per_allele = defaultdict(list)
            read_genes_mapped = defaultdict(set)
            snv_counts = defaultdict(int)
            cog_genes_covered = 0
            rep_alignments = defaultdict(lambda: defaultdict(set))
            with open(ref_fasta) as opff:
                for rec in SeqIO.parse(opff, 'fasta'):
                    cog_genes_for_bgc = set([x for x in cog_genes if comp_gene_info[x]['bgc_name'] == rec.id])
                    for i, g in enumerate(cog_genes_for_bgc):
                        ginfo = comp_gene_info[g]
                        gstart = ginfo['start']
                        gend = ginfo['end']

                        gene_length = gend-gstart+1
                        gene_covered_1 = 0
                        gene_covered_3 = 0
                        for pileupcolumn in bam_handle.pileup(contig=rec.id, start=gstart, stop=gend+1, stepper="nofilter", truncate=True):
                            pos_depth = 0
                            for pileupread in pileupcolumn.pileups:
                                read = pileupread.alignment
                                if pileupread.is_del or pileupread.is_refskip or not read.is_proper_pair: continue
                                if read.query_qualities[pileupread.query_position] < 20: continue
                                pos_depth += 1
                            if pos_depth >= 1:
                                gene_covered_1 += 1
                                if pos_depth >= 3:
                                    gene_covered_3 += 1

                        for read1_alignment, read2_alignment in read_pair_generator(bam_handle, region_string=rec.id, start=gstart, stop=gend):
                            read1_ascore = read1_alignment.tags[0][1]
                            read2_ascore = read2_alignment.tags[0][1]
                            combined_ascore = read1_ascore + read2_ascore

                            snvs = set([])
                            g_rep = cog_gene_to_rep[rec.id + '|' + g]
                            for b in read1_alignment.get_aligned_pairs(with_seq=True):
                                if b[-1] is not None and b[-1].islower():
                                    ref_pos, alt_al = b[1:]
                                    alt_al = alt_al.upper()
                                    snvs.add(str(rec.id) + '_|_' + str(ref_pos) + '_|_' + alt_al)
                                    snv_counts[str(rec.id) + '_|_' + str(ref_pos) + '_|_' + alt_al] += 1
                            for b in read2_alignment.get_aligned_pairs(with_seq=True):
                                if b[-1] is not None and b[-1].islower():
                                    ref_pos, alt_al = b[1:]
                                    alt_al = alt_al.upper()
                                    snvs.add(str(rec.id) + '_|_' + str(ref_pos) + '_|_' + alt_al)
                                    snv_counts[str(rec.id) + '_|_' + str(ref_pos) + '_|_' + alt_al] += 1
                            read_genes_mapped[read1_alignment.query_name].add(rec.id + '|' + g)
                            r1a = copy.deepcopy(read1_alignment)
                            r2a = copy.deepcopy(read2_alignment)
                            rep_alignments[rec.id + '|' + g][read1_alignment.query_name].add(tuple([r1a, r2a]))
                            read_ascores_per_allele[read1_alignment.query_name].append([g_rep.split('|')[1], g_rep.split('|')[0], combined_ascore, snvs, g])

                        gene_coverage_1 = gene_covered_1/float(gene_length)
                        gene_coverage_3 = gene_covered_3/float(gene_length)
                        if gene_coverage_1 < 0.90: continue
                        cog_genes_covered += 1

            if cog_genes_covered/float(len(cog_genes)) < 0.80: continue

            allele_reads = defaultdict(set)
            allele_reads_with_mismatch = defaultdict(set)
            multi_partitioned_reads = set([])
            for read in read_ascores_per_allele:
                top_score = -1000000
                top_score_grep = None
                for i, align in enumerate(sorted(read_ascores_per_allele[read], key=itemgetter(2), reverse=True)):
                    g_rep = align[1] + '|' + align[0]
                    if i == 0: top_score = align[2]; top_score_grep = g_rep
                    if align[2] == top_score:
                        g_map = read_genes_mapped[read].intersection(cog_rep_genes[g_rep])
                        g_map_prop = len(g_map)/float(len(cog_rep_genes[g_rep]))
                        if g_map_prop < 0.80: continue

                        allele_reads[g_rep].add(read)
                        for snv in align[3]:
                            if snv_counts[snv] >= 2:
                                allele_reads_with_mismatch[g_rep].add(read)

                    if g_rep != top_score_grep and align[2] == top_score and i > 0:
                        multi_partitioned_reads.add(read)

            for al in allele_reads:
                for r in allele_reads[al]:
                    for pa in rep_alignments[al][r]:
                        topaligns_handle.write(pa[0])
                        topaligns_handle.write(pa[1])
                        if not r in multi_partitioned_reads:
                            unialigns_handle.write(pa[0])
                            unialigns_handle.write(pa[1])
                outf.write('\t'.join([str(x) for x in [cog, al, len(allele_reads[al]), len(allele_reads_with_mismatch[al]),
                                                       len(allele_reads[al].difference(multi_partitioned_reads)),
                                                       len(allele_reads_with_mismatch[al].difference(multi_partitioned_reads))]]) + '\n')


        outf.close()
        topaligns_handle.close()
        unialigns_handle.close()
        bam_handle.close()

        os.system("samtools sort -@ %d %s -o %s" % (1, unialigns_file, unialigns_file_sorted))
        os.system("samtools index %s" % unialigns_file_sorted)

        os.system("samtools sort -@ %d %s -o %s" % (1, topaligns_file, topaligns_file_sorted))
        os.system("samtools index %s" % topaligns_file_sorted)
    except:
        pass

def parseGenbanksForGeneInfo(gbk, bgc_name):
    """
    :param gbk: AntiSMASH Genbank file
    :return: list of lists (genes) where each sub-list corresponding to a gene contains:
        1. locus tag
        2. start (smaller) coordinate
        3. end (larger) coordinate
        4. gene direction
        5. protein sequence
    each gene is presented in the same order it appears in the AntiSMASH genbank file
    """
    bgc_info = []
    domains = []
    core_positions = set([])
    full_sequence = ""
    with open(gbk) as ogbk:
        domain_feature_types = ['PFAM_domain', 'CDS_motif', 'aSDomain']
        for rec in SeqIO.parse(ogbk, 'genbank'):
            for feature in rec.features:
                if feature.type in domain_feature_types:
                    start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    aSDomain = "NA"
                    description = "NA"
                    try: aSDomain = feature.qualifiers.get('aSDomain')[0]
                    except: pass
                    try: description = feature.qualifiers.get('description')[0]
                    except: pass
                    domains.append({'start': start, 'end': end, 'type': feature.type, 'aSDomain': aSDomain, 'description': description})
                elif feature.type == 'protocluster':
                    detection_rule = feature.qualifiers.get('detection_rule')[0]
                    product = feature.qualifiers.get('product')[0]
                    contig_edge = feature.qualifiers.get('contig_edge')[0]
                    full_sequence = str(rec.seq)
                    bgc_info.append({'detection_rule': detection_rule, 'product': product, 'contig_edge': contig_edge, 'full_sequence': str(rec.seq)})
                elif feature.type == 'proto_core':
                    core_start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    core_end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    core_positions = core_positions.union(set(range(core_start, core_end+1)))

    if len(bgc_info) == 0:
        bgc_info = [{'detection_rule': 'NA', 'product': 'NA', 'contig_edge': 'NA', 'full_sequence': full_sequence}]

    genes = {}
    with open(gbk) as ogbk:
        for rec in SeqIO.parse(ogbk, 'genbank'):
            for feature in rec.features:
                if feature.type == "CDS":
                    lt = feature.qualifiers.get('locus_tag')[0]
                    prot_seq = feature.qualifiers.get('translation')[0]
                    start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    direction = str(feature.location).split('(')[1].split(')')[0]
                    product = feature.qualifiers.get('product')[0]
                    grange = set(range(start, end+1))
                    gene_domains = []
                    for d in domains:
                        drange = set(range(d['start'], d['end']+1))
                        if len(drange.intersection(grange)) > 0:
                            gene_domains.append(d)
                    nucl_seq = full_sequence[start:end]
                    core_overlap = False
                    if len(grange.intersection(core_positions)) > 0: core_overlap = True
                    if direction == '-':
                        nucl_seq = str(Seq(full_sequence[start:end]).reverse_complement())
                    genes[lt] = {'bgc_name': bgc_name, 'start': start, 'end': end, 'direction': direction,
                                 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
                                 'gene_domains': gene_domains, 'core_overlap': core_overlap}
    return([genes, bgc_info])

def parseOrthoFinderMatrix(orthofinder_matrix_file, all_gene_lts):
    """
    :param orthofinder_matrix: OrthoFinderV2 matrix Orthogroups.csv file, should also include singleton orthologroups
    :param all_gene_lts: Set of all the relevant gene locus tag identifiers found in BGC Genbanks
    :return: dictionary mapping gene locus tags to CoGs
    """
    gene_to_cog = {}
    cog_genes = defaultdict(set)
    cog_gene_counts = defaultdict(lambda: 'NA')
    with open(orthofinder_matrix_file) as ofm:
        for i, line in enumerate(ofm):
            if i == 0: continue
            line = line.strip()
            ls = line.split('\t')
            cog = ls[0]
            flag = False
            for sgs in ls[1:]:
                for g in sgs.split(', '):
                    if g in all_gene_lts:
                        flag = True
                        gene_to_cog[g] = cog
                        cog_genes[cog].add(g)
            if flag:
                gene_counts = []
                for sgs in ls[1:]:
                    gene_counts.append(len(sgs.split(', ')))
                cog_gene_counts[cog] = statistics.median(gene_counts)

    return([gene_to_cog, cog_genes, cog_gene_counts])

def bgmeta(bgc_specs_file, sequence_specs_file, orthofinder_matrix, outdir, cores):

    global CORES
    CORES = cores
    ### vet input files quickly
    sys.stderr.write('Checking if input files exist and are in valid formatting...\n')
    try:
        assert(os.path.isfile(sequence_specs_file))
        assert(os.path.isfile(orthofinder_matrix))
        assert(os.path.isfile(bgc_specs_file))
        with open(bgc_specs_file) as obsf:
            for line in obsf:
                line = line.strip()
                gbk, sample = line.split('\t')
                assert(os.path.isfile(gbk))
                with open(gbk) as ogbk:
                    for rec in SeqIO.parse(ogbk, 'genbank'):
                        pass
        with open(sequence_specs_file) as ossf:
            for line in ossf:
                line = line.strip()
                sample, frw_read, rev_read = line.split('\t')
                assert(os.path.isfile(frw_read) and os.path.isfile(rev_read))
    except:
        sys.stderr.write('Error validating input files. Exiting now...\n')
        sys.exit(1)
    sys.stderr.write('Input files check out, at least for now...\n')

    outdir = os.path.abspath(outdir) + '/'
    if not os.path.isdir(outdir):
        os.system('mkdir %s' % outdir)

    sample_index = defaultdict(int)
    edited_sample_names = []
    sample_naming = {}
    antismash_genbanks = []
    with open(bgc_specs_file) as obsf:
        for line in obsf:
            line = line.strip()
            gbk, sample = line.split('\t')
            try:
                assert (os.path.isfile(gbk))
                antismash_genbanks.append(gbk)
                if sample_index[sample] == 0:
                    edited_sample_names.append(sample)
                    sample_naming[sample] = sample
                else:
                    edited_sample_names.append(sample + '_' + str(sample_index[sample]+1))
                    sample_naming[sample + '_' + str(sample_index[sample]+1)] = sample
                sample_index[sample] += 1
            except:
                pass

    # parse and process genbanks
    global comp_gene_info, bgc_cog_genes

    bgc_genes = {}
    bgc_information = {}
    all_gene_lts = set([])
    bgc_ref_fasta = outdir + 'All_BGCs.fasta'
    bgc_ref_handle = open(bgc_ref_fasta, 'w')
    for i, gbk in enumerate(antismash_genbanks):
        bgc_name = edited_sample_names[i]
        gene_info, bgc_info = parseGenbanksForGeneInfo(gbk, bgc_name)
        comp_gene_info.update(gene_info)
        bgc_genes[bgc_name] = gene_info.keys()
        bgc_information[bgc_name] = bgc_info
        all_gene_lts = all_gene_lts.union(set(gene_info.keys()))
        bgc_ref_handle.write('>' + bgc_name + '\n' + bgc_info[0]['full_sequence'] + '\n')
    bgc_ref_handle.close()

    # parse orthofinder matrix
    gene_to_cog, bgc_cog_genes, cog_median_counts = parseOrthoFinderMatrix(orthofinder_matrix, all_gene_lts)

    nucl_seq_dir = os.path.abspath(outdir + 'Nucleotide_Sequences') + '/'
    nucl_cdhit_dir = os.path.abspath(outdir + 'Nucleotide_CD-HIT_Sequences') + '/'

    if not os.path.isdir(nucl_seq_dir): os.system('mkdir %s' % nucl_seq_dir)
    if not os.path.isdir(nucl_cdhit_dir): os.system('mkdir %s' % nucl_cdhit_dir)

    input_args = [[x, nucl_seq_dir, nucl_cdhit_dir] for x in bgc_cog_genes.keys()]
    p = multiprocessing.Pool(cores)
    p.map(defineKnownAlleleGroups, input_args)
    p.close()

    cog_gene_to_rep = {}
    for file in os.listdir(nucl_cdhit_dir):
        cog = file.split('.fasta')[0]
        if file.endswith('.fasta'): continue
        cluster = []
        rep = None
        with open(nucl_cdhit_dir + file) as off:
            for line in off:
                line = line.strip()
                ls = line.split()
                if line.startswith('>'):
                    if len(cluster) > 0:
                        for g in cluster:
                            cog_gene_to_rep[g] = rep
                    cluster = []
                    rep = None
                else:
                    gene_id = ls[2][1:-3]
                    cluster.append(gene_id)
                    if line.endswith('*'): rep = gene_id
        if len(cluster) > 0:
            if len(cluster) > 0:
                for g in cluster:
                    cog_gene_to_rep[g] = rep

    bgc_ref_concat_fasta = outdir + 'BGCs_Collapsed.fasta'
    #os.system("cd-hit-est -i %s -o %s -G 1 -g 1 -d 0 -n 10 -M 25000 -c 1.0 -aL 0.0 -aS 1.0 -T %d" % (bgc_ref_fasta, bgc_ref_concat_fasta, cores))
    #os.system("bowtie2-build %s %s" % (bgc_ref_concat_fasta, bgc_ref_concat_fasta.split('.fasta')[0]))

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

    samples = set([])
    cog_allele_representatives = set([])
    sample_allele_reads = defaultdict(lambda: defaultdict(int))
    sample_allele_unique_reads = defaultdict(lambda: defaultdict(int))
    sample_allele_novelty_reads = defaultdict(lambda: defaultdict(int))
    sample_allele_unique_novelty_reads = defaultdict(lambda: defaultdict(int))
    with open(sequence_specs_file) as ossf:
        for line in ossf:
            sample = line.strip().split('\t')[0]
            result_file = cog_read_dir + sample + '.txt'
            samples.add(sample)
            if not os.path.isfile(result_file): continue
            with open(result_file) as orf:
                for i, cog_al in enumerate(orf):
                    if i == 0: continue
                    cog_al = cog_al.strip()
                    cog, allele_representative, reads, reads_with_novelty, reads_uniquely_mapping, reads_uniquely_mapping_with_novelty = cog_al.split('\t')
                    car = cog + '|' + allele_representative
                    cog_allele_representatives.add(car)
                    sample_allele_reads[sample][car] = int(reads)
                    sample_allele_unique_reads[sample][car] = int(reads_uniquely_mapping)
                    sample_allele_novelty_reads[sample][car] = int(reads_with_novelty)
                    sample_allele_unique_novelty_reads[sample][car] = int(reads_uniquely_mapping_with_novelty)

    final_matrix_reads_file = outdir + 'Sample_by_OG_Allele_Read_Counts.matrix.txt'
    final_matrix_novelty_reads_file = outdir + 'Sample_by_OG_Allele_Novelty_Read_Counts.matrix.txt'
    final_matrix_unique_reads_file = outdir + 'Final_OG_Allele_Unique_Read_Counts.matrix.txt'
    final_matrix_unique_and_novelty_reads_file = outdir + 'Final_OG_Allele_Unique_and_Novelty_Read_Counts.matrix.txt'

    final_matrix_reads_handle = open(final_matrix_reads_file, 'w')
    final_matrix_novelty_reads_handle = open(final_matrix_novelty_reads_file, 'w')
    final_matrix_unique_reads_handle = open(final_matrix_unique_reads_file, 'w')
    final_matrix_unique_and_novelty_reads_handle = open(final_matrix_unique_and_novelty_reads_file, 'w')

    final_matrix_reads_handle.write('\t'.join(['Sample/OG_Allele'] + list(sorted(cog_allele_representatives))) + '\n')
    final_matrix_novelty_reads_handle.write('\t'.join(['Sample/OG_Allele'] + list(sorted(cog_allele_representatives))) + '\n')
    final_matrix_unique_reads_handle.write('\t'.join(['Sample/OG_Allele'] + list(sorted(cog_allele_representatives))) + '\n')
    final_matrix_unique_and_novelty_reads_handle.write('\t'.join(['Sample/OG_Allele'] + list(sorted(cog_allele_representatives))) + '\n')

    for s in samples:
        printlist_all = [s]
        printlist_uni = [s]
        printlist_nov = [s]
        printlist_uni_nov = [s]
        for c in cog_allele_representatives:
            printlist_all.append(str(sample_allele_reads[s][c]))
            printlist_uni.append(str(sample_allele_unique_reads[s][c]))
            printlist_nov.append(str(sample_allele_novelty_reads[s][c]))
            printlist_uni_nov.append(str(sample_allele_unique_novelty_reads[s][c]))
        final_matrix_reads_handle.write('\t'.join(printlist_all) + '\n')
        final_matrix_novelty_reads_handle.write('\t'.join(printlist_nov) + '\n')
        final_matrix_unique_reads_handle.write('\t'.join(printlist_uni) + '\n')
        final_matrix_unique_and_novelty_reads_handle.write('\t'.join(printlist_uni_nov) + '\n')

    final_matrix_reads_handle.close()
    final_matrix_novelty_reads_handle.close()
    final_matrix_unique_reads_handle.close()

    sys.exit(0)

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Program to explore lineage specific GCF instances across metagenomes.
    """)

    parser.add_argument('-b', '--bgc_specs_file', help='BGC specifications file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
    parser.add_argument('-i', '--sequence_specs_file', help="Sequencing data specifications file. Tab delimited: 1st column contains metagenomic sample name, whereas 2nd and 3rd columns contain full paths to forward and reverse reads, respectively.", required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
    parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)
    args = parser.parse_args()

    bgmeta(args.bgc_specs_file, args.sequence_specs_file, args.orthofinder_matrix, args.output_directory, args.cores)