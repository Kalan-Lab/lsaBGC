"""
Author: Rauf Salamzade
Kalan Lab - UW Madison, MMI
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import argparse
import random
import egglib
from operator import itemgetter
from ete3 import Tree
import multiprocessing
import statistics
from scipy.stats import fisher_exact

# define global variables
comp_gene_info = {}
bgc_cog_genes = {}
nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, codo_plo_dir = [None]*5
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
rscript_for_plotting = dir_path + 'plotCogConservation.R'

gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def parseCodonAlignmentStats(cog, codon_alignment_fasta,  phylogenetic_distances):
    domain_plot_file = codo_plo_dir + cog + '_domain.txt'
    position_plot_file = codo_plo_dir + cog + '_position.txt'
    popgen_plot_file = codo_plo_dir + cog + '_popgen.txt'
    plot_pdf_file = codo_plo_dir + cog + '.pdf'

    domain_plot_handle = open(domain_plot_file, 'w')
    position_plot_handle = open(position_plot_file, 'w')
    popgen_plot_handle = open(popgen_plot_file, 'w')

    seqs = []
    bgc_codons = defaultdict(list)
    num_codons = None
    samples = set([])
    gene_lengths = []
    gene_locs = defaultdict(dict)
    core_counts = defaultdict(int)
    products = set([])
    with open(codon_alignment_fasta) as ocaf:
        for rec in SeqIO.parse(ocaf, 'fasta'):
            sample_id, gene_id = rec.id.split('|')
            if comp_gene_info[gene_id]['core_overlap']: core_counts['core'] += 1
            else: core_counts['auxiliary'] += 1
            products.add(comp_gene_info[gene_id]['product'])
            real_pos = 1
            seqs.append(list(str(rec.seq)))
            codons = [str(rec.seq)[i:i+3] for i in range(0, len(str(rec.seq)), 3)]
            num_codons = len(codons)
            bgc_codons[rec.id] = codons
            samples.add(sample_id)
            for msa_pos, bp in enumerate(str(rec.seq)):
                if bp != '-':
                    gene_locs[gene_id][real_pos] = msa_pos+1
                    real_pos += 1
            gene_lengths.append(len(str(rec.seq).replace('-', '')))

    is_core = False
    if float(core_counts['core'])/sum(core_counts.values()) >= 0.8: is_core = True

    median_gene_length = statistics.median(gene_lengths)

    tot_phylogenetic_breadth = 0
    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            if i < j and s1 in phylogenetic_distances and s2 in phylogenetic_distances:
                tot_phylogenetic_breadth += phylogenetic_distances[s1][s2]

    variable_sites = set([])
    conserved_sites = set([])
    position_plot_handle.write('\t'.join(['pos', 'num_seqs', 'num_alleles', 'num_gaps', 'maj_allele_freq']) + '\n')
    for i, ls in enumerate(zip(*seqs)):
        al_counts = defaultdict(int)
        for al in ls:
            if al != '-': al_counts[al] += 1
        maj_allele_count = max(al_counts.values())
        tot_count = sum(al_counts.values())
        num_seqs = len(ls)
        num_alleles = len(al_counts.keys())
        num_gaps = num_seqs - tot_count
        maj_allele_freq = float(maj_allele_count)/tot_count
        position_plot_handle.write('\t'.join([str(x) for x in [i+1, num_seqs, num_alleles, num_gaps, maj_allele_freq]]) + '\n')
        if maj_allele_freq <= 0.90:
            variable_sites.add(i)
        else:
            conserved_sites.add(i)
    position_plot_handle.close()

    differential_domains = set([])
    domain_positions_msa = defaultdict(set)
    domain_min_position_msa = defaultdict(lambda: 100000)
    all_domains = set([])
    for gene in bgc_cog_genes[cog]:
        gene_start = comp_gene_info[gene]['start']
        gene_end = comp_gene_info[gene]['end']
        for domain in comp_gene_info[gene]['gene_domains']:
            domain_start = max(domain['start'], gene_start)
            domain_end = min(domain['end'], gene_end)
            domain_name = domain['aSDomain'] + '_|_' + domain['description']
            relative_start = domain_start-gene_start
            assert(len(gene_locs[gene])+3 >= (domain_end-gene_start))
            relative_end = min([len(gene_locs[gene]), domain_end-gene_start])
            domain_range = range(relative_start, relative_end)
            for pos in domain_range:
                msa_pos = gene_locs[gene][pos+1]
                domain_positions_msa[domain_name].add(msa_pos)
                if domain_min_position_msa[domain_name] > msa_pos:
                    domain_min_position_msa[domain_name] = msa_pos

            domain_variable = variable_sites.intersection(domain_positions_msa[domain_name])
            domain_conserved = conserved_sites.intersection(domain_positions_msa[domain_name])
            other_variable = variable_sites.difference(domain_positions_msa[domain_name])
            other_conserved = conserved_sites.difference(domain_positions_msa[domain_name])

            oddsratio, pvalue = fisher_exact([[len(domain_variable), len(domain_conserved)], [len(other_variable), len(other_conserved)]])
            if pvalue < 0.001: differential_domains.add(domain_name)
            all_domains.add(domain['type'] + '_|_' + domain['aSDomain'] + '_|_' + domain['description'])

    domain_plot_handle.write('\t'.join(['domain', 'domain_index', 'min_pos', 'max_pos']) + '\n')
    for i, dom in enumerate(sorted(domain_min_position_msa.items(), key=itemgetter(1))):
        tmp = []
        old_pos = None
        for j, pos in enumerate(sorted(domain_positions_msa[dom[0]])):
            if j == 0:
                old_pos = pos-1
            if pos-1 != old_pos:
                if len(tmp) > 0:
                    min_pos = min(tmp)
                    max_pos = max(tmp)
                    domain_plot_handle.write('\t'.join([str(x) for x in [dom[0], i, min_pos, max_pos]]) + '\n')
                tmp = []
            tmp.append(pos)
            old_pos = pos
        if len(tmp) > 0:
            min_pos = min(tmp)
            max_pos = max(tmp)
            domain_plot_handle.write('\t'.join([str(x) for x in [dom[0], i, min_pos, max_pos]]) + '\n')
    domain_plot_handle.close()

    a = egglib.io.from_fasta(codon_alignment_fasta, alphabet=egglib.alphabets.DNA)
    struct = egglib.get_structure(a)
    cs = egglib.stats.ComputeStats()
    cs.set_structure(struct)
    cs.add_stats('S', 'thetaW', 'Pi', 'D', 'lseff', 'nseff')
    cs_stats = cs.process_align(a)
    a.to_codons()
    cd = egglib.stats.CodingDiversity()
    cd.process(a, skipstop=True, max_missing=0.0, multiple_alleles=False, multiple_hits=False)
    effective_dn_positions = set([x.position for x in cd.positions_NS])
    effective_ds_positions = set([x.position for x in cd.positions_S])
    effective_dn_ds = "NA"
    if cd.num_pol_S > 0: effective_dn_ds =  float(cd.num_pol_NS) / cd.num_pol_S

    popgen_plot_handle.write('\t'.join(['pos', 'type', 'effective']) + '\n')
    core_codons = 0
    total_variable_codons = 0
    nonsynonymous_sites = 0
    synonymous_sites = 0
    for cod_index in range(0, num_codons):
        first_bp = (cod_index+1)*3
        second_bp = ((cod_index + 1)*3) + 1
        third_bp = ((cod_index + 1)*3) + 2
        cod_positions = set([first_bp, second_bp, third_bp])
        aa_count = defaultdict(int)
        cod_count = defaultdict(int)
        core = True
        for bgc in bgc_codons:
            cod = bgc_codons[bgc][cod_index]
            if '-' in cod or 'N' in cod:
                core = False
            else:
                aa_count[gencode[cod]] += 1
                cod_count[cod] += 1
        maj_allele_count = max(cod_count.values())
        tot_valid_codons = sum(cod_count.values())
        residues = len(aa_count.keys())
        maj_allele_freq = float(maj_allele_count)/tot_valid_codons

        nonsyn_flag = False; syn_flag = False
        if maj_allele_freq <= 0.9:
            total_variable_codons += 1
            if len(cod_count.keys()) > 1:
                if residues > 1: nonsynonymous_sites += 1; nonsyn_flag = True
                else: synonymous_sites += 1; syn_flag = True
        if core: core_codons += 1
        if nonsyn_flag or syn_flag:
            type = 'S'
            effective = 'False'
            if nonsyn_flag:
                type = 'NS'
                if len(cod_positions.intersection(effective_dn_positions)) > 0: effective = 'True'
            elif len(cod_positions.intersection(effective_ds_positions)): effective = 'True'
            popgen_plot_handle.write('\t'.join([str(x) for x in [first_bp, type, effective]]) + '\n')
    popgen_plot_handle.close()
    dn_ds = "NA"
    if synonymous_sites > 0: dn_ds =  float(nonsynonymous_sites) / synonymous_sites

    os.system('Rscript %s %s %s %s %s' % (rscript_for_plotting, domain_plot_file, position_plot_file, popgen_plot_file, plot_pdf_file))
    return (
        ['; '.join(products), is_core, median_gene_length, len(seqs), len(samples), cs_stats['D'], core_codons, total_variable_codons,
         nonsynonymous_sites, synonymous_sites,
         dn_ds, cd.num_codons_eff, cd.num_pol_NS, cd.num_pol_S,
         effective_dn_ds, tot_phylogenetic_breadth, '; '.join(differential_domains), '; '.join(all_domains)])

def create_msas(cog):
    cog_nucl_fasta = nucl_seq_dir + '/' + cog + '.fna'
    cog_prot_fasta = prot_seq_dir + '/' + cog + '.faa'
    cog_prot_msa = prot_alg_dir + '/' + cog + '.msa.faa'
    cog_codo_msa = codo_alg_dir + '/' + cog + '.msa.fna'

    cog_nucl_handle = open(cog_nucl_fasta, 'w')
    cog_prot_handle = open(cog_prot_fasta, 'w')
    for gene in bgc_cog_genes[cog]:
        cog_nucl_handle.write(
            '>' + comp_gene_info[gene]['bgc_name'] + '|' + gene + '\n' + comp_gene_info[gene]['nucl_seq'] + '\n')
        cog_prot_handle.write(
            '>' + comp_gene_info[gene]['bgc_name'] + '|' + gene + '\n' + comp_gene_info[gene]['prot_seq'] + '\n')
    cog_nucl_handle.close()
    cog_prot_handle.close()

    os.system('mafft --maxiterate 1000 --localpair %s > %s' % (cog_prot_fasta, cog_prot_msa))
    os.system('pal2nal.pl %s %s -output fasta > %s' % (cog_prot_msa, cog_nucl_fasta, cog_codo_msa))

def runCoinfinder(ortholog_matrix_file, species_phylogeny_file, samples_in_phylogeny, resdir, cores):
    input_reformatted_file = resdir + '/coinfinder_input.txt'
    input_reformatted_handle = open(input_reformatted_file, 'w')
    sample_names = []
    with open(ortholog_matrix_file) as omf:
        for i, line in enumerate(omf):
            line = line.strip('\n')
            ls = line.split('\t')
            if i == 0:
                sample_names = ls[1:-1]
            else:
                cog = ls[0]
                for j, val in enumerate(ls[1:-1]):
                    if val != '':
                        if sample_names[j] in samples_in_phylogeny:
                            input_reformatted_handle.write(cog + '\t' + sample_names[j] + '\n')
    input_reformatted_handle.close()

    coinfinder_result_prefix = resdir + 'coinfinder'
    os.system('coinfinder -i %s -p %s -a -d -m -t -x %d -o %s' % (input_reformatted_file, species_phylogeny_file, cores, coinfinder_result_prefix))

    coinfinder_results_file = resdir + '/coinfinder_pairs.tsv'
    coevolved_cogs = defaultdict(set)
    with open(coinfinder_results_file) as ocrf:
        for i, line in enumerate(ocrf):
            if i == 0: continue
            line = line.strip()
            ls = line.split('\t')
            cog1, cog2 = ls[:2]
            if cog1 in bgc_cog_genes.keys() and not cog2 in bgc_cog_genes.keys(): coevolved_cogs[cog1].add(cog2)
            if cog2 in bgc_cog_genes.keys() and not cog1 in bgc_cog_genes.keys(): coevolved_cogs[cog2].add(cog1)
    return(coevolved_cogs)

def getSpeciesRelationshipsFromPhylogeny(species_phylogeny, samples_in_gcf):
    samples_in_phylogeny = set([])
    t = Tree(species_phylogeny)
    for leaf in t:
        samples_in_phylogeny.add(str(leaf).strip('\n').lstrip('-'))

    pairwise_distances = defaultdict(lambda: defaultdict(float))
    for s1 in samples_in_gcf.intersection(samples_in_phylogeny):
        for s2 in samples_in_gcf.intersection(samples_in_phylogeny):
            try:
                s1_s2_dist = t.get_distance(s1, s2)
                pairwise_distances[s1][s2] = s1_s2_dist
            except:
                pass
    return ([pairwise_distances, samples_in_phylogeny])

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
                    bgc_info.append({'detection_rule': detection_rule, 'product': product, 'contig_edge': contig_edge})
                elif feature.type == 'proto_core':
                    core_start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    core_end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
                    core_positions = core_positions.union(set(range(core_start, core_end+1)))

            full_sequence = str(rec.seq)
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

def bgreport(bgc_specs_file, orthofinder_matrix, outdir, species_phylogeny, cores):

    ### vet input files quickly
    sys.stderr.write('Checking if input files exist and are in valid formatting...\n')
    try:
        if species_phylogeny:
            assert(os.path.isfile(species_phylogeny))
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
    global comp_gene_info
    bgc_genes = {}
    bgc_information = {}
    all_gene_lts = set([])
    cog_order_scores = defaultdict(int)
    for i, gbk in enumerate(antismash_genbanks):
        bgc_name = edited_sample_names[i]
        gene_info, bgc_info = parseGenbanksForGeneInfo(gbk, bgc_name)
        comp_gene_info.update(gene_info)
        bgc_genes[bgc_name] = gene_info.keys()
        bgc_information[bgc_name] = bgc_info
        all_gene_lts = all_gene_lts.union(set(gene_info.keys()))

    # parse orthofinder matrix
    global bgc_cog_genes
    gene_to_cog, bgc_cog_genes, cog_median_counts = parseOrthoFinderMatrix(orthofinder_matrix, all_gene_lts)

    # write the rest of the iTol track file for illustrating genes across BGC instances
    ref_cog_directions = {}
    for i, s in enumerate(bgc_genes):
        samp_genes = bgc_genes[s]
        cog_directions = {}
        cog_starts = {}
        for g in samp_genes:
            ginfo = comp_gene_info[g]
            if g in gene_to_cog:
                cog = gene_to_cog[g]
                cog_directions[cog] = ginfo['direction']
                cog_starts[cog] = ginfo['start']
        reverse_flag = False
        if i == 0:
            ref_cog_directions = cog_directions
        else:
            flip_support = 0
            keep_support = 0
            for c in ref_cog_directions:
                if not c in cog_directions: continue
                if cog_directions[c] == ref_cog_directions[c]: keep_support += 1
                else: flip_support += 1

            # reverse ordering
            if flip_support > keep_support:
                reverse_flag = True
        for c in sorted(cog_starts.items(), key=itemgetter(1), reverse=reverse_flag):
            cog_order_scores[c[0]] += c[1]

    # parse phylogenetic distances between isolates from newick file
    pairwise_distances, samples_in_phylogeny = getSpeciesRelationshipsFromPhylogeny(species_phylogeny, set(sample_naming.values()))

    # run coinfinder to find COGs under co-evolution/linked
    coinfinder_dir = os.path.abspath(outdir + 'CoinFinder/') + '/'
    if not os.path.isdir(coinfinder_dir): os.system('mkdir %s' % coinfinder_dir)
    coevolved_cogs = runCoinfinder(orthofinder_matrix, species_phylogeny, samples_in_phylogeny, coinfinder_dir, cores)

    global nucl_seq_dir; global prot_seq_dir; global prot_alg_dir; global codo_alg_dir; global codo_plo_dir
    nucl_seq_dir = os.path.abspath(outdir + 'Nucleotide_Sequences') + '/'
    prot_seq_dir = os.path.abspath(outdir + 'Protein_Sequences') + '/'
    prot_alg_dir = os.path.abspath(outdir + 'Protein_Alignments') + '/'
    codo_alg_dir = os.path.abspath(outdir + 'Codon_Alignments') + '/'
    codo_plo_dir = os.path.abspath(outdir + 'Codon_MSA_Plots') + '/'
    if not os.path.isdir(nucl_seq_dir): os.system('mkdir %s' % nucl_seq_dir)
    if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
    if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
    if not os.path.isdir(codo_alg_dir): os.system('mkdir %s' % codo_alg_dir)
    if not os.path.isdir(codo_plo_dir): os.system('mkdir %s' % codo_plo_dir)

    p = multiprocessing.Pool(cores)
    p.map(create_msas, list(bgc_cog_genes.keys()))

    cog_stats = {}
    for f in os.listdir(codo_alg_dir):
        cog = f.split('.msa.fna')[0]
        codon_alignment_fasta = codo_alg_dir + f
        parsed_stats = parseCodonAlignmentStats(cog, codon_alignment_fasta, pairwise_distances)
        cog_stats[cog] = parsed_stats

    final_output_handle = open(outdir + 'Ortholog_Group_Information.txt', 'w')
    final_output_handle.write('\t'.join(['cog', 'cog_order_index', 'cog_median_copy_count', 'annotation', 'is_core', 'median_gene_length', 'bgcs_with_cog', 'samples_with_cog', 'Tajimas_D', 'core_codons', 'total_variable_codons', 'nonsynonymous_codons', 'synonymous_codons', 'dn_ds', 'total_effective_codons', 'nonsynonymous_effective_codons', 'synonymous_effective_codons', 'effective_dn_ds', 'cog_phylogenetic_breadth',  'differentially_variable_domains', 'all_domains', 'coinfinder_coevolved_genes']) + '\n')
    for cog, order_rank in sorted(cog_order_scores.items(), key=itemgetter(1)):
        final_output_handle.write('\t'.join([str(x) for x in [cog, order_rank, cog_median_counts[cog]] + cog_stats[cog] + ['; '.join(coevolved_cogs[cog])]]) + '\n')
    final_output_handle.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Program to report statistics on ortholog groups found within BGCs of GCF.
    """)

    parser.add_argument('-b', '--bgc_specs_file', help='BGC specifications file. Tab delimited: 1st column contains path to AntiSMASH BGC Genbank and 2nd column contains sample name.', required=True)
    parser.add_argument('-m', '--orthofinder_matrix', help="OrthoFinder matrix.", required=True)
    parser.add_argument('-o', '--output_directory', help="Prefix for output files.", required=True)
    parser.add_argument('-s', '--species_phylogeny', help="The species phylogeny in Newick format.", required=True)
    parser.add_argument('-c', '--cores', type=int, help="The number of cores to use.", required=False, default=1)
    args = parser.parse_args()

    bgselect(args.bgc_specs_file, args.orthofinder_matrix, args.output_directory, args.species_phylogeny, args.cores)