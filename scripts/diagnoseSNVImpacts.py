import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import copy

novelty_results_file = sys.argv[1]
cod_alignment_file = sys.argv[2]

gene_codons = {}
with open(cod_alignment_file) as ocaf:
    for line in ocaf:
        line = line.strip()
        og, og_cod_msa = line.split('\t')
        with open(og_cod_msa) as ogcm:
            for rec in SeqIO.parse(ogcm, 'fasta'):
                codons = [str(rec.seq).replace('-', '')[i:i + 3] for i in range(0, len(str(rec.seq).replace('-', '')), 3)]
                gene = rec.id.split('|')[1]
                gene_codons[gene] = codons

sample_changes = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))

with open(novelty_results_file) as onrf:
    for i, line in enumerate(onrf):
        if i == 0: continue
        line = line.strip()
        gcf_id, sample, homolog_group, position_along_msa, alternate_allele, snv_count, reference_sample, reference_gene, reference_position, reference_allele = line.split('\t')
        sample_changes[sample][homolog_group][reference_gene][int(reference_position)] = reference_allele + '_|_' + alternate_allele 

for s in sample_changes:
    for hg in sample_changes[s]:
        for g in sample_changes[s][hg]:
            for ci in range(0, len(gene_codons[g])):
                ref_cod = gene_codons[g][ci]
                alt_cod = copy.deepcopy(ref_cod)
                
                first_bp_ind = (ci*3)+1
                second_bp_ind = (ci*3)+2
                third_bp_ind = (ci*3)+3

                alt_cod = ""
                if first_bp_ind in sample_changes[s][hg][g]:
                    ref_al, alt_al = sample_changes[s][hg][g][first_bp_ind].split('_|_')
                    assert(ref_al == ref_cod[0])
                    alt_cod += alt_al
                else: alt_cod += ref_cod[0]

                if second_bp_ind in sample_changes[s][hg][g]:
                    ref_al, alt_al = sample_changes[s][hg][g][second_bp_ind].split('_|_')
                    assert(ref_al == ref_cod[1])
                    alt_cod += alt_al
                else: alt_cod += ref_cod[1]

                if third_bp_ind in sample_changes[s][hg][g]:
                    ref_al, alt_al = sample_changes[s][hg][g][third_bp_ind].split('_|_')
                    assert(ref_al == ref_cod[2])
                    alt_cod += alt_al
                else: alt_cod += ref_cod[2]

                if ref_cod != alt_cod:
                    ref_aa = Seq(ref_cod).translate()
                    alt_aa = Seq(alt_cod).translate()

                    syn = True
                    if ref_aa != alt_aa: syn = False
                    
                    print('\t'.join([s, hg, g, str(first_bp_ind), str(ref_cod), str(alt_cod), str(ref_aa), str(alt_aa), str(syn)])) 
