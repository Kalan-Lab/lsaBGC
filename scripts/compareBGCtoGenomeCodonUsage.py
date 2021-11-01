import os
import sys
from Bio import SeqIO
import argparse
from collections import defaultdict
from scipy import stats, spatial

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: compareBGCtoGenomeCodonUsage.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--prokka_ffn', help="Path to predicted transcriptome for isolate's genome.", required=True, default=None)
    parser.add_argument('-b', '--bgc_genbanks', nargs='+', help="Path to BGC Genbanks for isolate. Each Genbank should reference CDS locus tag IDs matching predicted transcriptome file.", required=True, default=None)
    parser.add_argument('-o', '--output', help="Path to output file.", required=True, default=None)
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

    prokka_ffn = os.path.abspath(myargs.prokka_ffn)
    bgc_genbanks = myargs.bgc_genbanks
    output = os.path.abspath(myargs.output)

    print(bgc_genbanks)
    try:
        assert(os.path.isfile(prokka_ffn) and sum([1 for x in bgc_genbanks if os.path.isfile(x)]) == len(bgc_genbanks))
    except:
        raise RuntimeError('One or more input files do not exist. Exiting now ...')

    """
    START WORKFLOW
    """

    # parse Genbanks of BGCs
    gcf_lts = set([])
    for gbk in bgc_genbanks:
        with open(gbk) as ogbk:
            for rec in SeqIO.parse(ogbk, 'genbank'):
                scaffold = rec.id
                for feature in rec.features:
                    if not feature.type == 'CDS': continue
                    locus_tag = feature.qualifiers.get('locus_tag')[0]
                    gcf_lts.add(locus_tag)

    # parse Prokka predicted transcriptome
    cod_freq_dict_gcf = defaultdict(int)
    cod_freq_dict_background = defaultdict(int)
    all_cods = set([])
    with open(prokka_ffn) as offn:
        for rec in SeqIO.parse(off, 'fasta'):
            locus_tag = rec.id
            codon_seq = [str(rec.seq)[i:i + 3] for i in range(0, len(str(rec.seq)), 3)]
            for cod in codon_seq:
                if not(cod[0] in valid_bases and cod[1] in valid_bases and cod[2] in valid_bases): continue
                if locus_tag in gcf_lts:
                    cod_freq_dict_gcf[cod] += 1
                else:
                    cod_freq_dict_background[cod] += 1
                all_cods.add(cod)

    cod_order = []
    gcf_cod_freqs = []
    bkg_cod_freqs = []
    for cod in sorted(all_cods):
        cod_order.append(cod)
        gcf_cod_freqs.append(cod_freq_dict_gcf[cod])
        bkg_cod_freqs.append(cod_freq_dict_background[cod])

    rho, spm_pval = stats.spearmanr(gcf_cod_freqs, bkg_cod_freqs)
    cosine_distance = spatial.distance.cosine(gcf_cod_freqs, bkg_cod_freqs)

    output_handle = open(output, 'w')
    output_handle.write('Cosine_Distance\t%f\n' % round(cosine_distance, 3))
    output_handle.write('Spearman_Rho\t%f\n' % round(rho, 3))
    output_handle.write('Spearman_Pvalue\t%f\n' % round(spmm_pval, 3))
    output_handle.write('GCF_Codons\t%s\n' % ', '.join(cod_order))
    output_handle.write('GCF_Codon_Frequencies\t%s\n' % ', '.join([str(x) for x in gcf_cod_freqs]))
    output_handle.write('Background_Codon_Frequencies\t%s\n' % ', '.join([str(x) for x in gcf_cod_freqs]))
    output_handle.close()

    sys.exit(0)

if __name__ == '__main__':
    main()
