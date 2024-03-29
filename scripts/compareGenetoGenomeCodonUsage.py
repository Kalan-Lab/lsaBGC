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
	
	This program compares the codon distribution of a single gene (provided as a Genbank) to the codon usage of the 
	background genome (predicted gene/ORF sequences - e.g. *.ffn ending files from Prokka annotation). It will report
	the cosine distance and Spearman correlation between the two profiles. Only ORFs which are of length divisible 
	by 3 will be considered.
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--prokka_ffn', help="Path to predicted transcriptome for isolate's genome.", required=True, default=None)
    parser.add_argument('-i', '--gene_id', help="Path to BGC Genbanks for isolate. Each Genbank should reference CDS locus tag IDs matching predicted transcriptome file.", required=True, default=None)
    parser.add_argument('-o', '--output', help="Path to output file.", required=True, default=None)
    args = parser.parse_args()
    return args

valid_bases = set(['A', 'C', 'G', 'T'])
def main():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """
    myargs = create_parser()

    prokka_ffn = os.path.abspath(myargs.prokka_ffn)
    gene_id = myargs.gene_id
    output = os.path.abspath(myargs.output)

    try:
        assert(os.path.isfile(prokka_ffn))
    except:
        raise RuntimeError('One or more input files do not exist. Exiting now ...')

    """
    START WORKFLOW
    """

    # parse Prokka predicted transcriptome
    cod_freq_dict_gcf = defaultdict(int)
    cod_freq_dict_background = defaultdict(int)
    all_cods = set([])
    with open(prokka_ffn) as offn:
        for rec in SeqIO.parse(offn, 'fasta'):
            locus_tag = rec.id
            if not len(str(rec.seq))%3 == 0:
                print(locus_tag)
                continue
            codon_seq = [str(rec.seq)[i:i + 3] for i in range(0, len(str(rec.seq)), 3)]
            for cod in list(codon_seq):
                if not(len(cod) == 3 and cod[0] in valid_bases and cod[1] in valid_bases and cod[2] in valid_bases): continue
                if locus_tag == gene_id.strip():
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
    output_handle.write('Spearman_Pvalue\t%f\n' % round(spm_pval, 3))
    output_handle.write('GCF_Codons\t%s\n' % ', '.join(cod_order))
    output_handle.write('GCF_Codon_Frequencies\t%s\n' % ', '.join([str(x) for x in gcf_cod_freqs]))
    output_handle.write('Background_Codon_Frequencies\t%s\n' % ', '.join([str(x) for x in bkg_cod_freqs]))
    output_handle.close()

    sys.exit(0)

if __name__ == '__main__':
    main()
