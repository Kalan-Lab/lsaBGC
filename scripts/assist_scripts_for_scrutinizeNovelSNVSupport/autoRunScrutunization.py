import os
import sys

disco_dir = os.path.abspath(sys.argv[1]) + '/'
comp_list_dir = os.path.abspath(sys.argv[2]) + '/'
res_dir = os.path.abspath(sys.argv[3]) + '/'
prog_path = os.path.abspath(sys.argv[4]) 

for f in os.listdir(comp_list_dir):
    g = f.split('.')[0]
    g_disco_dir = disco_dir + g + '/'
    alg_pars_dir = g_disco_dir + 'Alignment_Parsing/'
    pot_nov_snps_file = g_disco_dir + 'Potentially_Novel_SNVSs.txt'
    rep_genes_file = g_disco_dir + 'GCF_Genes_Representatives.fasta'
    print('python ' + prog_path + ' -i ' + pot_nov_snps_file + ' -f ' + alg_pars_dir + ' -r ' + rep_genes_file + ' -o ' + res_dir + g + '.txt -b ' +  comp_list_dir + f)
