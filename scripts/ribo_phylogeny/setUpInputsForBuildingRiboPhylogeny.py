import os
import sys

sample_list = set([x.strip() for x in open(sys.argv[1]).readlines()])
prokka_dirs = sys.argv[2:]

proteome_dir = os.path.abspath('Proteomes') + '/'
prokka_dir = os.path.abspath('Prokka') + '/'

for pkd in prokka_dirs:
    pkd = os.path.abspath(pkd) + '/'
    pro_dir = pkd + '/Prokka_Proteomes/' 
    for f in os.listdir(pro_dir):
        s = f.split('.faa')[0]
        if s in sample_list:
            os.system('ln -s %s %s' % (pro_dir + f, proteome_dir))
    for s in os.listdir(pkd):
        if s in sample_list:
            os.system('ln -s %s %s' % (pkd + s, prokka_dir))

