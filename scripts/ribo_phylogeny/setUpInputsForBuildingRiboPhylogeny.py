import os
import sys

sample_list = set([x.strip() for x in open(sys.argv[1]).readlines()]) # text file with each line corresponding to sample name
prokka_dirs = sys.argv[2:] # a list of arguments which list the prokka subdirectories from lsaBGC-AutoProcess.py analysis (provided in order of preference)

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

