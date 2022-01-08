import os
import sys
from collections import defaultdict

disco_dir = os.path.abspath(sys.argv[1]) + '/' # disco directory from lsaBGC-AutoAnalyze.py
bowtie2_ref = os.path.abspath(sys.argv[2]) # mega reference
res_dir = os.path.abspath(sys.argv[3]) + '/'  # results directory for alignment of novel snv fastq to mega reference
for g in os.listdir(disco_dir):
    alg_dir = disco_dir + g + '/Alignment_Parsing/'
    if not os.path.isdir(alg_dir): continue
    g_out_dir = res_dir + g + '/' 
    os.system('mkdir %s' % g_out_dir)
    for f in os.listdir(alg_dir):
        if not f.endswith('snv_support.fastq.gz'): continue
        s = f.split('.snv_support.fastq.gz')[0] 
        sam_file = g_out_dir + s + '.sam'
        bam_file = g_out_dir + s + '.bam'
        bam_file_sorted = g_out_dir + s + '.sorted.bam'

        bowtie2_cmd = ['bowtie2', '--very-sensitive-local', '--no-unal', '-a', '-x', bowtie2_ref, '-U', alg_dir + f, '-S', sam_file, '-p', '40']
        samtools_view_cmd = ['samtools', 'view', '-h', '-Sb', sam_file, '>', bam_file]
        samtools_sort_cmd = ['samtools', 'sort', '-@', '40', bam_file, '-o', bam_file_sorted]
        samtools_index_cmd = ['samtools', 'index', bam_file_sorted]
        
        os.system(' '.join(bowtie2_cmd))
        os.system(' '.join(samtools_view_cmd))
        os.system(' '.join(samtools_sort_cmd))
        os.system(' '.join(samtools_index_cmd))
        os.system("rm -f %s %s" % (sam_file, bam_file))
