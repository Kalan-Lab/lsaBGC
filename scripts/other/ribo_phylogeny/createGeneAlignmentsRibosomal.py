import os
import sys
from Bio import SeqIO
from collections import defaultdict

hmmdir = os.path.abspath(sys.argv[1])+ '/'
annotdir = os.path.abspath(sys.argv[2]) + '/'
prokkdir = os.path.abspath(sys.argv[3]) + '/'
outdir = os.path.abspath(sys.argv[4]) + '/'

protdir = outdir + 'Protein_Seqs/'
algndir = outdir + 'Protein_Algn/'
cdsdir = outdir + 'CDS_Seqs/'
coddir = outdir + 'Codon_Algn/'

os.system('mkdir %s %s %s %s' % (protdir, algndir, cdsdir, coddir))

pep_seqs = {}
cds_seqs = {}
gene_to_strain = {}
all_samps = set([])
for f in os.listdir(annotdir):
        if not f.endswith('.faa'): continue
        s = f.split('.faa')[0]
	
        s_prokdir = prokkdir + s + '/'
        cdsf = [s_prokdir + f for f in os.listdir(s_prokdir) if f.endswith('.ffn')][0]
        pdsf = annotdir + f

        assert(os.path.isfile(cdsf))
        if os.path.isfile(cdsf):
            with open(cdsf) as ocf:
                for rec in SeqIO.parse(ocf, 'fasta'):
                    cds_seqs[rec.id] = str(rec.seq)
                    gene_to_strain[rec.id] = s
                    
            with open(pdsf) as opf:
                for rec in SeqIO.parse(opf, 'fasta'):
                    pep_seqs[rec.id] = str(rec.seq)
                    all_samps.add(s)

all_rps = set([])
best_hits = defaultdict(lambda: defaultdict(lambda: [None, 10000.0]))
for f in os.listdir(hmmdir):
    if not '_with_' in f: continue
    samp = f.split('_with_')[0]
    hmm = f.split('_with_')[1].split('.out')[0]
    all_rps.add(hmm)
    with open(hmmdir + f) as ohf:
        for line in ohf:
            line = line.strip()
            if line.startswith("#"): continue
            ls = line.split()
            hit = ls[0]
            evalue = float(ls[4])
            if best_hits[hmm][samp][1] > evalue and evalue < 1e-10: best_hits[hmm][samp] = [hit, evalue]

for i, rp in enumerate(best_hits):
    pepf = protdir + rp + '.faa'
    cdsf = cdsdir + rp + '.fna'
    pepaf = algndir + rp + '.msa.faa'
    cdsaf = coddir  + rp + '.msa.fna'

    opepf = open(pepf, 'w')
    ocdsf = open(cdsf, 'w') 
    for s in all_samps:
        bh = best_hits[rp][s][0]
        if bh != None:
            opepf.write('>' + s + '\n' + pep_seqs[bh] + '\n')
            ocdsf.write('>' + s + '\n' + cds_seqs[bh] + '\n')
    opepf.close()
    ocdsf.close()

    os.system('mafft-linsi %s > %s' % (pepf, pepaf))
    os.system('pal2nal.pl %s %s -output fasta > %s' % (pepaf, cdsf, cdsaf))
