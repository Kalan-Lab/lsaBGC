import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def is_genbank(fasta):
	try:
		with open(fasta) as of:
			SeqIO.parse(of, 'genbank')
		return True
	except:
		return False

