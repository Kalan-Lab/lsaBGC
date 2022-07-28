#!/usr/bin/env python

### Program: simpleProteinExtraction.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2021, Kalan-Lab
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import defaultdict
from lsaBGC import util
import subprocess
from operator import itemgetter

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: simpleProteinExtraction.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Simple script to extract proteins from full Genbank with CDS features.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--input_genbank', help='Path to input Genbank with CDS features.', required=True)
	parser.add_argument('-o', '--output_fasta', help='Path to output FASTA.', required=True)

	args = parser.parse_args()
	return args

def simpleProteinExtraction():
	"""
	Void function which runs primary workflow for program.
	"""

	myargs = create_parser()

	input_genbank_file = myargs.input_genbank
	output_fasta_file = myargs.output_fasta

	try:
		protein_dict = {}
		pid = 0
		with open(input_genbank_file) as oigf:
			for rec in SeqIO.parse(oigf, 'genbank'):
				for feature in rec.features:
					if feature.type == 'CDS':
						protein_dict[pid] = str(feature.qualifiers.get('translation')[0]).replace('*', '')
						pid += 1

		if len(protein_dict) > 0:
			outf = open(output_fasta_file, 'w')
			for p in protein_dict:
				outf.write('>' + str(p) + '\n' + protein_dict[p] + '\n')
			outf.close()
	except:
		raise RuntimeError('Issues processing Genbank and extraction proteins into output FASTA file.')

if __name__ == '__main__':
	simpleProteinExtraction()