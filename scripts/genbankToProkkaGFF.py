#!/usr/bin/env python

### Program: genbankToProkkaGFF.py
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
from Bio.Seq import Seq
from lsaBGC import util
import gzip

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: genbankToProkkaGFF.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Process NCBI Genbanks to create proteome + genbanks with specific locus tag.
	
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--genbank', help='Path to input full GenBank file.', required=True)
	parser.add_argument('-o', '--prokka_gff', help='Path to output file in Prokka GFF format.', required=True)

	args = parser.parse_args()
	return args


def reformatToProkkaGFF():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	input_genbank_file = os.path.abspath(myargs.genbank)
	prokka_gff = os.path.abspath(myargs.prokka_gff)

	try:
		assert(util.is_genbank(input_genbank_file))
	except:
		raise RuntimeError('Issue with input Genbank file from NCBI.')

	"""
	START WORKFLOW
	"""

	try:
		outf = open(prokka_gff, 'w')
		outf.write('##gff-version 3\n')
		with open(input_genbank_file) as oigf:
			for rec in SeqIO.parse(oigf, 'genbank'):
				outf.write('##sequence-region ' + str(rec.id) + ' 1 ' + str(len(rec.seq)) + '\n')
		
		with open(input_genbank_file) as oigf:
			for rec in SeqIO.parse(oigf, 'genbank'):
				for feature in rec.features:
					if not feature.type == 'CDS': continue
					start_coord = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
					end_coord = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					direction = str(feature.location).split('(')[1].split(')')[0]
					locus_tag = feature.qualifiers.get('locus_tag')[0]
					last_col = ';'.join(['ID=' + locus_tag, 'inference=ab initio prediction:p(y)rodigal', 'locus_tag=' + locus_tag, 'product=unannotated protein'])
					outf.write('\t'.join([rec.id, 'p(y)rodigal', 'CDS', str(start_coord), str(end_coord), '.', direction, '0', last_col]) + '\n')

		outf.write('##FASTA\n')
		with open(input_genbank_file) as oigf:
			for rec in SeqIO.parse(oigf, 'genbank'):
				outf.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
		outf.close()
	except:
		raise RuntimeError("Issue reformatting to a Prokka GFF file.")

if __name__ == '__main__':
	reformatToProkkaGFF()