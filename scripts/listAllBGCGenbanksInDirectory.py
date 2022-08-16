#!/usr/bin/env python

### Program: listAllBGCGenbanksInDirectory.py
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
import argparse
import glob


def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: listAllBGCGenbanksInDirectory.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Program to create BGC Genbanks listing file needed for lsaBGC-Ready.py. Provided a directory with BGC prediction 
	results (from antiSMASH, DeepBGC, or GECCO) for a set of samples, it will create a two-column, tab-delimited listing 
	file where the first column is the sample name and the second is the full path to an individual BGC genbank for the 
	sample. E.g. suppose the following setup:

	./AntiSMASH-General-Dir/Sample-Name-1/Sample-Name-1_Scaffold-1.region0001.gbk
	./AntiSMASH-General-Dir/Sample-Name-1/Sample-Name-1_Scaffold-5.region0007.gbk
	./AntiSMASH-General-Dir/Sample-Name-2/Sample-Name-1_Scaffold-1.region0002.gbk

	Then it will print to standard output the following:

	Sample-Name-1 <tab> /full-path-to/AntiSMASH-General-Dir/Sample-Name-1/Sample-Name-1_Scaffold-1.region0001.gbk
	Sample-Name-1 <tab> /full-path-to/AntiSMASH-General-Dir/Sample-Name-1/Sample-Name-1_Scaffold-5.region0007.gbk
	Sample-Name-2 <tab> /full-path-to/AntiSMASH-General-Dir/Sample-Name-1/Sample-Name-1_Scaffold-1.region0002.gbk
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--input_dir', help='Input directory which contains all BGC prediction results.',
						required=True)
	parser.add_argument('-p', '--bgc_prediction_software',
						help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO). Default is antiSMASH.',
						default='antiSMASH', required=False)
	parser.add_argument('-f', '--filter_incomplete', action='store_true',
						help='Filter out incomplete BGCs (those found on contig edges). Only works for antiSMASH predictions.',
						required=False, default=False)
	args = parser.parse_args()
	return args


def siftAndPrint():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	input_bgc_dir = os.path.abspath(myargs.input_dir) + '/'

	try:
		assert (os.path.isdir(input_bgc_dir))
	except:
		raise RuntimeError('Cannot find input directory with BGC predictions results.')

	filter_incomplete_flag = myargs.filter_incomplete
	bgc_prediction_software = myargs.bgc_prediction_software.upper()

	try:
		assert (bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
	except:
		raise RuntimeError('BGC prediction software option is not a valid option.')

	"""
	START WORKFLOW
	"""

	if bgc_prediction_software == 'ANTISMASH':
		for full_file_name in glob.glob(input_bgc_dir + "*/*region*.gbk"):
			sample = full_file_name.split('/')[-2]
			contig_edge_flag = False
			with open(full_file_name) as offn:
				for line in offn:
					line = line.strip()
					if '/contig_edge="True"' in line:
						contig_edge_flag = True

			if not filter_incomplete_flag or (filter_incomplete_flag and not contig_edge_flag):
				sample = sample.replace('#', '').replace('*', '_').replace(':', '_').replace(';', '_').replace(' ', '_').replace(':', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_").replace('-', '_').replace('(', '').replace(')', '').replace('/', '').replace('\\', '').replace('[', '').replace(']', '').replace(',', '')
				print(sample + '\t' + full_file_name)
	elif bgc_prediction_software == 'DEEPBGC':
		for full_file_name in glob.glob(input_bgc_dir + "*/*bgc.gbk"):
			sample = full_file_name.split('/')[-2]
			sample = sample.replace('#', '').replace('*', '_').replace(':', '_').replace(';', '_').replace(' ', '_').replace(':', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_").replace('-', '_').replace('(', '').replace(')', '').replace('/', '').replace('\\', '').replace('[', '').replace(']', '').replace(',', '')
			print(sample + '\t' + full_file_name)
	elif bgc_prediction_software == 'GECCO':
		for full_file_name in glob.glob(input_bgc_dir + "*/*_cluster_*.gbk"):
			sample = full_file_name.split('/')[-2]
			sample = sample.replace('#', '').replace('*', '_').replace(':', '_').replace(';', '_').replace(' ', '_').replace(':', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_").replace('-', '_').replace('(', '').replace(')', '').replace('/', '').replace('\\', '').replace('[', '').replace(']', '').replace(',', '')
			print(sample + '\t' + full_file_name)


if __name__ == '__main__':
	siftAndPrint()