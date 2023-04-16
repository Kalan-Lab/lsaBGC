#!/usr/bin/env python

### Program: visualize_BGC-Ome.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2022, Kalan-Lab
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
import shutil
import sys
import argparse
import subprocess
from collections import defaultdict
from Bio import SeqIO
from lsaBGC import util
import _pickle as cPickle
import traceback

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_BGCOME = lsaBGC_main_directory + '/lsaBGC/Rscripts/VisualizeBGCome.R'
gecco_pickle_weights_file_file = lsaBGC_main_directory + '/db/GECCO_PF_Weights.pkl'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: visualize_BGC-Ome.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	*******************************************************************************************************************
	Visualize the BGC-ome for a sample.	
	
	Visualization is performed using the gggenes library in R. 
	
	Coloring basis:
	
	- For antiSMASH predictions: Discrete. "rule-based" key domain for detection containing genes are gold, rest are grey.  
	- For DeepBGC predictions: Continuous. score by model for BGC-likeness. 
	- For GECCO predictions: Continuous. score by model for BGC-likeness.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-b', '--bgc_results_dir', help='Sample-specific results directory from antiSMASH, DeepBGC, or GECCO.', required=True)
	parser.add_argument('-p', '--bgc_prediction_software', help='Software used to predict BGCs (Options: antiSMASH, DeepBGC, GECCO) [Default is antiSMASH].', default='antismash', required=False)
	parser.add_argument('-o', '--output_pdf', help='Path to output PDF file [Default is ./Sample_BGC-ome.pdf].')
	parser.add_argument('-t', '--tmp_dir', help='Path to temporary dir. [Default is ./tmp_dir_$INDIR_BASE_NAME].', required=False, default=None)
	parser.add_argument('-l', '--length', type=int, help='Specify the height/length of the heatmap plot. [Default is 7].', required=False, default=7)
	parser.add_argument('-w', '--width', type=int, help='Specify the width of the heatmap plot. [Default is 10].', required=False, default=14)

	args = parser.parse_args()
	return args

def BGCome():
	myargs = create_parser()

	bgc_results_dir = os.path.abspath(myargs.bgc_results_dir) + '/'
	bgc_prediction_software = myargs.bgc_prediction_software.upper()
	output_pdf = myargs.output_pdf
	tmp_dir = myargs.tmp_dir
	length = myargs.length
	width = myargs.width
	curr_dir = os.path.abspath(os.getcwd()) + '/'

	try:
		assert (bgc_prediction_software in set(['ANTISMASH', 'DEEPBGC', 'GECCO']))
	except:
		sys.stderr.write('BGC prediction software option is not a valid option.\n')
		sys.exit(1)

	try:
		if tmp_dir == None:
			res_id = bgc_results_dir.split('/')[-2]
			tmp_dir = curr_dir + 'tmp_' + res_id + '/'
			assert(not os.path.isdir(tmp_dir))
			os.mkdir(tmp_dir)
			assert(os.path.isdir(tmp_dir))
		else:
			assert (not os.path.isdir(tmp_dir))
			os.mkdir(tmp_dir)
			assert (os.path.isdir(tmp_dir))
	except Exception as e:
		sys.stderr.write('Error: tmp directory already exists or unable to create it!\n')
		sys.exit(1)

	try:
		track_file = tmp_dir + 'CDS_Features.txt'
		if bgc_prediction_software == 'ANTISMASH':
			scaffold_gbks = defaultdict(list)
			for f in os.listdir(bgc_results_dir):
				if not '.region' in f or not f.endswith('.gbk'): continue
				scaff = f.split('.region')[0]
				gbk_path = bgc_results_dir + f
				scaffold_gbks[scaff].append(gbk_path)

			tf_handle = open(track_file, 'w')
			tf_handle.write('\t'.join(['gene_cluster', 'CDS_feature', 'CDS_start', 'CDS_end', 'CDS_dir', 'BGC_likeness']) + '\n')
			for scaff in scaffold_gbks:
				for bgc_gbk in scaffold_gbks[scaff]:
					product = 'unknown'
					bgc_start = 'NA'
					bgc_end = 'NA'
					contig_edge = ""
					try:
						with open(bgc_gbk) as obg:
							for line in obg:
								line = line.strip()
								if line.startswith('Orig. start '):
									bgc_start = line.split('::')[-1].strip()
								elif line.startswith('Orig. end '):
									bgc_end = line.split('::')[-1].strip()
								if '/contig_edge="True"' in line:
									contig_edge = "(near contig edge)"
						cds_features = []
						products = set([])
						with open(bgc_gbk) as obg:
							for rec in SeqIO.parse(obg, 'genbank'):
								for feat in rec.features:
									if feat.type == 'protocluster':
										try:
											products.add(feat.qualifiers.get('product')[0])
										except:
											pass
									elif feat.type == 'CDS':
										all_coords, start, end, direction, is_multi_part = util.parseCDSCoord(str(feat.location))
										rule_based_bgc_cds = '0.0'
										try:
											if 'rule-based-clusters' in feat.qualifiers.get('gene_functions')[0]:
												rule_based_bgc_cds = '1.0'
										except:
											pass
										if direction == '+':
											direction = '1'
										elif direction == '-':
											direction = '0'

										cds_features.append([str(x) for x in [start, end, direction, rule_based_bgc_cds]])

						if len(products) == 1:
							product = list(products)[0]
						elif len(products) == 2 and 'NRPS-like' in products and 'NRPS' in products:
							product = 'NRPS'
						else:
							product = '; '.join(sorted(products))

						bgc_name = bgc_gbk.split('/')[-1].split('.gbk')[0] + ' - ' + product + ' ' + bgc_start + ':' + bgc_end + ' ' + contig_edge

						for cds_it, cds_info in enumerate(cds_features):
							tf_handle.write('\t'.join([bgc_name, 'CDS_' + str(cds_it+1)] + cds_info) + '\n')

					except:
						raise RuntimeError()
			tf_handle.close()
		elif bgc_prediction_software == 'DEEPBGC':
			tf_handle = open(track_file, 'w')
			tf_handle.write('\t'.join(['gene_cluster', 'CDS_feature', 'BGC_likeness', 'CDS_start', 'CDS_end', 'CDS_dir']) + '\n')

			bgc_prediction_gbk_file = [bgc_results_dir + f for f in os.listdir(bgc_results_dir) if f.endswith('.bgc.gbk')][0]

			with open(bgc_prediction_gbk_file) as obpgf:
				for rec in SeqIO.parse(obpgf, 'genbank'):
					scaff = '_'.join(rec.id.split('_')[:-1])
					bgc_start = rec.id.split('_')[-1].split('-')[0]
					bgc_end = rec.id.split('_')[-1].split('-')[1].split('.')[0]
					products = set([])
					for feat in rec.features:
						if not feat.type == 'cluster': continue
						try:
							products.add(feat.qualifiers.get('product_class')[0])
						except:
							pass
					product = 'unknown'
					if len(products) == 1:
						product = list(products)[0]
					elif len(products) >= 2:
						product = 'multi-type'

					bgc_name = scaff+ ' - ' + product + ' - ' + bgc_start + ':' + bgc_end

					cds_it = 1
					for feat in rec.features:
						if not feat.type == 'CDS': continue
						deepbgc_score = "NA"
						try:
							deepbgc_score = feat.qualifiers.get('deepbgc_score')[0]
						except:
							pass
						start = min([int(x) for x in str(feat.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x) for x in str(feat.location)[1:].split(']')[0].split(':')])
						direction = str(feat.location).split('(')[1].split(')')[0]
						if direction == '+':
							direction = '1'
						elif direction == '-':
							direction = '0'
						tf_handle.write('\t'.join([bgc_name, 'CDS_' + str(cds_it), deepbgc_score, str(start), str(end), direction]) + '\n')
						cds_it += 1
			tf_handle.close()

		elif bgc_prediction_software == 'GECCO':
			scaffold_gbks = defaultdict(list)
			bgc_starts = defaultdict(lambda: 'NA')
			bgc_ends = defaultdict(lambda: 'NA')
			for f in os.listdir(bgc_results_dir):
				if f.endswith('.gbk'):
					gbk_path = bgc_results_dir + f
					scaff = f.split('_cluster_')[0]
					scaffold_gbks[scaff].append(gbk_path)
				elif f.endswith('.clusters.tsv'):
					with open(bgc_results_dir + f) as obf:
						for line in obf:
							line = line.strip()
							ls = line.split('\t')
							bgc_id = ls[1]
							bgc_starts[bgc_id] = ls[2]
							bgc_ends[bgc_id] = ls[3]

			tf_handle = open(track_file, 'w')
			tf_handle.write('\t'.join(['gene_cluster', 'CDS_feature', 'BGC_likeness', 'CDS_start', 'CDS_end', 'CDS_dir']) + '\n')
			for scaff in scaffold_gbks:
				for bgc_gbk in scaffold_gbks[scaff]:
					bgc_id = bgc_gbk.split('/')[-1].split('.gbk')[0]
					product = 'unknown'
					bgc_start = bgc_starts[bgc_id]
					bgc_end = bgc_ends[bgc_id]
					rec = SeqIO.read(bgc_gbk, 'genbank')
					product = 'unknown'
					try:
						product = rec.annotations['structured_comment']['GECCO-Data']['biosyn_class']
					except:
						try:
							product = rec.annotations['structured_comment']['GECCO-Data']['cluster_type']
						except:
							pass
					if product == 'Unknown':
						product = 'unknown'

					domains = []
					domain_weights = {}
					gecco_pfam_weights_pickle_handle = open(gecco_pickle_weights_file_file, "rb")
					gecco_pfam_weights = cPickle.load(gecco_pfam_weights_pickle_handle)
					rec = SeqIO.read(bgc_gbk, 'genbank')
					for feat in rec.features:
						if feat.type == 'misc_feature':
							start = feat.location.start + 1
							end = feat.location.end
							aSDomain = "NA"
							dom_weight = -7
							try:
								aSDomain = feat.qualifiers['standard_name'][0]
							except:
								pass
							try:
								dom_weight = gecco_pfam_weights[aSDomain]
							except:
								pass
							domain_weights[aSDomain + '|' + str(start + 1) + '|' + str(end)] = dom_weight
							domains.append({'start': start + 1, 'end': end, 'type': feat.type, 'aSDomain': aSDomain,
											 'is_multi_part': False})

					bgc_name = bgc_id + ' - ' + product + ' - ' + bgc_start + ':' + bgc_end
					cds_it = 1
					for feat in rec.features:
						if feat.type == "CDS":
							start = feat.location.start + 1
							end = feat.location.end
							direction = "-" if feat.location.strand == -1 else "+"

							if direction == '+':
								direction = '1'
							elif direction == '-':
								direction = '0'

							grange = set(range(start, end + 1))

							all_gene_domain_weights = []
							for d in domains:
								drange = set(range(d['start'], d['end'] + 1))
								if len(drange.intersection(grange)) > 0:
									if (d['aSDomain'] + '|' + str(d['start']) + '|' + str(d['end'])) in domain_weights:
										all_gene_domain_weights.append(domain_weights[d['aSDomain'] + '|' + str(d['start']) + '|' + str(d['end'])])

							max_domain_weight = 'NA'
							if len(all_gene_domain_weights) > 0:
								max_domain_weight = str(max(all_gene_domain_weights))
							tf_handle.write('\t'.join([bgc_name, 'CDS_' + str(cds_it), max_domain_weight, str(start), str(end), direction]) + '\n')
							cds_it += 1
			tf_handle.close()
		assert(os.path.isfile(track_file))

		# output PDF
		cmd = ['Rscript', RSCRIPT_FOR_BGCOME, track_file, str(length), str(width), output_pdf]
		try:
			subprocess.call(' '.join(cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
		except Exception as e:
			raise RuntimeError("Had issues running R plotting.")
		assert(os.path.isfile(output_pdf))
		shutil.rmtree(tmp_dir)
	except:
		shutil.rmtree(tmp_dir)
		raise RuntimeError(traceback.format_exc())

if __name__ == '__main__':
	BGCome()
