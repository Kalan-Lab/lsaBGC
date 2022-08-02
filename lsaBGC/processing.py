# Part of lsaBGC software suite
# Rauf Salamzade
# Kalan Lab, MMI, UW-Madison

# processing.py is a Python module featuring functions for preliminary analysis and processing necessary to use the
# core lsaBGC framework. These functions are primarily referenced in lsaBGC-Process.py and are involved with automating
# gene annotation, secondary metabolite annotation, and de novo ortholog group identification.

import os
import sys
import multiprocessing
import itertools
from lsaBGC import util
from Bio import SeqIO
from collections import defaultdict
import logging
import traceback
import subprocess

def readInAnnotationFilesForExpandedSampleSet(expansion_listing_file, logObject=None):
	"""
	Function to read in Prokka genbank and predicted proteome annotation paths from expansion listing file and load into dictionary with keys corresponding to sample IDs.

	:param expansion_listing_file: tab-delimited file with three columns: (1) sample ID (2) Prokka Genbank path (3) Prokka predicted proteome path.
	:param logObject: python logging object handler.

	:return sample_prokka_data: dictionary of dictionaries with primary keys as sample names and secondary keys as either "genbank" or "predicted_proteome", with final values being paths to corresponding files.
	"""
	sample_prokka_data = defaultdict(dict)
	try:
		with open(expansion_listing_file) as oalf:
			for line in oalf:
				line = line.strip()
				sample, genbank, predicted_proteome = line.split('\t')
				sample = util.cleanUpSampleName(sample)
				try:
					assert (util.is_genbank(genbank))
					assert (util.is_fasta(predicted_proteome))
					sample_prokka_data[sample]['genbank'] = genbank
					sample_prokka_data[sample]['predicted_proteome'] = predicted_proteome
				except Exception as e:
					if logObject:
						logObject.warning('Ignoring sample %s, because at least one of two Prokka annotation files does not seem to exist or be in the expected format.' % sample)
					else:
						sys.stderr.write('Ignoring sample %s, because at least one of two Prokka annotation files does not seem to exist or be in the expected format.' % sample)
		assert (len(sample_prokka_data) >= 1)
		return (sample_prokka_data)
	except Exception as e:
		if logObject:
			logObject.error("Input file listing the location of Prokka annotation files for samples leads to incorrect paths or something else went wrong with processing of it. Exiting now ...")
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def readInAssemblyListing(assembly_listing_file, logObject):
	"""
	Function to read in assembly paths from listing file and load into dictionary with keys corresponding to sample IDs.

	:param assembly_listing_file: tab-delimited file with first column corresponding to sample name and second to genomic assembly path.
	:param logObject: python logging object handler.

	:return sample_assembly_paths: dictionary with keys as sample names and values as genomic assembly paths.
	"""
	sample_assembly_paths = {}
	try:
		with open(assembly_listing_file) as oalf:
			for line in oalf:
				line = line.strip()
				sample, assembly = line.split('\t')
				sample = util.cleanUpSampleName(sample)
				try:
					assert (util.is_fasta(assembly))
					sample_assembly_paths[sample] = assembly
				except:
					logObject.info('Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.' % sample)
					sys.stderr.write('Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.\n' % sample)
		assert (len(sample_assembly_paths) >= 2)
		return (sample_assembly_paths)
	except Exception as e:
		logObject.error("Input file listing the location of assemblies for samples leads to incorrect path was provided. Exiting now ...")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def runProkka(sample_assemblies, prokka_outdir, prokka_proteomes, prokka_genbanks, prokka_load_code,
			  lineage, cpus, locus_tag_length, logObject, dry_run_flag=False, skip_annotation_flag=False):
	"""
	Void function to run Prokka based gene-calling and annotations.

	:param sample_assemblies: dictionary with keys as sample names and values as genomic assembly paths.
	:param prokka_outdir: full path to directory where Prokka results will be written.
	:param prokka_proteomes: full path to directory where Prokka generated predicted-proteome FASTA files will be moved after Prokka has run.
	:param prokka_genbanks: full path to directory where Prokka generated Genbank (featuring predicted CDS) files will be moved after Prokka has run.
	:param prokka_load_code: code to load conda environment for running Prokka
	:param lineage: name of the lineage of interest.
	:param cpus: number of cpus to use in multiprocessing Prokka cmds.
	:param locus_tag_length: length of locus tags to generate using unique character combinations.
	:param logObject: python logging object handler.
	:param dry_run_flag: flag which indicates commands should only be written to task file and not run. This can be used to manually parallelize across an HPC, if available.
	:param skip_annotation_flag: flag which indicates if Prokka should be run on fast mode with many annotation bells and whistles skipped.
	"""
	prokka_cmds = []
	try:
		task_file = prokka_outdir + 'prokka.cmds'
		if dry_run_flag:
			logObject.info("Dry run on: will just be writing Prokka commands to following: %s" % task_file)
			if os.path.isfile(task_file): os.system('rm -f %s' % task_file)
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = list(itertools.product(alphabet, repeat=locus_tag_length))
		for i, sample in enumerate(sample_assemblies):
			sample_assembly = sample_assemblies[sample]
			sample_outdir = prokka_outdir + sample + '/'
			prokka_cmd = None
			if skip_annotation_flag:
				prokka_cmd = [prokka_load_code, 'prokka', '--norrna', '--notrna', '--fast', '--cpus', '1', '--outdir',
							  sample_outdir, '--prefix', sample,
							  '--genus', lineage, '--locustag', ''.join(list(possible_locustags[i])),
							  sample_assembly, ';', 'mv', sample_outdir + sample + '.gbf',
							  prokka_genbanks + sample + '.gbk',
							  ';', 'mv', sample_outdir + sample + '.faa', prokka_proteomes]
			else:
				prokka_cmd = [prokka_load_code, 'prokka', '--cpus', '1', '--outdir', sample_outdir, '--prefix', sample,
							  '--genus', lineage, '--locustag', ''.join(list(possible_locustags[i])),
							  sample_assembly, ';', 'mv', sample_outdir + sample + '.gbf',
							  prokka_genbanks + sample + '.gbk',
							  ';', 'mv', sample_outdir + sample + '.faa', prokka_proteomes]
			prokka_cmds.append(prokka_cmd + [logObject])
			if dry_run_flag:
				task_handle = open(task_file, 'a+')
				task_handle.write(' '.join(prokka_cmd) + '\n')
				task_handle.close()

	except Exception as e:
		logObject.error("Problem with creating commands for running Prokka. Exiting now ...")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	if not dry_run_flag:
		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, prokka_cmds)
		p.close()

def appendSingletonHGsToPresenceMatrix(orthofinder_homolog_matrix, unassigned_orthofinder_homolog_matrix, result_file, logObject):
	"""
	Void function to append singleton HG instances (HGs represneted by only one protein) to the Orthogroups.csv hg by
	sample presence matrix.
	"""
	try:
		logObject.info('Beginning process to append singleton HGs to Orthogroups.csv homolog group by sample presence file to be written to the primary output directory.')

		result_handle = open(result_file, 'w')

		last_hg_identifier = None
		samples_og_header = []
		with open(orthofinder_homolog_matrix) as oohmf:
			for i, line in enumerate(oohmf):
				if i == 0:
					samples_og_header = line.strip('\n').split('\t')
				last_hg_identifier = line.strip('\n').split('\t')[0]
				result_handle.write(line)

		logObject.info("The last homolog group identifier included in the original Orthogroups.csv file was %s" % last_hg_identifier)

		assert(os.path.isfile(unassigned_orthofinder_homolog_matrix))

		samples_uog_header = []
		with open(unassigned_orthofinder_homolog_matrix) as ouohmf:
			for i, line in enumerate(ouohmf):
				if i == 0:
					samples_uog_header = line.strip('\n').split('\t')
				else:
					result_handle.write(line)
		result_handle.close()

		assert(os.path.isfile(result_file))

		for i, so in enumerate(samples_og_header):
			su = samples_uog_header[i]
			assert(su == so)

		logObject.info('Updated Orthogroups.csv file including unassigned/singleton homolog groups can be found at: %s' % result_file)
	except Exception as e:
		logObject.error('Issues with appending singleton homolog group instances to Orthogroups.csv. Now existing.')
		logObject.error(traceback.format_exc())
		raise RuntimeError('Issues with appending singleton homolog group instances to Orthogroups.csv. Now existing.')

def runAntiSMASH(prokka_genbanks_dir, antismash_outdir, antismash_load_code, cpus, logObject, dry_run_flag=False):
	"""
	Void function to run AntiSMASH based annotation of secondary metabolites / biosynthetic gene clusters.

	:param prokka_genbanks: full path to directory where Prokka generated Genbank (featuring predicted CDS) files are located.
	:param antismash_outdir: full path to directory where AntiSMASH results should be written to.
	:param antismash_load_code: code to load conda environment for running AntiSMASH.
	:param cpus: number of cpus to use in multiprocessing/threading AntiSMASH cmds.
	:param logObject: python logging object handler.
	:param dry_run_flag: flag which indicates commands should only be written to task file and not run. This can be used to manually parallelize across an HPC, if available.
	"""
	try:
		task_file = antismash_outdir + 'antismash.cmds'
		if dry_run_flag:
			logObject.info("Dry run on: will just be writing AntiSMASH commands to following: %s" % task_file)
			if os.path.isfile(task_file): os.system('rm -f %s' % task_file)

		asm_cpus = 1
		pool_size = 1
		if cpus < 4:
			asm_cpus = cpus
		else:
			asm_cpus = 4
			pool_size = int(cpus / 4)
		antismash_cmds = []
		for sample_gbk in os.listdir(prokka_genbanks_dir):
			sample = sample_gbk.split('.gbk')[0]
			sample_resdir = antismash_outdir + sample + '/'
			antismash_cmd = [antismash_load_code, 'antismash', '--taxon', 'bacteria', '--genefinding-tool', 'prodigal',
							 '--output-dir', sample_resdir, '--fullhmmer', '--asf', '--cb-general', '--cb-subclusters',
							 '--cb-knownclusters', '-c', str(asm_cpus),
							 prokka_genbanks_dir + sample_gbk]
			antismash_cmds.append(antismash_cmd + [logObject])
			if dry_run_flag:
				task_handle = open(task_file, 'a+')
				task_handle.write(' '.join(antismash_cmd) + '\n')
				task_handle.close()
	except Exception as e:
		logObject.error("Problem with creating commands for running AntiSMASH. Exiting now ...")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	if not dry_run_flag:
		p = multiprocessing.Pool(pool_size)
		p.map(util.multiProcess, antismash_cmds)
		p.close()

def runOrthoFinder(prokka_proteomes_dir, orthofinder_outdir, orthofinder_load_code, cpus, logObject, dry_run_flag=False):
	"""
	Void function to run AntiSMASH based annotation of secondary metabolites / biosynthetic gene clusters.

	:param prokka_genbanks: full path to directory where Prokka generated proteome FASTA files are located.
	:param orthofinder_outdir: full path to directory where OrthoFinder results should be written to.
	:param orthofinder_load_code: code to load conda environment for running OrthoFinder2.
	:param cpus: number of cpus to use in parallelizing OrthoFinder.
	:param logObject: python logging object handler.
	:param dry_run_flag: flag which indicates commands should only be written to task file and not run. This can be used to manually parallelize across an HPC, if available.
	"""
	task_file = orthofinder_outdir + 'orthofinder.cmd'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing OrthoFinder command to following: %s" % task_file)

	orthofinder_cmd = [orthofinder_load_code, 'orthofinder', '-f', prokka_proteomes_dir, '-t', str(cpus)]
	if dry_run_flag:
		task_handle = open(task_file, 'w')
		task_handle.write(' '.join(orthofinder_cmd) + '\n')
		task_handle.close()
	else:
		try:
			logObject.info('Running the following command: %s' % ' '.join(orthofinder_cmd))
			subprocess.call(' '.join(orthofinder_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran OrthoFinder!')
			tmp_orthofinder_dir = os.path.abspath(
				[prokka_proteomes_dir + 'OrthoFinder/' + f for f in os.listdir(prokka_proteomes_dir + 'OrthoFinder/') if f.startswith('Results')][0]) + '/'
			os.system('mv %s %s' % (tmp_orthofinder_dir, orthofinder_outdir))
		except Exception as e:
			logObject.error("Problem with running orthofinder cmd: %s." % ' '.join(orthofinder_cmd))
			logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

def extractBGCProteomes(sample_name, bgc_genbank, bgc_proteomes_outdir, logObject):
	"""
	Function to extract protein sequences from an AntiSMASH Genbank.

	:param sample_name: name of the sample.
	:param bgc_genbank: path to AntiSMASH genbank file
	:param bgc_proteomes_outdir: full path to directory where proteins from Genbank will be written to in FASTA format.
	:param logObject: python logging object handler.

	:return bgc_genbank_proteome: FASTA file with protein sequences from AntiSMASH Genbank.
	"""
	b = bgc_genbank.split('/')[-1].split('.gbk')[0]
	bgc_genbank_proteome = bgc_proteomes_outdir + sample_name + '_BGC_' + b + '.faa'
	try:
		bgc_genbank_proteome_handle = open(bgc_genbank_proteome, 'w')
		with open(bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == "CDS":
						lt = feature.qualifiers.get('locus_tag')[0]
						prot_seq = feature.qualifiers.get('translation')[0]
						prot_id = sample_name + '|' + b + '|' + lt
						bgc_genbank_proteome_handle.write('>' + prot_id + '\n' + str(prot_seq) + '\n')
		bgc_genbank_proteome_handle.close()
		return bgc_genbank_proteome
	except Exception as e:
		logObject.error("Had problems extracting protein sequences of CDS features in genbank %s" % bgc_genbank)
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
