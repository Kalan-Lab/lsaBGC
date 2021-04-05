import os
import sys
import multiprocessing
import itertools
from lsaBGC import util
from Bio import SeqIO
from collections import defaultdict
import subprocess

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
				sample = sample.replace(' ', '_').replace('|', '_').replace('"', '_').replace("'", '_')
				try:
					assert (util.is_fasta(assembly))
					sample_assembly_paths[sample] = assembly
				except:
					logObject.info(
						'Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.' % sample)
					sys.stderr.write(
						'Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.\n' % sample)
		assert (len(sample_assembly_paths) >= 2)
		return (sample_assembly_paths)
	except:
		raise RuntimeError(
			'Input file listing the location of assemblies for samples leads to incorrect path was provided. Exiting now ...')

def runProkka(sample_assemblies, prokka_outdir, prokka_proteomes, prokka_genbanks, prokka_load_code,
			  lineage, cores, locus_tag_length, logObject, dry_run_flag=False, skip_annotation_flag=False):
	"""
	Void function to run Prokka based gene-calling and annotations.

	:param sample_assemblies: dictionary with keys as sample names and values as genomic assembly paths.
	:param prokka_outdir: full path to directory where Prokka results will be written.
	:param prokka_proteomes: full path to directory where Prokka generated predicted-proteome FASTA files will be moved after Prokka has run.
	:param prokka_genbanks: full path to directory where Prokka generated Genbank (featuring predicted CDS) files will be moved after Prokka has run.
	:param prokka_load_code: code to load conda environment for running Prokka
	:param lineage: name of the lineage of interest.
	:param cores: number of cores to use in multiprocessing Prokka cmds.
	:param locus_tag_length: length of locus tags to generate using unique character combinations.
	:param logObject: python logging object handler.
	:param dry_run_flag: flag which indicates commands should only be written to task file and not run. This can be used to manually parallelize across an HPC, if available.
	:param skip_annotation_flag: flag which indicates if Prokka should be run on fast mode with many annotation bells and whistles skipped.
	"""
	task_file = prokka_outdir + 'prokka.cmds'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing Prokka commands to following: %s" % task_file)
		if os.path.isfile(task_file): os.system('rm -f %s' % task_file)
	alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
	possible_locustags = list(itertools.product(alphabet, repeat=locus_tag_length))
	prokka_cmds = []
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

	if not dry_run_flag:
		p = multiprocessing.Pool(cores)
		p.map(util.multiProcess, prokka_cmds)
		p.close()

def runAntiSMASHFromAssemblies(sample_assemblies, antismash_outdir, antismash_load_code, dry_run_flag, cores, logObject,
							   barebone=True):

	"""
	Function to run AntiSMASH biosynthetic gene cluster annotation on assemblies.



	"""
	task_file = antismash_outdir + 'antismash.cmds'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing AntiSMASH commands to following: %s" % task_file)
		if os.path.isfile(task_file): os.system('rm -f %s' % task_file)

	asm_cores = 1
	pool_size = 1
	if cores < 4:
		asm_cores = cores
	else:
		asm_cores = 4
		pool_size = int(cores / 4)
	antismash_cmds = []
	for sample in sample_assemblies:
		sample_resdir = antismash_outdir + sample + '/'
		antismash_cmd = []
		if barebone:
			antismash_cmd = [antismash_load_code, 'antismash', '--taxon', 'bacteria', '--genefinding-tool',
							 'prodigal-m', '--output-dir', sample_resdir, '-c', str(asm_cores),
							 sample_assemblies[sample]]
		else:
			antismash_cmd = [antismash_load_code, 'antismash', '--taxon', 'bacteria', '--genefinding-tool', 'prodigal',
							 '--output-dir', sample_resdir, '--fullhmmer', '--asf', '--cb-general', '--cb-subclusters',
							 '--cb-knownclusters', '--cf-create-clusters', '-c', str(asm_cores),
							 sample_assemblies[sample]]
		antismash_cmds.append(antismash_cmd + [logObject])
		if dry_run_flag:
			task_handle = open(task_file, 'a+')
			task_handle.write(' '.join(antismash_cmd) + '\n')
			task_handle.close()

	if not dry_run_flag:
		p = multiprocessing.Pool(pool_size)
		p.map(multiProcess, antismash_cmds)
		p.close()


def runAntiSMASH(prokka_genbanks_dir, antismash_outdir, antismash_load_code, dry_run_flag, cores, logObject):
	task_file = antismash_outdir + 'antismash.cmds'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing AntiSMASH commands to following: %s" % task_file)
		if os.path.isfile(task_file): os.system('rm -f %s' % task_file)

	asm_cores = 1
	pool_size = 1
	if cores < 4:
		asm_cores = cores
	else:
		asm_cores = 4
		pool_size = int(cores / 4)
	antismash_cmds = []
	for sample_gbk in os.listdir(prokka_genbanks_dir):
		sample = sample_gbk.split('.gbk')[0]
		sample_resdir = antismash_outdir + sample + '/'
		antismash_cmd = [antismash_load_code, 'antismash', '--taxon', 'bacteria', '--genefinding-tool', 'prodigal',
						 '--output-dir', sample_resdir, '--fullhmmer', '--asf', '--cb-general', '--cb-subclusters',
						 '--cb-knownclusters', '--cf-create-clusters', '-c', str(asm_cores),
						 prokka_genbanks_dir + sample_gbk]
		antismash_cmds.append(antismash_cmd + [logObject])
		if dry_run_flag:
			task_handle = open(task_file, 'a+')
			task_handle.write(' '.join(antismash_cmd) + '\n')
			task_handle.close()

	if not dry_run_flag:
		p = multiprocessing.Pool(pool_size)
		p.map(multiProcess, antismash_cmds)
		p.close()

def runOrthoFinder(prokka_proteomes_dir, orthofinder_outdir, orthofinder_load_code, dry_run_flag, cores, logObject):
	task_file = orthofinder_outdir + 'orthofinder.cmd'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing OrthoFinder command to following: %s" % task_file)

	orthofinder_cmd = [orthofinder_load_code, 'orthofinder', '-f', prokka_proteomes_dir, '-t', str(cores)]
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
				[prokka_proteomes_dir + f for f in os.listdir(prokka_proteomes_dir) if f.startswith('Results')][
					0]) + '/'
			os.system('mv %s %s' % (tmp_orthofinder_dir, orthofinder_outdir))
		except:
			logObject.warning('Had an issue running: %s' % ' '.join(orthofinder_cmd))
			raise RuntimeError('Had an issue running: %s' % ' '.join(orthofinder_cmd))


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
	except:
		logObject.error("Had problems extracting protein sequences of CDS features in genbank %s" % bgc_genbank)
		raise RuntimeError("Had problems extracting protein sequences of CDS features in genbank %s" % bgc_genbank)
