import os
import sys
from Bio import SeqIO
import logging
from logging.handlers import QueueHandler, QueueListener
import subprocess
import multiprocessing
from collections import defaultdict
from operator import itemgetter
import itertools

def multiProcess(input):
	input_cmd = input[:-1]
	logObject = input[-1]
	logObject.info('Running the following command: %s' % ' '.join(input_cmd))
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=sys.stderr, executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(input_cmd))
	except:
		logObject.warning('Had an issue running: %s' % ' '.join(input_cmd))

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
			subprocess.call(' '.join(orthofinder_cmd), shell=True, stdout=sys.stderr, stderr=sys.stderr)
			logObject.info('Successfully ran OrthoFinder!')
			tmp_orthofinder_dir = os.path.abspath([prokka_proteomes_dir + f for f in os.listdir(prokka_proteomes_dir) if f.startswith('Results')][0]) + '/'
			os.system('mv %s %s' % (tmp_orthofinder_dir, orthofinder_outdir))
		except:
			logObject.warning('Had an issue running: %s' % ' '.join(orthofinder_cmd))
			raise RuntimeError('Had an issue running: %s' % ' '.join(orthofinder_cmd))

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
		pool_size = int(cores/4)
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

def runProkka(sample_assemblies, prokka_outdir, prokka_proteomes, prokka_genbanks, prokka_load_code, dry_run_flag, lineage, cores, logObject):
	task_file = prokka_outdir + 'prokka.cmds'
	if dry_run_flag:
		logObject.info("Dry run on: will just be writing Prokka commands to following: %s" % task_file)
		if os.path.isfile(task_file): os.system('rm -f %s' % task_file)
	alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
	possible_locustags = list(itertools.product(alphabet, repeat=3))
	prokka_cmds = []
	for i, sample in enumerate(sample_assemblies):
		sample_assembly = sample_assemblies[sample]
		sample_outdir = prokka_outdir + sample + '/'
		prokka_cmd = [prokka_load_code, 'prokka', '--cpus', '1', '--outdir', sample_outdir, '--prefix', sample, '--genus', lineage, '--locustag', ''.join(list(possible_locustags[i])),
						sample_assembly, ';', 'mv', sample_outdir + sample + '.gbf', prokka_genbanks + sample + '.gbk', ';', 'mv',
						sample_outdir + sample + '.faa', prokka_proteomes]
		prokka_cmds.append(prokka_cmd + [logObject])
		if dry_run_flag:
			task_handle = open(task_file, 'a+')
			task_handle.write(' '.join(prokka_cmd) + '\n')
			task_handle.close()

	if not dry_run_flag:
		p = multiprocessing.Pool(cores)
		p.map(multiProcess, prokka_cmds)
		p.close()

def is_fasta(fasta):
	try:
		with open(fasta) as of:
			SeqIO.parse(of, 'fasta')
		return True
	except:
		return False

def createLoggerObject(log_file):
	"""
	Function which creates logger object.
	:param log_file: path to log file.
	:return: logging logger object.
	"""

	logger = logging.getLogger('task_logger')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	# create console handler with a higher log level
	ch = logging.StreamHandler()
	ch.setLevel(logging.ERROR)
	# create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)
	# add the handlers to the logger
	logger.addHandler(fh)
	logger.addHandler(ch)

	#q = multiprocessing.Queue()
	#queue_listner = QueueListener(q, ch)
	#queue_listner.start()

	#logger.handlers[0].stream = sys.stderr
	return logger

def readInAssemblyListing(assembly_listing_file, logObject):
	sample_assembly_paths = {}
	try:
		with open(assembly_listing_file) as oalf:
			for line in oalf:
				line = line.strip()
				sample, assembly = line.split('\t')
				sample = sample.replace(' ', '_').replace('|', '_').replace('"', '_').replace("'", '_')
				try:
					assert (is_fasta(assembly))
					sample_assembly_paths[sample] = assembly
				except:
					logObject.info('Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.' % sample)
					sys.stderr.write('Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.\n' % sample)
		assert (len(sample_assembly_paths) >= 2)
		return (sample_assembly_paths)
	except:
		raise RuntimeError('Input file listing the location of assemblies for samples leads to incorrect path was provided. Exiting now ...')

def closeLoggerObject(logObject):
	"""
	Function which closes/terminates loggerObject.
	:param logObject: logging logger object to close
	"""
	handlers = logObject.handlers[:]
	for handler in handlers:
		handler.close()
		logObject.removeHandler(handler)

def logParameters(parameter_names, parameter_values):
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		sys.stderr.write(pn + ': ' + str(pv) + '\n')

def logParametersToFile(parameter_file, parameter_names, parameter_values):
	parameter_handle = open(parameter_file, 'w')
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		parameter_handle.write(pn + ': ' + str(pv) + '\n')
	parameter_handle.close()

def logParametersToObject(logObject, parameter_names, parameter_values):
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		logObject.info(pn + ': ' + str(pv))