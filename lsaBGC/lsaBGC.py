import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from ete3 import Tree
import logging
import subprocess
import multiprocessing
from collections import defaultdict
from operator import itemgetter
import itertools
import statistics
import random

def multiProcess(input):
	input_cmd = input[:-1]
	logObject = input[-1]
	logObject.info('Running the following command: %s' % ' '.join(input_cmd))
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(input_cmd))
	except:
		logObject.warning('Had an issue running: %s' % ' '.join(input_cmd))

def readInBGCGenbanksPerGCF(gcf_specs_file, logObject):
	sample_index = defaultdict(int)
	bgc_gbk = {}
	bgc_sample = {}
	sample_bgcs = defaultdict(set)
	bgc_genes = {}
	all_genes = set([])
	comp_gene_info = {}
	with open(gcf_specs_file) as obsf:
		for i, line in enumerate(obsf):
			line = line.strip()
			try:
				assert (len(line.split('\t')) == 2)
			except:
				logObject.error("More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
				raise RuntimeError("More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
			sample, gbk = line.split('\t')
			try:
				assert (is_genbank(gbk))
				bgc_id = sample
				if sample_index[sample] > 0:
					bgc_id = sample + '_' + str(sample_index[sample] + 1)
				sample_index[sample] += 1

				bgc_genes_full, bgc_info_full = parseGenbanks(gbk, bgc_id)

				comp_gene_info.update(bgc_genes_full)
				all_genes = all_genes.union(bgc_genes[gbk])
				bgc_gbk[bgc_id] = gbk
				bgc_genes[bgc_id] = set(bgc_genes_full.keys())
				bgc_sample[bgc_id] = sample
				sample_bgcs[sample].add(bgc_id)

				logObject.info("Incorporating genbank %s for sample %s into analysis." % (gbk, sample))
			except:
				logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis.")
				raise RuntimeWarning("Unable to validate %s as Genbank. Skipping ...")

	return [bgc_gbk, bgc_genes, comp_gene_info, all_genes, bgc_sample, sample_bgcs]

def createItolBGCSeeTrack(bgc_genes, gene_to_cog, cog_to_color, comp_gene_info, dataset_label, outdir, logObject):
	track_file = outdir + 'BGCs_Visualization.iTol.txt'
	track_handle = open(track_file, 'w')

	logObject.info("Writing track file to: %s" % track_file)
	logObject.info("Track will have label: %s" % dataset_label)

	# write header for iTol track file
	track_handle.write('DATASET_DOMAINS\n')
	track_handle.write('SEPARATOR TAB\n')
	track_handle.write('DATASET_LABEL\t%s\n' % dataset_label)
	track_handle.write('COLOR\t#000000\n')
	track_handle.write('BORDER_WIDTH\t1\n')
	track_handle.write('BORDER_COLOR\t#000000\n')
	track_handle.write('SHOW_DOMAIN_LABELS\t0\n')
	track_handle.write('DATA\n')

	# write the rest of the iTol track file for illustrating genes across BGC instances
	ref_cog_directions = {}
	for i, bgc in enumerate(bgc_genes):
		curr_bgc_genes = bgc_genes[bgc]
		#genes[lt] = {'bgc_name': bgc_name, 'start': start, 'end': end, 'direction': direction,
		#			 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
		#			 'gene_domains': gene_domains, 'core_overlap': core_overlap}

		last_gene_end = max([comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
		printlist = [bgc, str(last_gene_end)]
		cog_directions = {}
		for lt in curr_bgc_genes:
			ginfo = comp_gene_info[lt]
			cog = 'singleton'
			if lt in gene_to_cog:
				cog = gene_to_cog[lt]
			shape = 'None'
			if ginfo['direction'] == '+':
				shape = 'TR'
			elif ginfo['direction'] == '-':
				shape = 'TL'
			gstart = ginfo['start']
			gend = ginfo['end']
			cog_color = "#dbdbdb"
			if cog in cog_to_color:
				cog_color = cog_to_color[cog]
			gene_string = '|'.join([str(x) for x in [shape, gstart, gend, cog_color, cog]])
			printlist.append(gene_string)
			if cog != 'singleton':
				cog_directions[cog] = ginfo['direction']
		if i == 0:
			ref_cog_directions = cog_directions
			track_handle.write('\t'.join(printlist) + '\n')
		else:
			flip_support = 0
			keep_support = 0
			for c in ref_cog_directions:
				if not c in cog_directions: continue
				if cog_directions[c] == ref_cog_directions[c]:
					keep_support += 1
				else:
					flip_support += 1

			# flip the genbank visual if necessary, first BGC processed is used as reference guide
			if flip_support > keep_support:
				flip_printlist = printlist[:2]
				bgc_end = int(printlist[1])
				for gene_string in printlist[2:]:
					gene_info = gene_string.split('|')
					new_shape = None
					if gene_info[0] == 'TR':
						new_shape = 'TL'
					elif gene_info[0] == 'TL':
						new_shape = 'TR'
					new_gstart = int(bgc_end) - int(gene_info[2])
					new_gend = int(bgc_end) - int(gene_info[1])
					new_gene_info = '|'.join([new_shape, str(new_gstart), str(new_gend)] + gene_info[-2:])
					flip_printlist.append(new_gene_info)
				track_handle.write('\t'.join(flip_printlist) + '\n')
			else:
				track_handle.write('\t'.join(printlist) + '\n')
	track_handle.close()

def constructCodonAlignments(bgc_sample, cog_genes, comp_gene_info, outdir, cores, logObject):
	nucl_seq_dir = os.path.abspath(outdir + 'Nucleotide_Sequences') + '/'
	prot_seq_dir = os.path.abspath(outdir + 'Protein_Sequences') + '/'
	prot_alg_dir = os.path.abspath(outdir + 'Protein_Alignments') + '/'
	codo_alg_dir = os.path.abspath(outdir + 'Codon_Alignments') + '/'
	if not os.path.isdir(nucl_seq_dir): os.system('mkdir %s' % nucl_seq_dir)
	if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
	if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
	if not os.path.isdir(codo_alg_dir): os.system('mkdir %s' % codo_alg_dir)

	all_samples = set([bgc_sample.values()])
	try:
		inputs = []
		for cog in cog_genes:
			sample_counts = defaultdict(int)
			gene_data = {}
			for gene in cog:
				gene_info = comp_gene_info[gene]
				bgc_id = gene_info['bgc_name']
				sample_id = bgc_sample[bgc_id]
				sample_counts[sample_id] += 1
			samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
			# check that cog is single-copy-core
			if len(samples_with_single_copy.symmetric_difference(all_samples)) > 0: continue
			logObject.info
			cog, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir])

		p = multiprocessing.Pool(cores)
		p.map(create_msas, )
	except:

def create_msas(input):
	cog, = input
	cog_nucl_fasta = nucl_seq_dir + '/' + cog + '.fna'
	cog_prot_fasta = prot_seq_dir + '/' + cog + '.faa'
	cog_prot_msa = prot_alg_dir + '/' + cog + '.msa.faa'
	cog_codo_msa = codo_alg_dir + '/' + cog + '.msa.fna'

	cog_nucl_handle = open(cog_nucl_fasta, 'w')
	cog_prot_handle = open(cog_prot_fasta, 'w')
	for gene in bgc_cog_genes[cog]:
		cog_nucl_handle.write(
			'>' + comp_gene_info[gene]['bgc_name'] + '|' + gene + '\n' + comp_gene_info[gene]['nucl_seq'] + '\n')
		cog_prot_handle.write(
			'>' + comp_gene_info[gene]['bgc_name'] + '|' + gene + '\n' + comp_gene_info[gene]['prot_seq'] + '\n')
	cog_nucl_handle.close()
	cog_prot_handle.close()

	os.system('mafft --maxiterate 1000 --localpair %s > %s' % (cog_prot_fasta, cog_prot_msa))
	os.system('pal2nal.pl %s %s -output fasta > %s' % (cog_prot_msa, cog_nucl_fasta, cog_codo_msa))
def parseSpeciesPhylogeny(species_phylogeny, sample_bgcs, outdir, logObject):
	try:
		number_of_added_leaves = 0
		t = Tree(species_phylogeny)
		for node in t.traverse('postorder'):
			if node.name in sample_bgcs and len(sample_bgcs[node.name]) > 1:
				for bgc_id in sample_bgcs[node.name]:
					if bgc_id == node.name: continue
					node.add_sister(name=bgc_id)
					sister_node = t.search_nodes(name=bgc_id)[0]
					sister_node.dist = node.dist
					number_of_added_leaves += 1
		edited_species_phylogeny = outdir + 'species.expanded.nwk'
		t.write(format=1, outfile=edited_species_phylogeny)
		logObject.info("New phylogeny with an additional %d leafs to reflect samples with multiple BGCs can be found at: %s." % (number_of_added_leaves, edited_species_phylogeny))
	except:
		logObject.error("Had difficulties properly editing species phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")
		raise RuntimeError("Had difficulties properly editing species phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")

def readInBGCGenbanksComprehensive(bgc_specs_file, logObject):
	bgc_sample = {}
	bgc_product = {}
	bgc_genes = {}
	all_genes = set([])
	with open(bgc_specs_file) as obsf:
		for i, line in enumerate(obsf):
			line = line.strip()
			try:
				assert (len(line.split('\t')) == 2)
			except:
				logObject.error(
					"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
				raise RuntimeError(
					"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
			sample, gbk = line.split('\t')
			try:
				assert (is_genbank(gbk))
				bgc_genes_full, bgc_info_full = parseGenbanks(gbk, gbk)
				bgc_product[gbk] = [x['product'] for x in bgc_info_full]
				bgc_genes[gbk] = set(bgc_genes_full.keys())
				all_genes = all_genes.union(bgc_genes[gbk])
				bgc_sample[gbk] = sample
				logObject.info("Incorporating genbank %s for sample %s into analysis." % (gbk, sample))
			except:
				logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis.")
				raise RuntimeWarning("Unable to validate %s as Genbank. Skipping ...")

	return [bgc_sample, bgc_product, bgc_genes, all_genes]


def parseGenbanks(gbk, bgc_name):
	"""
    :param gbk: AntiSMASH Genbank file
    :return:
    """
	bgc_info = []
	domains = []
	core_positions = set([])
	full_sequence = ""
	with open(gbk) as ogbk:
		domain_feature_types = ['PFAM_domain', 'CDS_motif', 'aSDomain']
		for rec in SeqIO.parse(ogbk, 'genbank'):
			for feature in rec.features:
				if feature.type in domain_feature_types:
					start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					aSDomain = "NA"
					description = "NA"
					try:
						aSDomain = feature.qualifiers.get('aSDomain')[0]
					except:
						pass
					try:
						description = feature.qualifiers.get('description')[0]
					except:
						pass
					domains.append({'start': start, 'end': end, 'type': feature.type, 'aSDomain': aSDomain,
									'description': description})
				elif feature.type == 'protocluster':
					detection_rule = feature.qualifiers.get('detection_rule')[0]
					product = feature.qualifiers.get('product')[0]
					contig_edge = feature.qualifiers.get('contig_edge')[0]
					full_sequence = str(rec.seq)
					bgc_info.append({'detection_rule': detection_rule, 'product': product, 'contig_edge': contig_edge,
									 'full_sequence': str(rec.seq)})
				elif feature.type == 'proto_core':
					core_start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					core_end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					core_positions = core_positions.union(set(range(core_start, core_end + 1)))

	if len(bgc_info) == 0:
		bgc_info = [{'detection_rule': 'NA', 'product': 'NA', 'contig_edge': 'NA', 'full_sequence': full_sequence}]

	genes = {}
	with open(gbk) as ogbk:
		for rec in SeqIO.parse(ogbk, 'genbank'):
			for feature in rec.features:
				if feature.type == "CDS":
					lt = feature.qualifiers.get('locus_tag')[0]
					prot_seq = feature.qualifiers.get('translation')[0]
					start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
					direction = str(feature.location).split('(')[1].split(')')[0]
					product = feature.qualifiers.get('product')[0]
					grange = set(range(start, end + 1))
					gene_domains = []
					for d in domains:
						drange = set(range(d['start'], d['end'] + 1))
						if len(drange.intersection(grange)) > 0:
							gene_domains.append(d)
					nucl_seq = full_sequence[start:end]
					core_overlap = False
					if len(grange.intersection(core_positions)) > 0: core_overlap = True
					if direction == '-':
						nucl_seq = str(Seq(full_sequence[start:end]).reverse_complement())
					genes[lt] = {'bgc_name': bgc_name, 'start': start, 'end': end, 'direction': direction,
								 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
								 'gene_domains': gene_domains, 'core_overlap': core_overlap}
	return ([genes, bgc_info])

def calculateBGCPairwiseRelations(bgc_genes, gene_to_cog, avg_nz_cog_counts, outdir, logObject):
	pair_relations_txt_file = outdir + 'bgc_pair_relationships.txt'
	pair_relations_mci_file = outdir + 'bgc_pair_relationships.mci'
	pair_relations_tab_file = outdir + 'bgc_pair_relationships.tab'

	bgc_cogs = defaultdict(set)
	for bgc in bgc_genes:
		for gene in bgc_genes[bgc]:
			if gene in gene_to_cog:	bgc_cogs[bgc].add(gene_to_cog[gene])

	try:
		prf_handle = open(pair_relations_txt_file, 'w')
		pairwise_relations = defaultdict(lambda: defaultdict(float))
		for i, bgc1 in enumerate(bgc_cogs):
			bgc1_cogs = set([x for x in bgc_cogs[bgc1] if avg_nz_cog_counts[x] < 2])
			for j, bgc2 in enumerate(bgc_cogs):
				if i < j:
					bgc2_cogs = set([x for x in bgc_cogs[bgc2] if avg_nz_cog_counts[x] < 2])
					overlap_metric = float(len(bgc1_cogs.intersection(bgc2_cogs)))/float(min([len(bgc1_cogs), len(bgc2_cogs)]))
					overlap_metric_scaled = 100.00*overlap_metric
					if overlap_metric_scaled > 0:
						pairwise_relations[bgc1][bgc2] = overlap_metric_scaled
						pairwise_relations[bgc2][bgc1] = overlap_metric_scaled
						prf_handle.write('%s\t%s\t%f\n' % (bgc1, bgc2, overlap_metric_scaled))
		prf_handle.close()
		logObject.info("Calculated pairwise relations and wrote to: %s" % pair_relations_txt_file)
	except:
		logObject.error("Problem with creating relations file between pairs of BGCs.")
		raise RuntimeError("Problem with creating relations file between pairs of BGCs.")

	mcxload_cmd = ['mcxload', '-abc', pair_relations_txt_file, '--stream-mirror', '-write-tab', pair_relations_tab_file,
	               '-o', pair_relations_mci_file]

	logObject.info('Running the following command: %s' % ' '.join(mcxload_cmd))
	try:
		subprocess.call(' '.join(mcxload_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(mcxload_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(mcxload_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(mcxload_cmd))
	logObject.info('Converted format of pair relationship file via mxcload.')

	return bgc_cogs, pairwise_relations, pair_relations_txt_file, pair_relations_mci_file, pair_relations_tab_file

def runMCLAndReportGCFs(mip, outdir, sf_handle, pairwise_relations, pair_relations_mci_file, pair_relations_tab_file, bgc_cogs, bgc_product, bgc_sample, inflation_testing, cores, logObject):
	relations_mcl_file = outdir + 'mcl.' + str(mip).replace('.', '_') + '.out'
	mcxdump_out_file = outdir + 'final_mcl.' + str(mip).replace('.', '_') + '.out'
	mcl_cmd = ['mcl', pair_relations_mci_file, '-I', str(mip), '-o', relations_mcl_file, '-te', str(cores)]
	mcxdump_cmd = ['mcxdump', '-icl', relations_mcl_file, '-tabr', pair_relations_tab_file, '-o', mcxdump_out_file]

	logObject.info('Running MCL and MCXDUMP with inflation parameter set to %f' % mip)
	logObject.info('Running the following command: %s' % ' '.join(mcl_cmd))
	try:
		subprocess.call(' '.join(mcl_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
		assert(os.path.isfile(relations_mcl_file) and os.path.getsize(relations_mcl_file) > 100)
		logObject.info('Successfully ran: %s' % ' '.join(mcl_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(mcl_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(mcl_cmd))

	sys.stderr.write('Dumping results in human-readable format ...')
	logObject.info('Running the following command: %s' % ' '.join(mcxdump_cmd))
	try:
		subprocess.call(' '.join(mcxdump_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
		assert(os.path.isfile(mcxdump_out_file) and os.path.getsize(mcxdump_out_file) > 100)
		logObject.info('Successfully ran: %s' % ' '.join(mcl_cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(mcxdump_cmd))
		raise RuntimeError('Had an issue running: %s' % ' '.join(mcxdump_cmd))

	logObject.info('Successfully ran MCL and MCXDUMP with inflation parameter set to %f!' % mip)

	clustered_bgcs = set([])
	with open(mcxdump_out_file) as omo:
		for j, gcf in enumerate(omo):
			gcf = gcf.strip()
			gcf_mems = gcf.split()
			if len(gcf_mems) < 2: continue
			diffs = set([])
			samp_counts = defaultdict(int)
			samp_ogs = defaultdict(set)
			products = set([])
			for a, bgc1 in enumerate(gcf_mems):
				samp_counts[bgc1] += 1
				for prod in bgc_product[bgc1]: products.add(prod)
				samp_ogs[bgc1] = samp_ogs[bgc1].union(bgc_cogs[bgc1])
				clustered_bgcs.add(bgc1)
				for b, bgc2 in enumerate(gcf_mems):
					if a < b:
						diffs.add(pairwise_relations[bgc1][bgc2])
			multi_same_sample = 0
			num_ogs = []
			for si, s in enumerate(samp_counts):
				if samp_counts[s] > 1:
					multi_same_sample += 1
				if si == 0:
					scc = samp_ogs[s]
				else:
					scc = scc.intersection(samp_ogs[s])
				num_ogs.append(len(samp_ogs[s]))
			gcf_stats = ['GCF_' + str(j+1), len(gcf_mems), multi_same_sample, len(scc), statistics.mean(num_ogs),
						 statistics.stdev(num_ogs), min(diffs), max(diffs), '; '.join(products)]
			if inflation_testing:
				gcf_stats = [mip] + gcf_stats
			sf_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
	singleton_bgcs = set([])
	for bgc in bgc_cogs:
		if not bgc in clustered_bgcs: singleton_bgcs.add(bgc)
	singleton_stats = ['singletons', len(singleton_bgcs)] + (['NA'] * 7)
	if inflation_testing:
		singleton_stats = [mip] + singleton_stats
	sf_handle.write('\t'.join([str(x) for x in singleton_stats]) + '\n')

	if not inflation_testing:
		logObject.info("Writing list of BGCs for each GCF, which will be used as input for downstream programs in the suite!")
		gcf_listing_dir = outdir + 'GCF_Listings/'
		if not os.path.isdir(gcf_listing_dir): os.system('mkdir %s' % gcf_listing_dir)
		with open(mcxdump_out_file) as omo:
			for j, gcf in enumerate(omo):
				gcf = gcf.strip()
				gcf_mems = gcf.split()
				if len(gcf_mems) < 2: continue
				outf_list = open(gcf_listing_dir + 'GCF_' + str(j + 1) + '.txt', 'w')
				for bgc in gcf_mems:
					sname = bgc_sample[bgc]
					gbkpath = bgc
					outf_list.write('%s\t%s\n' % (sname, gbkpath))
				outf_list.close()
		logObject.info("Successfully wrote lists of BGCs for each GCF.")
	return

def parseOrthoFinderMatrix(orthofinder_matrix_file, all_gene_lts):
	"""
    :param orthofinder_matrix: OrthoFinderV2 matrix Orthogroups.csv file, should also include singleton orthologroups
    :param all_gene_lts: Set of all the relevant gene locus tag identifiers found in BGC Genbanks
    :return: dictionary mapping gene locus tags to CoGs
    """
	gene_to_cog = {}
	cog_genes = defaultdict(set)
	cog_median_gene_counts = defaultdict(lambda: 'NA')
	with open(orthofinder_matrix_file) as ofm:
		for i, line in enumerate(ofm):
			if i == 0: continue
			line = line.strip()
			ls = line.split('\t')
			cog = ls[0]
			flag_in_bgc = False
			for sgs in ls[1:]:
				for g in sgs.split(', '):
					if g in all_gene_lts:
						flag_in_bgc = True
						gene_to_cog[g] = cog
						cog_genes[cog].add(g)
			if flag_in_bgc:
				gene_counts = []
				for sgs in ls[1:]:
					gene_counts.append(len(sgs.split(', ')))
				cog_median_gene_counts[cog] = statistics.median(gene_counts)
	return ([gene_to_cog, cog_genes, cog_median_gene_counts])

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
					logObject.info(
						'Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.' % sample)
					sys.stderr.write(
						'Ignoring sample %s, because assembly at: does not seem to exist or be a FASTA.\n' % sample)
		assert (len(sample_assembly_paths) >= 2)
		return (sample_assembly_paths)
	except:
		raise RuntimeError('Input file listing the location of assemblies for samples leads to incorrect path was provided. Exiting now ...')

def runProkka(sample_assemblies, prokka_outdir, prokka_proteomes, prokka_genbanks, prokka_load_code, dry_run_flag,
			  lineage, cores, logObject):
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
		prokka_cmd = [prokka_load_code, 'prokka', '--cpus', '1', '--outdir', sample_outdir, '--prefix', sample,
					  '--genus', lineage, '--locustag', ''.join(list(possible_locustags[i])),
					  sample_assembly, ';', 'mv', sample_outdir + sample + '.gbf', prokka_genbanks + sample + '.gbk',
					  ';', 'mv',
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

def runAntiSMASH(prokka_genbanks_dir, antismash_outdir, antismash_load_code, dry_run_flag, cores, logObject, assemblies_draft=False):
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
		genefinding_tool = 'prodigal'
		if assemblies_draft: genefinding_tool = 'prodigal-m'
		antismash_cmd = [antismash_load_code, 'antismash', '--taxon', 'bacteria', '--genefinding-tool', genefinding_tool,
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


def is_fasta(fasta):
	try:
		with open(fasta) as of:
			SeqIO.parse(of, 'fasta')
		return True
	except:
		return False


def is_genbank(fasta):
	try:
		with open(fasta) as of:
			SeqIO.parse(of, 'genbank')
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

	# q = multiprocessing.Queue()
	# queue_listner = QueueListener(q, ch)
	# queue_listner.start()

	# logger.handlers[0].stream = sys.stderr
	return logger


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


def assignColorsToCOGs(cogs):
	"""
    :param cogs: set of CoGs.
    :return: dictionary mapping each CoG to a hex color value.
    """

	# read in list of colors
	dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
	colors_file = dir_path + 'colors_200.txt'
	colors = []
	with open(colors_file) as ocf:
		colors = [x.strip() for x in ocf.readlines()]
	random.shuffle(colors)

	cog_to_color = {}
	for i, c in enumerate(set(cogs)):
		cog_to_color[c] = colors[i]
	return(cog_to_color)

