import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import itemgetter
import traceback
import copy

class BGC:
	def __init__(self, bgc_genbank, bgc_id, prediction_method='ANTISMASH'):
		self.bgc_genbank = bgc_genbank
		self.bgc_id = bgc_id
		self.prediction_method = prediction_method.upper()
		self.gene_information = None
		self.cluster_information = None

	def parseGECCO(self, comprehensive_parsing=True, flank_size=2000):
		domains = []
		full_sequence = ""
		domain_evalues = {}
		with open(self.bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				full_sequence = str(rec.seq)
				for feature in rec.features:
					if feature.type == 'misc_feature':
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						asDomain = "NA"
						description = "NA"
						evalue = 1e10
						try:
							aSDomain = feature.qualifiers.get('standard_name')[0]
						except:
							pass
						try:
							description = feature.qualifiers.get('function')[0]
						except:
							pass
						try:
							evalue = [float(x.split(': ')[1].strip()) for x in feature.qualifiers.get('note') if x.startswith('e-value')][0]
						except:
							pass
						domain_evalues[asDomain + '|' + str(start) + '|' + str(end)] = evalue
						domains.append({'start': start + 1, 'end': end, 'type': feature.type, 'aSDomain': aSDomain, 'description': description})

		product = "NA"
		with open(self.bgc_genbank) as ogbk:
			for line in ogbk:
				line = line.strip()
				if line.startswith('biosyn_class'):
					ls = line.split()
					product = ' '.join(ls[2:])
		bgc_info = [{'prediction_method': self.prediction_method, 'detection_rule': 'NA', 'product': product, 'contig_edge': 'NA', 'full_sequence': full_sequence}]

		# determine top 10% of domains with lowest e-values
		num_total_domains = len(domain_evalues)
		core_domains = set([])
		for i, d in enumerate(sorted(domain_evalues.items(), key=itemgetter(1))):
			if i <= num_total_domains*0.1:
				core_domains.add(d[0])

		# sys.stderr.write('Processing %s\n' % self.bgc_genbank)
		genes = {}
		core_genes = set([])
		gene_order = {}
		with open(self.bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == "CDS":
						lt = feature.qualifiers.get('locus_tag')[0]
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						direction = str(feature.location).split('(')[1].split(')')[0]

						try:
							product = feature.qualifiers.get('product')[0]
						except:
							product = "hypothetical protein"

						grange = set(range(start, end + 1))

						gene_domains = []
						core_overlap = False
						for d in domains:
							drange = set(range(d['start'], d['end'] + 1))
							if len(drange.intersection(grange)) > 0:
								gene_domains.append(d)
								if (d['asDomain'] + '|' + str(d['start']) + '|' + str(d['end'])) in core_domains:
									core_overlap = True
									core_genes.add(lt)

						gene_order[lt] = start

						prot_seq, nucl_seq, nucl_seq_with_flanks, relative_start, relative_end, gene_domains = [None] * 6
						if comprehensive_parsing:
							prot_seq = feature.qualifiers.get('translation')[0]

							flank_start = start - flank_size
							flank_end = end + flank_size

							if flank_start < 1: flank_start = 1

							if flank_end >= len(full_sequence): flank_end = None
							if end >= len(full_sequence): end = None

							if end:
								nucl_seq = full_sequence[start - 1:end]
							else:
								nucl_seq = full_sequence[start - 1:]
								end = len(full_sequence)

							if flank_end:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:flank_end]
							else:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:]

							gene_length = end - start

							relative_start = nucl_seq_with_flanks.find(nucl_seq)
							relative_end = relative_start + gene_length

							if direction == '-':
								nucl_seq = str(Seq(nucl_seq).reverse_complement())
								nucl_seq_with_flanks = str(Seq(nucl_seq_with_flanks).reverse_complement())
								relative_start = nucl_seq_with_flanks.find(nucl_seq)
								relative_end = relative_start + gene_length

						genes[lt] = {'bgc_name': self.bgc_id, 'start': start, 'end': end, 'direction': direction,
									 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
									 'nucl_seq_with_flanks': nucl_seq_with_flanks, 'gene_domains': gene_domains,
									 'core_overlap': core_overlap, 'relative_start': relative_start,
									 'relative_end': relative_end}

		number_of_core_gene_groups = 0
		tmp = []
		for lt in sorted(gene_order.items(), key=itemgetter(1), reverse=True):
			if lt[0] in core_genes:
				tmp.append(lt[0])
			elif len(tmp) > 0:
				number_of_core_gene_groups += 1
				tmp = []
		if len(tmp) > 0:
			number_of_core_gene_groups += 1

		for i, pc in enumerate(bgc_info):
			bgc_info[i]['count_core_gene_groups'] = number_of_core_gene_groups

		self.gene_information = genes
		self.cluster_information = bgc_info


	def parseDeepBGC(self, comprehensive_parsing=True, flank_size=2000):
		domains = []
		full_sequence = ""
		domain_score = {}
		product = "NA|NA"
		with open(self.bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				full_sequence = str(rec.seq)
				for feature in rec.features:
					if feature.type == 'PFAM_domain':
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						asDomain = "NA"
						description = "NA"
						deepbgc_score = 0.0
						try:
							aSDomain = feature.qualifiers.get('db_xref')[0]
						except:
							pass
						try:
							description = feature.qualifiers.get('description')[0]
						except:
							pass
						try:
							deepbgc_score = float(feature.qualifiers.get('deepbgc_score')[0])
						except:
							pass
						domain_score[asDomain + '|' + str(start) + '|' + str(end)] = deepbgc_score
						domains.append({'start': start + 1, 'end': end, 'type': feature.type, 'aSDomain': aSDomain, 'description': description})
					elif feature.type == 'cluster':
						product_class = "NA"
						product_activity = "NA"
						try:
							product_activity = feature.qualifiers.get('product_activity')[0]
						except:
							pass
						try:
							product_classes = []
							product_class_scores = feature.qualifiers.get('product_class_score')[0]
							for pci in product_class_scores.split(','):
								pc, pcs = pci.split('=')
								if float(pcs) >= 0.5:
									product_classes.append(pc)
							product_class = '-'.join(sorted(product_classes))
						except:
							pass
						product = 'PC:' + product_class + '_PA:' + product_activity

		bgc_info = [{'prediction_method': self.prediction_method, 'detection_rule': 'NA', 'product': product, 'contig_edge': 'NA', 'full_sequence': full_sequence}]

		# determine top 10% of domains with lowest e-values
		num_total_domains = len(domain_score)
		core_domains = set([])
		for i, d in enumerate(sorted(domain_score.items(), key=itemgetter(1), reverse=True)):
			if i <= num_total_domains*0.1:
				core_domains.add(d[0])

		# sys.stderr.write('Processing %s\n' % self.bgc_genbank)
		genes = {}
		core_genes = set([])
		gene_order = {}
		with open(self.bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == "CDS":
						lt = feature.qualifiers.get('locus_tag')[0]
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						direction = str(feature.location).split('(')[1].split(')')[0]

						try:
							product = feature.qualifiers.get('product')[0]
						except:
							product = "hypothetical protein"

						grange = set(range(start, end + 1))

						gene_domains = []
						core_overlap = False
						for d in domains:
							drange = set(range(d['start'], d['end'] + 1))
							if len(drange.intersection(grange)) > 0:
								gene_domains.append(d)
								if (d['asDomain'] + '|' + str(d['start']) + '|' + str(d['end'])) in core_domains:
									core_overlap = True
									core_genes.add(lt)

						gene_order[lt] = start

						prot_seq, nucl_seq, nucl_seq_with_flanks, relative_start, relative_end, gene_domains = [None] * 6
						if comprehensive_parsing:
							prot_seq = feature.qualifiers.get('translation')[0]

							flank_start = start - flank_size
							flank_end = end + flank_size

							if flank_start < 1: flank_start = 1

							if flank_end >= len(full_sequence): flank_end = None
							if end >= len(full_sequence): end = None

							if end:
								nucl_seq = full_sequence[start - 1:end]
							else:
								nucl_seq = full_sequence[start - 1:]
								end = len(full_sequence)

							if flank_end:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:flank_end]
							else:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:]

							gene_length = end - start

							relative_start = nucl_seq_with_flanks.find(nucl_seq)
							relative_end = relative_start + gene_length

							if direction == '-':
								nucl_seq = str(Seq(nucl_seq).reverse_complement())
								nucl_seq_with_flanks = str(Seq(nucl_seq_with_flanks).reverse_complement())
								relative_start = nucl_seq_with_flanks.find(nucl_seq)
								relative_end = relative_start + gene_length

						genes[lt] = {'bgc_name': self.bgc_id, 'start': start, 'end': end, 'direction': direction,
									 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
									 'nucl_seq_with_flanks': nucl_seq_with_flanks, 'gene_domains': gene_domains,
									 'core_overlap': core_overlap, 'relative_start': relative_start,
									 'relative_end': relative_end}

		number_of_core_gene_groups = 0
		tmp = []
		for lt in sorted(gene_order.items(), key=itemgetter(1), reverse=True):
			if lt[0] in core_genes:
				tmp.append(lt[0])
			elif len(tmp) > 0:
				number_of_core_gene_groups += 1
				tmp = []
		if len(tmp) > 0:
			number_of_core_gene_groups += 1

		for i, pc in enumerate(bgc_info):
			bgc_info[i]['count_core_gene_groups'] = number_of_core_gene_groups

		self.gene_information = genes
		self.cluster_information = bgc_info

	def parseAntiSMASH(self, comprehensive_parsing=True, flank_size=2000):
		bgc_info = []
		domains = []
		core_positions = set([])
		full_sequence = ""
		with open(self.bgc_genbank) as ogbk:
			domain_feature_types = ['PFAM_domain', 'CDS_motif', 'aSDomain']
			for rec in SeqIO.parse(ogbk, 'genbank'):
				full_sequence = str(rec.seq)
				for feature in rec.features:
					if comprehensive_parsing and feature.type in domain_feature_types:
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
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
						domains.append({'start': start + 1, 'end': end, 'type': feature.type, 'aSDomain': aSDomain,
										'description': description})
					elif feature.type == 'protocluster':
						detection_rule = feature.qualifiers.get('detection_rule')[0]
						try:
							product = feature.qualifiers.get('product')[0]
						except:
							product = "NA"
						contig_edge = feature.qualifiers.get('contig_edge')[0]
						bgc_info.append(
							{'detection_rule': detection_rule, 'product': product, 'contig_edge': contig_edge,
							 'full_sequence': str(rec.seq)})
					elif feature.type == 'proto_core':
						core_start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						core_end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						core_positions = core_positions.union(set(range(core_start + 1, core_end + 1)))

		if len(bgc_info) == 0:
			bgc_info = [
				{'detection_rule': 'NA', 'product': 'NA', 'contig_edge': 'NA', 'full_sequence': full_sequence}]

		# sys.stderr.write('Processing %s\n' % self.bgc_genbank)
		genes = {}
		core_genes = set([])
		gene_order = {}
		with open(self.bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == "CDS":
						lt = feature.qualifiers.get('locus_tag')[0]
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						direction = str(feature.location).split('(')[1].split(')')[0]

						try:
							product = feature.qualifiers.get('product')[0]
						except:
							product = "hypothetical protein"

						rule_based_bgc_cds = False
						try:
							if 'rule-based-clusters' in feature.qualifiers.get('gene_functions')[0]:
								rule_based_bgc_cds = True
						except:
							pass

						grange = set(range(start, end + 1))
						core_overlap = False
						if len(grange.intersection(core_positions)) > 0 and rule_based_bgc_cds:
							core_overlap = True
							core_genes.add(lt)

						gene_order[lt] = start

						prot_seq, nucl_seq, nucl_seq_with_flanks, relative_start, relative_end, gene_domains = [None] * 6
						if comprehensive_parsing:
							prot_seq = feature.qualifiers.get('translation')[0]
							gene_domains = []
							for d in domains:
								drange = set(range(d['start'], d['end'] + 1))
								if len(drange.intersection(grange)) > 0:
									gene_domains.append(d)

							flank_start = start - flank_size
							flank_end = end + flank_size

							if flank_start < 1: flank_start = 1

							if flank_end >= len(full_sequence): flank_end = None
							if end >= len(full_sequence): end = None

							if end:
								nucl_seq = full_sequence[start - 1:end]
							else:
								nucl_seq = full_sequence[start - 1:]
								end = len(full_sequence)

							if flank_end:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:flank_end]
							else:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:]

							gene_length = end - start

							relative_start = nucl_seq_with_flanks.find(nucl_seq)
							relative_end = relative_start + gene_length

							if direction == '-':
								nucl_seq = str(Seq(nucl_seq).reverse_complement())
								nucl_seq_with_flanks = str(Seq(nucl_seq_with_flanks).reverse_complement())
								relative_start = nucl_seq_with_flanks.find(nucl_seq)
								relative_end = relative_start + gene_length

						genes[lt] = {'bgc_name': self.bgc_id, 'start': start, 'end': end, 'direction': direction,
									 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
									 'nucl_seq_with_flanks': nucl_seq_with_flanks, 'gene_domains': gene_domains,
									 'core_overlap': core_overlap, 'relative_start': relative_start,
									 'relative_end': relative_end}

		number_of_core_gene_groups = 0
		tmp = []
		for lt in sorted(gene_order.items(), key=itemgetter(1), reverse=True):
			if lt[0] in core_genes:
				tmp.append(lt[0])
			elif len(tmp) > 0:
				number_of_core_gene_groups += 1
				tmp = []
		if len(tmp) > 0:
			number_of_core_gene_groups += 1

		for i, pc in enumerate(bgc_info):
			bgc_info[i]['count_core_gene_groups'] = number_of_core_gene_groups

		self.gene_information = genes
		self.cluster_information = bgc_info

	def parseGenbanks(self, comprehensive_parsing=True, flank_size=2000):
		"""
		Function to parse an AntiSMASH, DeepBGC, or GECCO Genbank file.
		"""
		if self.prediction_method.upper() == 'ANTISMASH':
			self.parseAntiSMASH(comprehensive_parsing=comprehensive_parsing, flank_size=flank_size)
		elif self.prediction_method.upper() == 'DEEPBGC':
			self.parseDeepBGC(comprehensive_parsing=comprehensive_parsing, flank_size=flank_size)
		elif self.prediction_method.upper() == 'GECCO':
			self.parseGECCO(comprehensive_parsing=comprehensive_parsing, flank_size=flank_size)
		else:
			raise RuntimeError("Unable to parse file because BGC prediction method is not an accepted option!")

	def refineGenbank(self, refined_genbank_file, first_bg, second_bg):
		"""
		Function to prune and update coordinates of BGC Genbank
		"""
		try:
			rgf_handle = open(refined_genbank_file, 'w')
			start_coord = min([self.gene_information[first_bg]['start'], self.gene_information[first_bg]['end'], self.gene_information[second_bg]['start'], self.gene_information[second_bg]['end']])
			end_coord = max([self.gene_information[first_bg]['start'], self.gene_information[first_bg]['end'], self.gene_information[second_bg]['start'], self.gene_information[second_bg]['end']])
			pruned_coords = set(range(start_coord, end_coord+1))
			with open(self.bgc_genbank) as ogbk:
				recs = [x for x in SeqIO.parse(ogbk, 'genbank')]
				try:
					assert(len(recs) == 1)
				except Exception as e:
					raise RuntimeError(traceback.format_exc())
				for rec in recs:
					original_seq = str(rec.seq)
					filtered_seq = ""
					if end_coord == len(original_seq):
						filtered_seq = original_seq[start_coord-1:]
					else:
						filtered_seq = original_seq[start_coord-1:end_coord]

					new_seq_object = Seq(filtered_seq)

					updated_rec = copy.deepcopy(rec)
					updated_rec.seq = new_seq_object

					updated_features = []
					for feature in rec.features:
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])+1
						end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])

						feature_coords = set(range(start, end+1))
						if len(feature_coords.intersection(pruned_coords)) > 0:
							updated_start = start - start_coord + 1
							updated_end = end - start_coord + 1

							if end > end_coord:
								if feature.type == 'CDS':
									continue
								else:
									updated_end = end_coord - start_coord + 1
							if start < start_coord:
								if feature.type == 'CDS':
									continue
								else:
									updated_start = 1

							strand = 1
							if '(-)' in str(feature.location):
								strand = -1

							updated_location = FeatureLocation(updated_start-1, updated_end, strand=strand)
							updated_feature = copy.deepcopy(feature)

							updated_feature.location = updated_location
							updated_features.append(updated_feature)
					updated_rec.features = updated_features
					SeqIO.write(updated_rec, rgf_handle, 'genbank')
			rgf_handle.close()
		except Exception as e:
			raise RuntimeError(traceback.format_exc())

