import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter

FLANK_SIZE = 500

class BGC:
  def __init__(self, bgc_genbank, bgc_id, comprehensive_parsing=True):
	self.bgc_genbank = bgc_genbank
	self.bgc_id = bgc_id
	self.gene_information = None
	self.cluster_information = None
	self.parseGenbanks(comprehensive_parsing=comprehensive_parsing)

	def parseGenbanks(self, comprehensive_parsing=True):
		"""
		Function to parse an AntiSMASH Genbank file. Values of the object
		"""
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
						core_positions = core_positions.union(set(range(core_start, core_end + 1)))

		if len(bgc_info) == 0:
			bgc_info = [{'detection_rule': 'NA', 'product': 'NA', 'contig_edge': 'NA', 'full_sequence': full_sequence}]

		sys.stderr.write('Processing %s\n' % self.bgc_genbank)
		genes = {}
		core_genes = set([])
		gene_order = {}
		with open(self.bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == "CDS":
						lt = feature.qualifiers.get('locus_tag')[0]
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						direction = str(feature.location).split('(')[1].split(')')[0]

						try:
							product = feature.qualifiers.get('product')[0]
						except:
							product = "hypothetical protein"

						grange = set(range(start, end + 1))
						core_overlap = False
						if len(grange.intersection(core_positions)) > 0:
							core_overlap = True
							core_genes.add(lt)

						gene_order[lt] = start

						prot_seq, nucl_seq, nucl_seq_with_flanks, relative_start, relative_end, gene_domains = [None]*6
						if comprehensive_parsing:
							prot_seq = feature.qualifiers.get('translation')[0]
							gene_domains = []
							for d in domains:
								drange = set(range(d['start'], d['end'] + 1))
								if len(drange.intersection(grange)) > 0:
									gene_domains.append(d)

							flank_start = start - FLANK_SIZE
							flank_end = end + FLANK_SIZE
							if flank_start < 0: flank_start = 0
							if flank_end >= len(full_sequence): flank_end = None
							if end >= len(full_sequence): end = None
							if end:
								nucl_seq = full_sequence[start:end]
							else:
								nucl_seq = full_sequence[start:]
								end = len(full_sequence)
							if flank_end:
								nucl_seq_with_flanks = full_sequence[flank_start:flank_end]
							else:
								nucl_seq_with_flanks = full_sequence[flank_start:]

							gene_length = end - start

							relative_start = start - flank_start
							relative_end = relative_start + gene_length

							if direction == '-':
								nucl_seq = str(Seq(nucl_seq).reverse_complement())
								nucl_seq_with_flanks = str(Seq(nucl_seq_with_flanks).reverse_complement())
								relative_end = len(nucl_seq_with_flanks) - relative_start
								relative_start = relative_end - gene_length

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



