import os
import sys
from lsaBGC.classes.BGC import BGC
from lsaBGC import util as lBu
from collections import defaultdict

class Pan:
  def __init__(self, bgc_genbanks_listing, lineage_name='Unnamed lineage'):
	self.bgc_genbanks_listing = bgc_genbanks_listing
	self.lineage_name = lineage_name
	self.bgc_comp_info = None

	def readInBGCGenbanks(self, logObject=None, comprehensive_parsing=True):
		"""
		:param logObject (optional): logging object handler.
		:param comprehensive_parsing (optional): flag specifying whether to perform comprehensive extraction of information from Genbanks. default is True.
		:return:

		"""
		sample_index = defaultdict(int)
		bgc_gbk = {}
		bgc_sample = {}
		sample_bgcs = defaultdict(set)
		bgc_genes = {}
		all_genes = set([])
		comp_gene_info = {}
		bgc_info = {}
		with open(self.bgc_genbanks_listing) as obsf:
			for i, line in enumerate(obsf):
				line = line.strip()
				try:
					assert (len(line.split('\t')) == 2)
				except:
					logObject.error(
						"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
								i + 1))
					raise RuntimeError(
						"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
								i + 1))
				sample, gbk = line.split('\t')
				try:
					assert (lBu.is_genbank(gbk))
					bgc_id = sample
					if sample_index[sample] > 0:
						bgc_id = sample + '_' + str(sample_index[sample] + 1)
					sample_index[sample] += 1

					BGC_Object = BGC(gbk, gbk, comprehensive_parsing=comprehensive_parsing)
					comp_gene_info.update(BGC_Object.gene_information)
					bgc_info[bgc_id] = BGC_Object.cluster_information
					bgc_genes[bgc_id] = set(BGC_Object.gene_information)
					all_genes = all_genes.union(BGC_Object.gene_information.keys())
					bgc_gbk[bgc_id] = gbk
					bgc_sample[bgc_id] = sample
					sample_bgcs[sample].add(bgc_id)

					logObject.info("Incorporating genbank %s for sample %s into analysis." % (gbk, sample))
				except:
					logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis.")
					raise RuntimeWarning("Unable to validate %s as Genbank. Skipping ...")
