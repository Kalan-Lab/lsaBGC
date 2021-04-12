import os
import sys
import logging
import traceback
import statistics
import random
from ete3 import Tree
from operator import itemgetter
from collections import defaultdict
from lsaBGC.classes.Pan import Pan

class GCF(Pan):
  def __init__(self, bgc_genbanks_listing, gcf_id='GCF_X', lineage_name='Unnamed lineage'):
    super().__init__(bgc_genbanks_listing, lineage_name=lineage_name)
    self.gcf_id = gcf_id

    #######
    ## Variables not set during initialization
    #######

    # General variables
    self.hg_to_color = None

  def modifyPhylogenyForSamplesWithMultipleBGCs(self, input_phylogeny, result_phylogeny):
    try:
      number_of_added_leaves = 0
      t = Tree(input_phylogeny)
      for node in t.traverse('postorder'):
        if node.name in self.sample_bgcs and len(self.sample_bgcs[node.name]) > 1:
          og_node_name = node.name
          node.name = node.name + '_INNERNODE'
          for bgc_id in self.sample_bgcs[og_node_name]:
            # if bgc_id == node.name: continue
            node.add_child(name=bgc_id)
            child_node = t.search_nodes(name=bgc_id)[0]
            child_node.dist = 0
            if bgc_id != og_node_name: number_of_added_leaves += 1
      t.write(format=1, outfile=result_phylogeny)
      if self.logObject:
          self.logObject.info("New phylogeny with an additional %d leafs to reflect samples with multiple BGCs can be found at: %s." % (number_of_added_leaves, result_phylogeny))
    except Exception as e:
      if self.logObject:
        self.logObject.error("Had difficulties properly editing phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")
        self.logObject.error(traceback.format_exc())
      raise RuntimeError(traceback.format_exc())

  def assignColorsToHGs(self, gene_to_hg, bgc_genes):
    """
	:param gene_to_hg: gene to HG relationship.
	:param bgc_genes:  set of genes per HG.
	:return: dictionary mapping each HG to a hex color value.
	"""

    hg_bgc_counts = defaultdict(int)
    for b in bgc_genes:
      for g in bgc_genes[b]:
        if g in gene_to_hg:
          hg_bgc_counts[gene_to_hg[g]] += 1

    hgs = set([])
    for c in hg_bgc_counts:
      if hg_bgc_counts[c] > 1:
        hgs.add(c)

    # read in list of colors
    dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
    colors_file = dir_path + 'colors_200.txt'
    colors = []
    with open(colors_file) as ocf:
      colors = [x.strip() for x in ocf.readlines()]
    random.shuffle(colors)

    hg_to_color = {}
    for i, c in enumerate(set(hgs)):
      hg_to_color[c] = colors[i]
    self.hg_to_color = hg_to_color

  def createItolBGCSeeTrack(track_file, bgc_genes, gene_to_cog, cog_to_color, comp_gene_info, dataset_label, logObject):
    track_handle = open(track_file, 'w')

    logObject.info("Writing iTol track file to: %s" % track_file)
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
    bgc_gene_counts = defaultdict(int)
    for bgc in bgc_genes:
      bgc_gene_counts[bgc] = len(bgc_genes[bgc])

    for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
      bgc = item[0]
      curr_bgc_genes = bgc_genes[bgc]
      last_gene_end = max([comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
      printlist = [bgc, str(last_gene_end)]
      cog_directions = {}
      cog_lengths = defaultdict(list)
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
          cog_lengths[cog].append(gend - gstart)
      if i == 0:
        ref_cog_directions = cog_directions
        track_handle.write('\t'.join(printlist) + '\n')
      else:
        flip_support = 0
        keep_support = 0
        for c in ref_cog_directions:
          if not c in cog_directions: continue
          cog_weight = statistics.mean(cog_lengths[c])
          if cog_directions[c] == ref_cog_directions[c]:
            keep_support += cog_weight
          else:
            flip_support += cog_weight

        # flip the genbank visual if necessary, first BGC processed is used as reference guide
        if flip_support > keep_support:
          flip_printlist = printlist[:2]
          for gene_string in printlist[2:]:
            gene_info = gene_string.split('|')
            new_shape = None
            if gene_info[0] == 'TR':
              new_shape = 'TL'
            elif gene_info[0] == 'TL':
              new_shape = 'TR'
            new_gstart = int(last_gene_end) - int(gene_info[2])
            new_gend = int(last_gene_end) - int(gene_info[1])
            new_gene_info = '|'.join([new_shape, str(new_gstart), str(new_gend)] + gene_info[-2:])
            flip_printlist.append(new_gene_info)
          track_handle.write('\t'.join(flip_printlist) + '\n')
        else:
          track_handle.write('\t'.join(printlist) + '\n')
    track_handle.close()
