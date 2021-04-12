import os
import sys
import logging
import traceback
import statistics
import random
import subprocess
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
    """
    Function which takes in an input phylogeny and produces a replicate resulting phylogeny with samples/leafs which
    have multiple BGC instances for a GCF expanded.

    :param input_phylogeny: input newick phylogeny file
    :result result_phylogeny: resulting newick phylogeny file
    """
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

  def createItolBGCSeeTrack(self, result_track_file):
    """
    Function to create a track file for visualizing BGC gene architecture across a phylogeny in the interactive tree
    of life (iTol)

    :param result_track_file: The path to the resulting iTol track file for BGC gene visualization.
    """
    try:
      track_handle = open(result_track_file, 'w')

      self.logObject.info("Writing iTol track file to: %s" % result_track_file)
      self.logObject.info("Track will have label: %s" % self.gcf_id)

      # write header for iTol track file
      track_handle.write('DATASET_DOMAINS\n')
      track_handle.write('SEPARATOR TAB\n')
      track_handle.write('DATASET_LABEL\t%s\n' % self.gcf_id)
      track_handle.write('COLOR\t#000000\n')
      track_handle.write('BORDER_WIDTH\t1\n')
      track_handle.write('BORDER_COLOR\t#000000\n')
      track_handle.write('SHOW_DOMAIN_LABELS\t0\n')
      track_handle.write('DATA\n')

      # write the rest of the iTol track file for illustrating genes across BGC instances
      ref_hg_directions = {}
      bgc_gene_counts = defaultdict(int)
      for bgc in self.bgc_genes:
        bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])

      for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
        bgc = item[0]
        curr_bgc_genes = self.bgc_genes[bgc]
        last_gene_end = max([self.comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
        printlist = [bgc, str(last_gene_end)]
        hg_directions = {}
        hg_lengths = defaultdict(list)
        for lt in curr_bgc_genes:
          ginfo = self.comp_gene_info[lt]
          hg = 'singleton'
          if lt in self.gene_to_hg:
            hg = self.gene_to_hg[lt]
          shape = 'None'
          if ginfo['direction'] == '+':
            shape = 'TR'
          elif ginfo['direction'] == '-':
            shape = 'TL'
          gstart = ginfo['start']
          gend = ginfo['end']
          hg_color = "#dbdbdb"
          if hg in self.hg_to_color:
            hg_color = self.hg_to_color[hg]
          gene_string = '|'.join([str(x) for x in [shape, gstart, gend, hg_color, hg]])
          printlist.append(gene_string)
          if hg != 'singleton':
            hg_directions[hg] = ginfo['direction']
            hg_lengths[hg].append(gend - gstart)
        if i == 0:
          ref_hg_directions = hg_directions
          track_handle.write('\t'.join(printlist) + '\n')
        else:
          flip_support = 0
          keep_support = 0
          for c in ref_hg_directions:
            if not c in hg_directions: continue
            hg_weight = statistics.mean(hg_lengths[c])
            if hg_directions[c] == ref_hg_directions[c]:
              keep_support += hg_weight
            else:
              flip_support += hg_weight

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
    except Exception as e:
      if self.logObject:
        self.logObject.error(
          "Had difficulties properly editing phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")
        self.logObject.error(traceback.format_exc())
      raise RuntimeError(traceback.format_exc())

  def visualizeGCFViaR(self, gggenes_track_file, heatmap_track_file, phylogeny_file, result_pdf_file):

    if not os.path.isfile(gggenes_track_file) or not os.path.isfile(heatmap_track_file):
      gggenes_track_handle = open(gggenes_track_file, 'w')
      heatmap_track_handle = open(heatmap_track_file, 'w')
      if self.logObject:
        self.logObject.info("Writing gggenes input file to: %s" % gggenes_track_file)
        self.logObject.info("Writing heatmap input file to: %s" % heatmap_track_file)
      # write header for track files
      gggenes_track_handle.write('label\tgene\tstart\tend\tforward\tog\tog_color\n')
      heatmap_track_handle.write('label\tog\tog_presence\tog_count\n')

      ref_hg_directions = {}

      bgc_gene_counts = defaultdict(int)
      for bgc in self.bgc_genes:
        bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])

      tree_obj = Tree(phylogeny_file)
      bgc_weights = defaultdict(int)
      for leaf in tree_obj:
        bgc_weights[str(leaf).strip('\n').lstrip('-')] += 1

      bgc_hg_presence = defaultdict(lambda: defaultdict(lambda: 'Absent'))
      hg_counts = defaultdict(int)
      for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
        bgc = item[0]
        curr_bgc_genes = bgc_genes[bgc]
        last_gene_end = max([comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
        printlist = []
        hg_directions = {}
        hg_lengths = defaultdict(list)
        for lt in curr_bgc_genes:
          ginfo = comp_gene_info[lt]
          hg = 'singleton'
          if lt in gene_to_hg:
            hg = gene_to_hg[lt]

          gstart = ginfo['start']
          gend = ginfo['end']
          forward = "FALSE"
          if ginfo['direction'] == '+': forward = "TRUE"

          hg_color = '"#dbdbdb"'
          if hg in hg_to_color:
            hg_color = '"' + hg_to_color[hg] + '"'

          gene_string = '\t'.join([str(x) for x in [bgc, lt, gstart, gend, forward, hg, hg_color]])
          printlist.append(gene_string)
          if hg != 'singleton':
            bgc_hg_presence[bgc][hg] = hg
            hg_counts[hg] += bgc_weights[bgc]
            hg_directions[hg] = ginfo['direction']
            hg_lengths[hg].append(gend - gstart)
        if i == 0:
          ref_hg_directions = hg_directions
          gggenes_track_handle.write('\n'.join(printlist) + '\n')
        else:
          flip_support = 0
          keep_support = 0
          for c in ref_hg_directions:
            if not c in hg_directions: continue
            hg_weight = statistics.mean(hg_lengths[c])
            if hg_directions[c] == ref_hg_directions[c]:
              keep_support += hg_weight
            else:
              flip_support += hg_weight

          # flip the genbank visual if necessary, first BGC processed is used as reference guide
          if flip_support > keep_support:
            flip_printlist = []
            for gene_string in printlist:
              gene_info = gene_string.split('\t')
              new_forward = 'TRUE'
              if gene_info[4] == 'TRUE': new_forward = 'FALSE'
              new_gstart = int(last_gene_end) - int(gene_info[3])
              new_gend = int(last_gene_end) - int(gene_info[2])
              new_gene_string = '\t'.join([str(x) for x in
                                           [gene_info[0], gene_info[1], new_gstart, new_gend, new_forward,
                                            gene_info[-2], gene_info[-1]]])
              flip_printlist.append(new_gene_string)
            gggenes_track_handle.write('\n'.join(flip_printlist) + '\n')
          else:
            gggenes_track_handle.write('\n'.join(printlist) + '\n')
      gggenes_track_handle.close()

      for bgc in bgc_hg_presence:
        for hg in hg_counts:
          heatmap_track_handle.write('\t'.join([bgc, hg, bgc_hg_presence[bgc][hg], str(hg_counts[hg])]) + '\n')
      heatmap_track_handle.close()

    rscript_plot_cmd = ["Rscript", RSCRIPT_FOR_BGSEE, phylogeny_file, gggenes_track_file, heatmap_track_file,
                        result_pdf_file]
    logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
    try:
      subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=sys.stderr, stderr=sys.stderr,
                      executable='/bin/bash')
      logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
    except:
      logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
      raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
    logObject.info('Plotting completed!')
