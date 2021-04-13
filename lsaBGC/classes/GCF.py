import os
import sys
import logging
import traceback
import statistics
import random
import subprocess
import multiprocessing
from ete3 import Tree
from Bio import SeqIO
from operator import itemgetter
from collections import defaultdict
from lsaBGC.classes.Pan import Pan

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2])
RSCRIPT_FOR_BGSEE = lsaBGC_main_directory + '/lsaBGC/Rscripts/bgSee.R'

class GCF(Pan):
  def __init__(self, bgc_genbanks_listing, gcf_id='GCF_X', lineage_name='Unnamed lineage'):
    super().__init__(bgc_genbanks_listing, lineage_name=lineage_name)
    self.gcf_id = gcf_id

    #######
    ## Variables not set during initialization
    #######

    # General variables
    self.hg_to_color = None

    # Sequence and alignment directories
    self.nucl_seq_dir = None
    self.prot_seq_dir = None
    self.prot_alg_dir = None
    self.codo_alg_dir = None

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
        self.logObject.error("Had difficulties creating iTol track for visualization of BGC gene architecture.")
        self.logObject.error(traceback.format_exc())
      raise RuntimeError(traceback.format_exc())

  def visualizeGCFViaR(self, gggenes_track_file, heatmap_track_file, phylogeny_file, result_pdf_file):
    """
    Function to create tracks for visualization of gene architecture of BGCs belonging to GCF and run Rscript bgSee.R
    to produce automatic PDFs of plots. In addition, bgSee.R also produces a heatmap to more easily identify homolog
    groups which are conserved across isolates found to feature GCF.

    :param gggenes_track_file: Path to file with gggenes track information (will be created/written to by function, if it doesn't exist!)
    :param heatmap_track_file: Path to file for heatmap visual component (will be created/written to by function, if it doesn't exist!)
    :param phylogeny_file: Phylogeny to use for visualization.
    :param result_pdf_file: Path to PDF file where plots from bgSee.R will be written to.
    """
    try:
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
          curr_bgc_genes = self.bgc_genes[bgc]
          last_gene_end = max([self.comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
          printlist = []
          hg_directions = {}
          hg_lengths = defaultdict(list)
          for lt in curr_bgc_genes:
            ginfo = self.comp_gene_info[lt]
            hg = 'singleton'
            if lt in self.gene_to_hg:
              hg = self.gene_to_hg[lt]

            gstart = ginfo['start']
            gend = ginfo['end']
            forward = "FALSE"
            if ginfo['direction'] == '+': forward = "TRUE"

            hg_color = '"#dbdbdb"'
            if hg in self.hg_to_color:
              hg_color = '"' + self.hg_to_color[hg] + '"'

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
    except Exception as e:
      if self.logObject:
        self.logObject.error("Had difficulties creating tracks for visualization of BGC gene architecture along phylogeny using R libraries.")
        self.logObject.error(traceback.format_exc())
      raise RuntimeError(traceback.format_exc())

    rscript_plot_cmd = ["Rscript", RSCRIPT_FOR_BGSEE, phylogeny_file, gggenes_track_file, heatmap_track_file,
                        result_pdf_file]
    if self.logObject:
      self.logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
    try:
      subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=sys.stderr, stderr=sys.stderr,
                      executable='/bin/bash')
      self.logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
    except Exception as e:
      if self.logObject:
        self.logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
        self.logObject.error(traceback.format_exc())
      raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))

    if self.logObject:
      self.logObject.info('Plotting completed (I think successfully)!')

  def constructCodonAlignments(self, outdir, cores=1, only_scc=False):
    """
    Function to automate construction of codon alignments. This function first extracts protein and nucleotide sequnces
    from BGC Genbanks, then creates protein alignments for each homolog group using MAFFT, and finally converts those
    into codon alignments using PAL2NAL.

    :param outdir: Path to output/workspace directory. Intermediate files (like extracted nucleotide and protein
                   sequences, protein and codon alignments, will be writen to respective subdirectories underneath this
                   one).
    :param cores: Number of cores/threads to use when fake-parallelizing jobs using multiprocessing.
    :param only_scc: Whether to construct codon alignments only for homolog groups which are found to be core and in
                     single copy for samples with the GCF. Note, if working with draft genomes and the BGC is fragmented
                     this should be able to still identify SCC homolog groups across the BGC instances belonging to the
                     GCF.
    """

    nucl_seq_dir = os.path.abspath(outdir + 'Nucleotide_Sequences') + '/'
    prot_seq_dir = os.path.abspath(outdir + 'Protein_Sequences') + '/'
    prot_alg_dir = os.path.abspath(outdir + 'Protein_Alignments') + '/'
    codo_alg_dir = os.path.abspath(outdir + 'Codon_Alignments') + '/'

    if not os.path.isdir(nucl_seq_dir): os.system('mkdir %s' % nucl_seq_dir)
    if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
    if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
    if not os.path.isdir(codo_alg_dir): os.system('mkdir %s' % codo_alg_dir)

    all_samples = set(self.bgc_sample.values())
    try:
      inputs = []
      for hg in self.hg_genes:
        # if len(self.hg_genes[hg]) < 2: continue
        sample_counts = defaultdict(int)
        gene_sequences = {}
        for gene in self.hg_genes[hg]:
          gene_info = self.comp_gene_info[gene]
          bgc_id = gene_info['bgc_name']
          sample_id = self.bgc_sample[bgc_id]
          nucl_seq = gene_info['nucl_seq']
          prot_seq = gene_info['prot_seq']
          sample_counts[sample_id] += 1
          gid = sample_id + '|' + gene
          if only_scc:
            gid = sample_id
          gene_sequences[gid] = tuple([nucl_seq, prot_seq])
        samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
        # check that cog is single-copy-core
        if only_scc and len(samples_with_single_copy.symmetric_difference(all_samples)) > 0:
          continue
        elif only_scc and self.logObject:
          self.logObject.info('Homolog group %s detected as SCC across samples (not individual BGCs).' % hg)
        inputs.append([hg, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir])

      p = multiprocessing.Pool(cores)
      p.map(self.create_codon_msas, inputs)

      self.nucl_seq_dir = nucl_seq_dir
      self.prot_seq_dir = prot_seq_dir
      self.prot_alg_dir = prot_alg_dir
      self.codo_alg_dir = codo_alg_dir

    except Exception as e:
      if self.logObject:
        self.logObject.error("Issues with create protein/codon alignments of SCC homologs for BGC.")
        self.logObject.error(traceback.format_exc())
      raise RuntimeError(traceback.format_exc())

  def create_codon_msas(self, inputs):

    """
    Helper function which is to be called from the constructCodonAlignments() function to parallelize construction
    of codon alignments for each homolog group of interest in the GCF.

    :param inputs: list of inputs passed in by constructCodonAlignments().
    """
    hg, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir = inputs

    hg_nucl_fasta = nucl_seq_dir + '/' + hg + '.fna'
    hg_prot_fasta = prot_seq_dir + '/' + hg + '.faa'
    hg_prot_msa = prot_alg_dir + '/' + hg + '.msa.faa'
    hg_codo_msa = codo_alg_dir + '/' + hg + '.msa.fna'

    hg_nucl_handle = open(hg_nucl_fasta, 'w')
    hg_prot_handle = open(hg_prot_fasta, 'w')
    for s in gene_sequences:
      hg_nucl_handle.write('>' + s + '\n' + str(gene_sequences[s][0]) + '\n')
      hg_prot_handle.write('>' + s + '\n' + str(gene_sequences[s][1]) + '\n')
    hg_nucl_handle.close()
    hg_prot_handle.close()

    mafft_cmd = ['mafft', '--maxiterate', '1000', '--localpair', hg_prot_fasta, '>', hg_prot_msa]
    pal2nal_cmd = ['pal2nal.pl', hg_prot_msa, hg_nucl_fasta, '-output', 'fasta', '>', hg_codo_msa]

    if self.logObject:
      self.logObject.info('Running mafft with the following command: %s' % ' '.join(mafft_cmd))
    try:
      subprocess.call(' '.join(mafft_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                      executable='/bin/bash')
      if self.logObject:
        self.logObject.info('Successfully ran: %s' % ' '.join(mafft_cmd))
    except Exception as e:
      if self.logObject:
        self.logObject.error('Had an issue running: %s' % ' '.join(mafft_cmd))
        self.logObject.error(traceback.format_exc())
      raise RuntimeError('Had an issue running: %s' % ' '.join(mafft_cmd))

    if self.logObject:
      self.logObject.info('Running PAL2NAL with the following command: %s' % ' '.join(pal2nal_cmd))
    try:
      subprocess.call(' '.join(pal2nal_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                      executable='/bin/bash')
      if self.logObject:
        self.logObject.info('Successfully ran: %s' % ' '.join(pal2nal_cmd))
    except Exception as e:
      if self.logObject:
        self.logObject.error('Had an issue running: %s' % ' '.join(pal2nal_cmd))
        self.logObject.error(traceback.format_exc())
      raise RuntimeError('Had an issue running: %s' % ' '.join(pal2nal_cmd))

    if self.logObject:
      self.logObject.info('Achieved codon alignment for homolog group %s' % hg)

  def constructBGCPhylogeny(self, output_alignment, output_phylogeny):
    """
    Function to create phylogeny based on codon alignments of SCC homolog groups for GCF.

    :param output_alignment: Path to output file for concatenated SCC homolog group alignment.
    :param output_phylogeny: Path to output file for approximate maximum-likelihood phylogeny produced by FastTree2 from
                             concatenated SCC homolog group alignment.
    """
    try:
      bgc_sccs = defaultdict(lambda: "")
      fasta_data = []
      fasta_data_tr = []

      for f in os.listdir(self.codo_alg_dir):
        cog_align_msa = self.codo_alg_dir + f
        # concatenate gene alignments
        with open(cog_align_msa) as opm:
          for rec in SeqIO.parse(opm, 'fasta'):
            bgc_sccs['>' + rec.id] += str(rec.seq).upper()

      for b in bgc_sccs:
        fasta_data.append([b] + list(bgc_sccs[b]))

      for i, ls in enumerate(zip(*fasta_data)):
        if i == 0:
          fasta_data_tr.append(ls)
        else:
          n_count = len([x for x in ls if x == '-'])
          if (float(n_count) / len(ls)) < 0.1:
            fasta_data_tr.append(list(ls))

      scc_handle = open(output_alignment, 'w')

      for rec in zip(*fasta_data_tr):
        scc_handle.write(rec[0] + '\n' + ''.join(rec[1:]) + '\n')
      scc_handle.close()
    except Exception as e:
      if self.logObject:
        self.logObject.error('Had issues with creating concatenated alignment of the SCC homolog groups.')
        self.logObject.error(traceback.format_exc())
      raise RuntimeError(traceback.format_exc())

      # use FastTree2 to construct phylogeny
      fasttree_cmd = ['fasttree', '-nt', output_alignment, '>', output_phylogeny]
      logObject.info('Running FastTree2 with the following command: %s' % ' '.join(fasttree_cmd))
      try:
        subprocess.call(' '.join(fasttree_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                        executable='/bin/bash')
        logObject.info('Successfully ran: %s' % ' '.join(fasttree_cmd))
      except Exception as e:
        if self.logObject:
          self.logObject.error('Had an issue running: %s' % ' '.join(fasttree_cmd))
          self.logObject.error(traceback.format_exc())
        raise RuntimeError('Had an issue running: %s' % ' '.join(fasttree_cmd))
