library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(phytools)
library(aplot)
library(dplyr)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

phylo.tree_file <- args[1]
heatmap_file <- args[2]
# 'hg', 'hg_order', 'hg_start', 'hg_end', 'hg_direction', 'label', 'median_consensus_difference'
species.gcf.count_file <- args[3]
# 'label', 'isolates_with_gcf'
popgen_file <- args[4]
# 'hg', 'hg_order', 'hg_start', 'hg_end', 'samples_with_hg', 'tajimas_d', 'beta_rd'
annotation_file <- args[5]
# 'hg', 'hg_order', 'hg_start', 'hg_end', 'manual_annotation'
pdf_file_phylogeny <- args[6]
pdf_file_popgen <- args[7]

# ggplot2 objects for first figure: (1) tree/phylogeny, (2) heatmap of consensus difference,
#                                   (3) bargraph for count of samples with gcf per species
# ggplot2 objects for second figure: (1) manual annotation track (2) bar graph (log-10-scale) for number of samples,
#                                    (3) bar graph for beta-rd around 1.0, (4) bar graph for tajima's d around 0.0

phylo.tree <- read.tree(phylo.tree_file)
phylo.tree <- midpoint.root(phylo.tree)

heatmap.data <- read.table(heatmap_file, header=T, sep='\t')
species.gcf.count.data <- read.table(species.gcf.count_file, header=T, sep='\t')
popgen.data <- read.table(popgen_file, header=T, sep='\t')
annotation.data <- read.table(annotation_file, header=T, sep='\t')

tree.labels <- phylo.tree$tip.label

pdf(pdf_file_phylogeny, height=10, width=30)
gg_tr <- ggtree(phylo.tree)
gg_gn <- ggplot(heatmap.data, aes(xmin = hg_start, xmax = hg_end, y = label, fill = median_consensus_difference,
								forward = hg_direction)) +
	     geom_gene_arrow(show.legend=F) + theme_void() + scale_fill_gradient(low='black', high='white') + ylab("") +
	     xlab("")
gg_br <- ggplot(species.gcf.count.data, aes(y = label, x = log(isolates_with_gcf, 10))) + theme_classic() +
	     geom_col(stat='identity', fill='black') + xlab("Isolates with GCF") + ylab("") +
	     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg_gn %>% insert_left(gg_tr, width=0.4) %>% insert_right(gg_br, width=0.1)
dev.off()

colors <- c('#f03a3a', '#b80f0f', '#702f2f', '#636262', '#1a77ad', '#1c7534', '#d49c11')
names(colors) <- c('MGE', 'MGE - Phage', 'Addiction System', 'Other', 'Other - Regulatory', 'Transport',
				   'Overlaps Protocore Biosynthesis Region')

pdf(pdf_file_popgen, height=20, width=20)
gg_annot <- ggplot(annotation.data, aes(xmin = hg_start, xmax = hg_end, ymin = 0, ymax = 0, fill = manual_annotation)) +
	        geom_segment(color='black') + theme_classic()
gg_br_ns <- ggplot(popgen.data, aes(xmin = hg_start, xmax = hg_end, ymin = 0, ymax = samples_with_hg)) + theme_classic() +
			geom_segment(fill='black')
gg_br_rd <- ggplot(popgen.data, aes(xmin = hg_start, xmax = hg_end, ymin = 0, ymax = beta_rd)) + theme_classic() +
			geom_segment(fill='black')
gg_br_td <- ggplot(popgen.data, aes(xmin = hg_start, xmax = hg_end, ymin = 0, ymax = tajimas_d)) + theme_classic() +
			geom_segment(fill='black')
plot_grid(gg_annot, gg_br_ns, gg_br_rd, gg_br_td, ncol=1, align = 'v', axis = 'l', rel_heights=c(1,2,2,2))
dev.off()