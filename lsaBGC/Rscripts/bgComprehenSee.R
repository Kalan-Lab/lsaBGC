library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(phytools)
library(aplot)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

phylo.tree_file <- args[1]
heatmap.data_file <- args[2]
track_file <- args[3]
pdf_file <- args[4]

phylo.tree <- read.tree(phylo.tree_file)
phylo.tree <- midpoint.root(phylo.tree)
heatmap.data <- read.table(heatmap.data_file, header=T, sep='\t')
track.data <- read.table(track_file, header=T, sep='\t')

tree.labels <- phylo.tree$tip.label

og_colors <- c('#FFFFFF', '#949292', '#403f3f')
names(og_colors) <- c('0', '1', 'Multi')

pdf(pdf_file, height=30, width=30)
gg_tr <- ggtree(phylo.tree)#  + ggplot2::xlim(NA, 1)
gg_tr <- gg_tr %<+% track.data + geom_tippoint(aes(color=detection_method), show.legend=T, size=3)
gg_hm <- ggplot(heatmap.data, aes(x = og, y = label, fill=as.factor(og_copy))) +
         theme_classic() + scale_fill_manual('Copy Count', values=og_colors) +
         xlab("Homolog Group IDs") + ylab("BGC IDs") + geom_tile(color='white') +
		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg_hm %>% insert_left(gg_tr, width=0.3)
dev.off()
