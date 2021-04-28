library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(tidyverse)
library(phytools)
library(aplot)

args = commandArgs(trailingOnly=TRUE)

phylo.tree_file <- args[1]
genes.data_file <- args[2]
heatmap.data_file <- args[3]
pdf_file <- args[4]

phylo.tree <- read.tree(phylo.tree_file)
phylo.tree <- midpoint.root(phylo.tree)
genes.data <- as.tibble(read.table(genes.data_file, header=T, sep='\t'))
heatmap.data <- as.tibble(read.table(heatmap.data_file, header=T, sep='\t'))

og_colors <- c(genes.data$og_color, '#FFFFFF')
names(og_colors) <- c(genes.data$og, 'Absent')

pdf(pdf_file, height=30, width=40)
gg_tr <- ggtree(phylo.tree) + geom_tiplab(align=TRUE) + theme_void()  + ggplot2::xlim(0, 0.5)
gg_gn <- ggplot(genes.data, aes(xmin = start, xmax = end, y = label, fill = og, forward = forward)) +
          geom_gene_arrow(show.legend=F) + theme_void() + scale_fill_manual(values=og_colors)
gg_hm <- ggplot(heatmap.data, aes(x=reorder(og, -og_count), y=label, fill=og_presence)) +
         theme_classic() + scale_fill_manual(values=og_colors) +#values=c("#FFFFFF", "#000000")) +
         xlab("Homolog Group IDs") + ylab("BGC IDs") + geom_tile(color='white', show.legend=F) + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gg_hm %>% insert_left(gg_gn, width=1.0) %>% insert_left(gg_tr, width=0.4)

dev.off()