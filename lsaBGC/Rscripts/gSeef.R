library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(phytools)
library(aplot)
library(dplyr)
library(grid)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

phylo.tree_file <- args[1]
heatmap.data_file <- args[2]
png_file <- args[3]
legend_png_file <- args[4]

phylo.tree <- read.tree(phylo.tree_file)
heatmap.data <- read.table(heatmap.data_file, header=T, sep='\t')

tree.labels <- phylo.tree$tip.label

heatmap.data.select <- distinct(heatmap.data[c("annotation", "colors")])
gcf_colors <- c(heatmap.data.select$colors)
names(gcf_colors) <- c(heatmap.data.select$annotation)
print(gcf_colors)
png(png_file, height=10, width=20, units='in', res=300)
gg_tr <- ggtree(phylo.tree)
gg_hm <- ggplot(heatmap.data, aes(x = reorder(gcf, gcf_index), y = label, fill = annotation)) +
         theme_classic() + scale_fill_manual(values=gcf_colors) +
         xlab("GCF IDs") + ylab("") + geom_tile(color='white', show.legend=F) +
		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg_hm %>% insert_left(gg_tr, width=0.4)
dev.off()

png(legend_png_file, height=10, width=5, units='in', res=300)
my_hist <- ggplot(heatmap.data, aes(x=annotation, y=1, fill = annotation)) + geom_bar(stat='identity') + scale_fill_manual(values=gcf_colors)
legend <- cowplot::get_legend(my_hist)
grid.newpage()
grid.draw(legend)
dev.off()