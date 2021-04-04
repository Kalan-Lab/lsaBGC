library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(tidyverse)
library(cowplot)
library(phytools)
args = commandArgs(trailingOnly=TRUE)

### following function taken/adapted from Thomas Hackl's blog: https://thackl.github.io/ggtree-composite-plots
# overwrite the default expand for continuous scales
scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
    scale_y_continuous(expand=expand, ...)
}

tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree")) {
    stop("not a ggtree object")
  }
  left_join(select(data, label), select(ggtree$data, label, y)) %>% pull(y)
}

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

pdf(pdf_file, height=10, width=20)

gg_tr <- ggtree(phylo.tree) + geom_tiplab(align=TRUE, size=0) + scale_x_continuous(expand=expand_scale(0.2)) + scale_y_tree()
gg_gn <- ggplot(genes.data, aes(xmin = start, xmax = end, y = tree_y(gg_tr, genes.data), fill = og, forward = forward)) +
          geom_gene_arrow(show.legend=F) + theme_void() + scale_fill_manual(values=og_colors)
plot_grid(gg_tr, gg_gn, rel_widths=c(1,4))

#label\tog\tog_presence\tog_count

gg_hm <- ggplot(heatmap.data, aes(x=reorder(og, -og_count), y=reorder(label, tree_y(gg_tr, heatmap.data)), fill=og_presence)) +
         theme_classic() + scale_fill_manual(values=og_colors) +#values=c("#FFFFFF", "#000000")) +
         xlab("Homolog Group IDs") + ylab("BGC IDs") + geom_tile(color='white', show.legend=F) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(gg_hm)
dev.off()