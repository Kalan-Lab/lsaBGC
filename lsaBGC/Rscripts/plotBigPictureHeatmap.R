library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(phytools)
library(aplot)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

phylo.tree_file <- args[1]
heatmap.data_file <- args[2]
track_file <- args[3]
pdf_file <- args[4]

phylo.tree <- read.tree(phylo.tree_file)
phylo.tree <- midpoint.root(phylo.tree)
heatmap.data <- read.table(heatmap.data_file, header=T, sep='\t')

if (track_file != "None") {
  track.data <- read.table(track_file, header=F, sep='\t')
  colnames(track.data) <- c('name', 'population')
}

tree.labels <- phylo.tree$tip.label

pdf(pdf_file, height=20, width=40)
gg_tr <- ggtree(phylo.tree)

if (track_file != "None" ){
  gg_tr <- gg_tr %<+% track.data + geom_tippoint(aes(color=as.factor(population)), show.legend=F, size=3)
}

gg_hm <- ggplot(heatmap.data, aes(x = as.numeric(Homolog_Group_Order), y = label, fill=log(Difference_to_Consensus_Sequence+1e-5,10))) +
  geom_tile() + scale_fill_gradient("similarity to consensus\nsequence (blue=\nmore similar)", low = "darkblue", high = "white", na.value="grey") + theme_minimal() +
  facet_grid(. ~ GCF_Order, scales = "free_x", space='free') +
  theme(text = element_text(size=20), axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color='red') 
gg_hm %>% insert_left(gg_tr, width=0.2)
dev.off()
