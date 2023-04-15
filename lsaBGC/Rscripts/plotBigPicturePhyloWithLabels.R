library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(phytools)
library(aplot)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

phylo.tree_file <- args[1]
track_file <- args[2]
pdf_file <- args[3]

phylo.tree <- read.tree(phylo.tree_file)
phylo.tree <- midpoint.root(phylo.tree)

if (track_file != "None") {
  track.data <- read.table(track_file, header=F, sep='\t')
  colnames(track.data) <- c('name', 'population')
}

tree.labels <- phylo.tree$tip.label

gg_tr <- ggtree(phylo.tree) + geom_tiplab()

if (track_file != "None" ){
  gg_tr <- gg_tr %<+% track.data + geom_tippoint(aes(color=as.factor(population)), show.legend=T, size=3)
}

pdf(pdf_file, height=20, width=40)
gg_tr
dev.off()