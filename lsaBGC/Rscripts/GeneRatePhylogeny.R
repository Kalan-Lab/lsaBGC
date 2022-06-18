library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(phytools)
library(aplot)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

phylo_file <- args[1]
track_file <- args[2]
pdf_file <- args[3]

phylo.tree <- read.tree(phylo_file)
phylo.tree <- midpoint.root(phylo.tree)
track_data <- read.table(track_file, header=F, sep='\t')

if (track_data[1,1] == "name") {
  col.names(track_data) <- track_data[,1]
} else {
  col.names(track_data) <- c('name', 'type')
}


pdf(pdf_file, height=10, width=10)
p <- ggtree(phylo.tree, layout="circular")
p <- p %<+% track_data + geom_tippoint(aes(color=as.factor(type)))
print(p)
dev.off()
