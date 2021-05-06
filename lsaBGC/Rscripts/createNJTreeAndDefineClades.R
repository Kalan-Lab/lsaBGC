()# Title     : createNJTreeAndDefineClades.R
# Objective : TODO
# Created by: rauf
# Created on: 4/27/21

library(ape)
library(phytools)
library(treestructure)
library(ggtree)

args = commandArgs(trailingOnly=TRUE)

dat <- read.table(args[1], header=T, sep='\t', row.names=1)
d <- as.dist(dat)
njt <- nj(d)
mid.njt <- midpoint.root(njt)
write.tree(mid.njt, file=args[2])

s <- trestruct(mid.njt)

pdf(args[3], height=10, width=7)
plot(s)  + ggtree::geom_tiplab()
dev.off()

structureData <- as.data.frame(s)

write.table(structureData, file=args[4], col.names=F, row.names=F, quote=F, sep='\t')