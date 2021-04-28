# Title     : createNJTreeAndDefineClades.R
# Objective : TODO
# Created by: rauf
# Created on: 4/27/21

library(ape)
library(phytools)
library(treestructure)

args = commandArgs(trailingOnly=TRUE)

dat <- read.table(args[1], header=T, sep='\t', row.names=1)
d <- as.dist(dat)
njt <- nj(d)
mid.njt <- midpoint.root(njt)
write.tree(mid.njt, file=args[2])

s <- trestruct(mid.njt)
structureData <- as.data.frame(s)

write.table(structureData, file=args[3], col.names=F, row.names=F, quote=F, sep='\t')