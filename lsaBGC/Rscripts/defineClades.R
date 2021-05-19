# Title     : createNJTreeAndDefineClades.R
# Objective : TODO
# Created by: rauf
# Created on: 4/27/21

library(ape)
library(phytools)
library(phylogram)
library(dynamicTreeCut)

args = commandArgs(trailingOnly=TRUE)

tree <- read.tree(args[1])
mid.njt <- midpoint.root(tree)
print(mid.njt)
hc <- as.hclust(mid.njt)
clusters <- cutreeDynamic(hc, minClusterSize=5, method="tree")

hc.labels <- labels(hc)
results <- cbind(hc.labels, clusters)f
write.table(results, file=args[2], col.names=F, row.names=F, quote=F, sep='\t')