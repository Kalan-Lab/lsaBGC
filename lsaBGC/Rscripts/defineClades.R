# Title     : createNJTreeAndDefineClades.R
# Objective : TODO
# Created by: rauf
# Created on: 4/27/21

library(ape)
library(phytools)
library(phylogram)
library(dynamicTreeCut)
library(ggtree)

args = commandArgs(trailingOnly=TRUE)

tree <- read.tree(args[1])

dist <- as.matrix(as.dist(read.table(args[2], header=T, sep='\t', row.names=1)))

mid.njt <- midpoint.root(tree)
mid.njt <- force.ultrametric(mid.njt, method=c("nnls"))
mid.njt <-multi2di(mid.njt)

hc <- as.hclust(mid.njt)
clusters <- cutreeDynamic(hc, distM=dist, minClusterSize=1, deepSplit=4, method="hybrid", useMedoids=T)

hc.labels <- mid.njt$tip.label
results <- cbind(hc.labels, clusters)
write.table(results, file=args[3], col.names=F, row.names=F, quote=F, sep='\t')

pdf(args[4], height=10, width=14)
results.data.frame <- as.data.frame(cbind(hc.labels, clusters))
colnames(results.data.frame) <- c("node", "clade")
rownames(results.data.frame) <- results.data.frame$node
results.data.frame <-  results.data.frame[,2, drop=FALSE]
print(head(results.data.frame))
circ <- ggtree(midpoint.root(tree))

p1 <- gheatmap(circ, results.data.frame)
print(p1)
dev.off()