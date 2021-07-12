# Title     : TODO
# Objective : TODO
# Created by: rauf
# Created on: 7/9/21

library(ape)
library(phytools)
library(phylogram)
library(dynamicTreeCut)
library(ggtree)

args = commandArgs(trailingOnly=TRUE)

tree_file <- args[1]
distance_matrix <- args[2]
deep_split_value <- args[3]
txt_output_file <- args[4]
pdf_output_file <- args[5]

tree <- read.tree(tree_file)

distM <- as.matrix(as.dist(read.table(distance_matrix, header=T, sep='\t', row.names=1)))

mid.njt <- midpoint.root(tree)
mid.njt <- force.ultrametric(mid.njt, method=c("nnls"))
mid.njt <-multi2di(mid.njt)

clusters <- cutreeDynamic(hc, distM=distM, minClusterSize=1, deepSplit=as.numeric(deep_split_value), method="hybrid")
hc.labels <- mid.njt$tip.label
results <- cbind(hc.labels, clusters)
write.table(results, file=txt_output_file, col.names=F, row.names=F, quote=F, sep='\t')

pdf(pdf_output_file, height=10, width=14)
results.data.frame <- as.data.frame(cbind(hc.labels, clusters))
colnames(results.data.frame) <- c("node", "clade")
rownames(results.data.frame) <- results.data.frame$node
results.data.frame <-  results.data.frame[,2, drop=FALSE]
print(head(results.data.frame))
circ <- ggtree(midpoint.root(tree))

p1 <- gheatmap(circ, results.data.frame)
print(p1)
dev.off()
