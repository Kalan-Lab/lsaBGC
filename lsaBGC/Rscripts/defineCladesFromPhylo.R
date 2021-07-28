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

tree_file <- args[1]
deep_split_value <- args[2]
txt_output_file <- args[3]
pdf_output_file <- args[4]

tree <- read.tree(tree_file)

mid.njt <- midpoint.root(tree)
mid.njt <- force.ultrametric(mid.njt, method=c("nnls"))
mid.njt <-multi2di(mid.njt)

hc <- as.hclust(mid.njt)
distM <- cophenetic(tree)
clusters <- cutreeDynamic(hc, distM=distM, minClusterSize=1, deepSplit=as.numeric(deep_split_value), method="hybrid")#, cutHeight=0.6)

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
