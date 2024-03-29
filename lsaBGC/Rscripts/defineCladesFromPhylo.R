# Title     : defineCladesFromPhylo.R
# Objective : TODO
# Created by: rauf
# Created on: 4/27/21

library(ape)
library(phytools)
library(ggtree)

args = commandArgs(trailingOnly=TRUE)

tree_file <- args[1]
k <- args[2]
txt_output_file <- args[3]
pdf_output_file <- args[4]

tree <- read.tree(tree_file)
mid.njt <- midpoint.root(tree)
mid.njt <- force.ultrametric(mid.njt, method=c("nnls"))
mid.njt <-multi2di(mid.njt)
hc <- as.hclust(mid.njt)
clusters <- cutree(hc, k=as.numeric(k))
hc.labels <- mid.njt$tip.label
results <- cbind(hc.labels, clusters)
write.table(results, file=txt_output_file, col.names=F, row.names=F, quote=F, sep='\t')

results.data.frame <- as.data.frame(cbind(hc.labels, clusters))
colnames(results.data.frame) <- c("name", "type")

pdf(pdf_output_file, height=10, width=14)
p <- ggtree(tree, layout="circular")
p <- p %<+% results.data.frame + geom_tippoint(aes(color=as.factor(type)))
print(p)
dev.off()
