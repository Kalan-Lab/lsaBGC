library(scatterpie)
library(ggplot2)
library(RColorBrewer)
library(ggforce)

args = commandArgs(trailingOnly=TRUE)

# read in plotting inputs
data <- read.table(args[1], header=T, sep='\t')
#print(data)
annot_classes <- scan(args[2], what="", sep="\n")

#print(annot_classes)
nb.cols <- length(annot_classes)-1
mycolors <- colorRampPalette(brewer.pal(nb.cols, "Set2"))(nb.cols)
mycolors <- c("#808080", mycolors)
names(mycolors) <- annot_classes

npages <- length(unique(data$Parameters))

# GCF     Parameters      Samples MultiBGCSamples SCCExists       CoreGeneClusters
pdf(args[3], height=5, width=7)
for(i in 1:npages){
  print(ggplot() + geom_scatterpie(aes(x=Samples, y=MultiBGCSamples), alpha=0.4, data=data, cols=annot_classes, show.legend=F) +
    coord_equal() + scale_fill_manual(values=mycolors) +  facet_wrap_paginate(~Parameters, ncol = 1, nrow = 1, page=i) +
    theme_classic() + xlab("Number of Samples with GCF") + ylab("Number of Samples with >= 2 BGCs for GCF") +
    geom_point(aes(x=Samples, y=MultiBGCSamples, shape=as.factor(SCCExists)), size=1.3, show.legend=F, data=data, alpha=0.5))
}
dev.off()