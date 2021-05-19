# Title     : createNJTreeAndDefineClades.R
# Objective : TODO
# Created by: rauf
# Created on: 4/27/21

library(ape)
library(phytools)

args = commandArgs(trailingOnly=TRUE)

dat <- read.table(args[1], header=T, sep='\t', row.names=1)
d <- as.dist(dat)
njt <- nj(d)
mid.njt <- midpoint.root(njt)
write.tree(mid.njt, file=args[2])