# Title     : TODO
# Objective : TODO
# Created by: rauf
# Created on: 7/12/21
library(ggplot2)
library(gggenes)

args = commandArgs(trailingOnly=TRUE)
pop.gene.file <- args[1]
pdf.file <- args[2]

pop.gene.data <- read.table(file=pop.gene.file, sep='\t', header=T)

pdf(pdf.file, height=10, width=20)
ggplot(pop.gene.data, aes(xmin = gene_start, xmax = gene_stop, y = gcf_id, forward = hg_consensus_direction)) +
  geom_gene_arrow(aes(linetype=is_core_to_bgc, fill=Tajimas_D, alpha=proportion_of_samples_with_hg)) + theme_classic() +
  scale_fill_gradient(low = "#e84a5d", high="#3763ad", na.value = "grey50", guide = "colourbar",  aesthetics = "fill") + ggtitle("Tajima's D")


ggplot(pop.gene.data, aes(xmin = gene_start, xmax = gene_stop, y = gcf_id, forward = hg_consensus_direction)) +
  geom_gene_arrow(aes(linetype=is_core_to_bgc, fill=median_beta_rd, alpha=proportion_of_samples_with_hg)) + theme_classic() +
  scale_fill_gradient(low = "#e84a5d", high="#3763ad", na.value = "grey50", guide = "colourbar",  aesthetics = "fill") + ggtitle("Median Beta-RD")


dev.off()