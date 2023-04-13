library(ggplot2)
library(gggenes)

args = commandArgs(trailingOnly=TRUE)
input.file <- args[1]
height <- as.numeric(args[2])
width <- as.numeric(args[3])
pdf.file <- args[4]

input.data <- read.table(file=input.file, sep='\t', header=T)
# tf_handle.write('\t'.join(['gene_cluster', 'CDS_feature', 'BGC_likeness', 'CDS_start', 'CDS_end', 'CDS_dir']) + '\n')

pdf(pdf.file, height=height, width=width)
ggplot(input.data, aes(xmin=CDS_start, xmax = CDS_end, y = "", forward = CDS_dir)) +
  geom_gene_arrow(aes(fill=BGC_likeness)) + theme_genes() +
  scale_fill_gradient("BGC likeness", low='#a7c4f2', high='#1f498c', breaks =c(0.00, 0.25, 0.50, 0.75, 1.00),
                      labels=c("0.00", "0.25", "0.50", "0.75", "1.00"), limits=c(0.0,1.0), na.value='grey50',
                      guide='colourbar', aesthetics='fill') + facet_wrap(~gene_cluster, ncol=1) +
  theme(legend.position="bottom") + theme(strip.text = element_text(size = 10, color = "black")) + xlab("") + ylab("")
dev.off()