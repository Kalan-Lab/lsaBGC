library(ggplot2)
library(gggenes)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
tajimd.file <- args[1]
betard.file <- args[2]
propmc.file <- args[3]
conser.file <- args[4]
pdf.file <- args[5]

tajimd.data <- read.table(file=tajimd.file, sep='\t', header=T)
betard.data <- read.table(file=betard.file, sep='\t', header=T)
propmc.data <- read.table(file=propmc.file, sep='\t', header=T)
conser.data <- read.table(file=conser.file, sep='\t', header=T)

unique_gcfs <- unique(tajimd.data$GCF)

pdf(pdf.file, height=10, width=20)
for (gcf in unique_gcfs) {
  tajimd.data.gcf <- tajimd.data[tajimd.data$GCF == gcf, ]
  betard.data.gcf <- betard.data[betard.data$GCF == gcf, ]
  propmc.data.gcf <- propmc.data[propmc.data$GCF == gcf, ]
  conser.data.gcf <- conser.data[conser.data$GCF == gcf, ]

  print(gcf)
  g1 <- ggplot(tajimd.data.gcf, aes(xmin=start, xmax = end, y = GCF, forward = con_dir)) +
    geom_gene_arrow(aes(linetype=core, fill=Tajimas_D, alpha=proportion)) + theme_classic() +
    scale_fill_gradient(low='#e84a5d', high='#3763ad', na.value='grey50', guide='colourbar', aesthetics='fill') +
    ggtitle("Tajima's D") + theme(legend.position="bottom")

  g2 <- ggplot(betard.data.gcf, aes(xmin=start, xmax = end, y = GCF, forward = con_dir)) +
    geom_gene_arrow(aes(linetype=core, fill=Beta_RD, alpha=proportion)) + theme_classic() +
    scale_fill_gradient(low='#e84a5d', high='#3763ad', na.value='grey50', guide='colourbar', aesthetics='fill') +
    ggtitle("Beta-RD") + theme(legend.position="bottom")

  g3 <- ggplot(propmc.data.gcf, aes(xmin=start, xmax = end, y = GCF, forward = con_dir)) +
    geom_gene_arrow(aes(linetype=core, fill=Prop_Multi_Copy, alpha=proportion)) + theme_classic() +
    scale_fill_gradient(low='#e84a5d', high='#3763ad', na.value='grey50', guide='colourbar', aesthetics='fill') +
    ggtitle("Proportion of Samples with HG in Multi-Copy Across Genome") + theme(legend.position="bottom")

  g4 <- ggplot(conser.data.gcf, aes(x=reorder(HG, order), y=count, fill=population, color=core)) + theme_classic() +
    scale_color_manual(values=c('#FFFFFF', '#000000')) + geom_bar(stat='identity')+ scale_fill_brewer(palette='Set3') +
    ggtitle(paste0("Number of Samples with HG in this Specific GCF Context - ", gcf)) + theme(legend.position='bottom',
                                                                              axis.text.x = element_text(angle = 45,
                                                                                                         vjust = 1,
                                                                                                         hjust = 1))

  print(plot_grid(g4, g1, g2, g3, nrow=4, align='v'))
}

dev.off()