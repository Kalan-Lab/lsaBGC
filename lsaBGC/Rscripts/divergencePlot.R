library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
divergence.report <- args[1]
result.pdf <- args[2]

dat <- read.table(divergence.report, header=T, sep='\t')
dat.filt.90 <- dat[dat$gcf_content_sim >= 0.90,]
dat.filt.75 <- dat[dat$gcf_content_sim >= 0.75,]
dat.filt.50 <- dat[dat$gcf_content_sim >= 0.50,]

pdf(result.pdf, height=10, width=20)

g1 <- ggplot(dat.filt.50, aes(x=gcf_id, y=beta_rd)) + geom_boxplot(fill='grey') + theme_classic() + geom_hline(yintercept=1.0, linetype=2) + ggtitle("Jaccard Index >= 90")
g2 <- ggplot(dat.filt.50, aes(x=gcf_id, y=beta_rd)) + geom_boxplot(fill='grey') + theme_classic() + geom_hline(yintercept=1.0, linetype=2) + ggtitle("Jaccard Index >= 75")
g3 <- ggplot(dat.filt.50, aes(x=gcf_id, y=beta_rd)) + geom_boxplot(fill='grey') + theme_classic() + geom_hline(yintercept=1.0, linetype=2) + ggtitle("Jaccard Index >= 50")

plot_grid(g1, g2, g3, rel_heights=c(0.3, 0.3, 0.3), align='v', axis='b', ncol=1)

dev.off()

