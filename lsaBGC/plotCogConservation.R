library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

domains.dat <- read.table(args[1], header=T, sep='\t')
# domain  domain_index    min_pos max_pos
positions.dat <- read.table(args[2], header=T, sep='\t')
# pos     num_seqs        num_alleles     num_gaps        maj_allele_freq
popgen.dat <- read.table(args[3], header=T, sep='\t')
# pos     type

max_pos <- max(positions.dat$pos)
colors <- c("NS" = "#db2c1d", "S" = "#999897")
names(colors) <- c('NS', 'S')

maj_allele_gg <- ggplot(positions.dat, aes(x=pos, y=(1-maj_allele_freq))) + geom_line() + theme_bw() + xlab("") +
                 ylab("1.0 - Major Allele Frequency") + xlim(0, max_pos+1)
allele_count_gg <- ggplot(positions.dat, aes(x=pos, y=num_alleles)) + geom_bar(stat='identity', fill='black') +
                 theme_classic() + xlab("") + ylab("Count of Alleles") + xlim(0, max_pos+1)+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
site_coverage_gg <- ggplot(positions.dat, aes(x=pos, y=(num_seqs-num_gaps))) + geom_bar(stat='identity', fill='black') +
                 theme_classic() + xlab("") + ylab("Site Coverage") + xlim(0, max_pos+1)+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

dnds_gg <- ggplot(popgen.dat, aes(x=pos, y=0, color=type)) + geom_point(show.legend=F, alpha=0.7) + theme_void() +
           xlab("") + ylab("Variable Positions") + xlim(0, max_pos+1) + scale_color_manual(values=colors)

doms_gg <- ggplot(domains.dat, aes(color=domain, x= min_pos, xend=max_pos, y=reorder(domain, domain_index), yend=reorder(domain, domain_index))) +
           geom_segment(show.legend=F, size=5) + xlab("Position Along MSA") + ylab("Domains") + scale_color_brewer(palette='Set1') + xlim(0, max_pos+1) +
           geom_label(aes(label=domain),vjust = 0, nudge_y = 0.1, show.legend=F) + theme_grey() +
           theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

pdf(args[4], height=10, width=21)

if (nrow(popgen.dat) > 0 && nrow(doms_gg) > 0) {
  plot_grid(site_coverage_gg,
            allele_count_gg,
            dnds_gg,
            maj_allele_gg,
            doms_gg,
            rel_heights=c(0.85,0.85,0.25,4,2.5),
            align='v',
            axis='b',
            ncol=1)
} else if (nrow(popgen.dat) > 0) {
    plot_grid(site_coverage_gg,
            allele_count_gg,
            dnds_gg,
            maj_allele_gg,
            rel_heights=c(0.85,0.85,0.25,4),
            align='v',
            axis='b',
            ncol=1)
} else {
    plot_grid(site_coverage_gg,
            allele_count_gg,
            maj_allele_gg,
            doms_gg,
            rel_heights=c(0.85,0.85,4,2.5),
            align='v',
            axis='b',
            ncol=1)
}
dev.off()