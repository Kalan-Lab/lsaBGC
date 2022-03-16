library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

pairwise_distance_file <- args[1]
sample_info_file <- args[2]
pdf_file <- args[3]

pwdist_data <- read.table(pairwise_distance_file, header=T, sep='\t', row.names=1)
sainfo_data <- read.table(sample_info_file, header=T, sep='\t')

print(pwdist_data)
d <- as.dist(xtabs(pw_distance ~ sample2 + sample1, data=pwdist_data))
pdf(pdf_file, height=10, width=10)

print(d)
mds.cmdscale <- as.data.frame(cmdscale(as.matrix(d), k=2))
print(mds.cmdscale)
mds.cmdscale$names <- rownames(mds.cmdscale)

plot.df <- merge(mds.cmdscale, sainfo_data, by.x="names", by.y="sample_id")
print(plot.df)
ggplot(plot.df, aes(V1, V2, color=sample_depth, label=names)) +
  geom_text(aes(label=names), check_overlap = TRUE, size=2.2,
            hjust = "center", vjust = "bottom", nudge_x = 0) +
  labs(x="", y="", title="") + scale_fill_gradient() + theme_bw()

dev.off()
