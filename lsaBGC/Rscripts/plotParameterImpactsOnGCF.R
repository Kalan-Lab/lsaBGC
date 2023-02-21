library(scatterpie)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(cowplot)
args = commandArgs(trailingOnly=TRUE)

# read in plotting inputs
data <- read.table(args[1], header=T, sep='\t')
annot_classes <- scan(args[2], what="", sep="\n")
largest_gcfs_data <- read.table(args[3], header=T, sep='\t')
overview_data <- read.table(args[4], header=T, sep='\t')
alluvian_data <- read.table(args[5], header=T, sep='\t')

# define color scheme for plots
nb.cols <- length(annot_classes)-1
annot_colors <- colorRampPalette(brewer.pal(nb.cols, "Set2"))(nb.cols)
annot_colors <- c("#808080", annot_colors)
names(annot_colors) <- annot_classes
core_colors <- c("#737574", "#222422", "#ed665f")
names(core_colors) <- c("Not Applicable", "Core Exists", "Core DNE")

# open pdf for the report
pdf(args[6], height=8, width=11)

# Title page + overview scattermap
repname <- ggplot()+theme_void() + ggtitle("Parameter Impact on\nGCF Delineations") + theme(plot.title = element_text(hjust = 0.5, size=70, face='bold'))
scatter_plot <- ggplot(overview_data, aes(x=NumberClusters, y=SingletonToClustered)) +
  geom_jitter(aes(color=as.factor(Inflation), shape=as.factor(JaccardSim)), alpha=0.7, size=2) + theme_classic() +
  ggtitle("Tradeoff between Cluster Granularity vs. Singletons Incurred") + xlab("Number of GCFs") +
  ylab("Ratio of BGCs which are\nsingletons to those in GCFs") +
  guides(shape=guide_legend(title="Jaccard Similarity Cutoff"), color=guide_legend(title="MCL Inflation"))+
  theme(legend.position = "bottom")
print(plot_grid(repname, scatter_plot, rel_heights=c(1,2), ncol=1))

# Overview Heatmaps
gheat1 <- ggplot(overview_data, aes(x=as.factor(Inflation), y=as.factor(JaccardSim), fill=SingletonToClustered)) +
  geom_tile() + theme_classic() +  scale_fill_gradient(low = "lightblue", high = "darkblue") + xlab("MCL Inflation") +
  ylab("Jaccard Similarity Cutoff") + ggtitle("Ratio of BGCs which are\nsingletons to those in GCFs")+ theme(legend.title = element_blank())
gheat2 <- ggplot(overview_data, aes(x=as.factor(Inflation), y=as.factor(JaccardSim), fill=NumberClusters)) +
  geom_tile() + theme_classic() +   scale_fill_gradient(low = "lightgreen", high = "darkgreen") + xlab("MCL Inflation") +
  ylab("Jaccard Similarity Cutoff") + ggtitle("Number of GCFs")+ theme(legend.title = element_blank())
gheat3 <- ggplot(overview_data, aes(x=as.factor(Inflation), y=as.factor(JaccardSim), fill=MixedAnnotationProportion)) +
  geom_tile() + theme_classic() +   scale_fill_gradient(low = "yellow", high = "gold") + xlab("MCL Inflation") +
  ylab("Jaccard Similarity Cutoff") + ggtitle("Proportion of GCFs which feature\nmultiple BGC annotation categories\n(including unannotated)")+ theme(legend.title = element_blank())
gheat4 <- ggplot(overview_data, aes(x=as.factor(Inflation), y=as.factor(JaccardSim), fill=MixedAnnotationWoHypotheticalProportion)) +
  geom_tile() + theme_classic() +   scale_fill_gradient(low = "pink", high = "darkred") + xlab("MCL Inflation") +
  ylab("Jaccard Similarity Cutoff") + ggtitle("Proportion of GCFs which feature\nmultiple BGC annotation categories\n(not including unannotated)")+ theme(legend.title = element_blank())

gheatmaps_top <- plot_grid(gheat1, gheat2)
gheatmaps_bottom <- plot_grid(gheat3, gheat4)
print(plot_grid(gheatmaps_top, gheatmaps_bottom, rel_heights=c(1,1), ncol=1))

# Boxplot views
boxpl1 <- ggplot(data, aes(x=as.factor(Inflation), y=CoreSize, fill=as.factor(Inflation))) + geom_boxplot(show.legend=F) + facet_wrap(~JaccardSim, ncol=length(data$JaccardSim)) + theme_classic() + xlab("MCL Inflation") + ylab("") + ggtitle("Number of Core Homolog Groups per GCF") + scale_fill_brewer(palette="Set2")# + scale_y_log10()
boxpl2 <- ggplot(data, aes(x=as.factor(Inflation), y=Samples, fill=as.factor(Inflation))) + geom_boxplot(show.legend=F) + facet_wrap(~JaccardSim, ncol=length(data$JaccardSim)) + theme_classic() + xlab("MCL Inflation") + ylab("") + ggtitle("Number of Samples with GCF") + scale_fill_brewer(palette="Set2") #+ scale_y_log10()
boxpl3 <- ggplot(data, aes(x=as.factor(Inflation), y=AvgGeneCount, fill=as.factor(Inflation))) + geom_boxplot(show.legend=F) + facet_wrap(~JaccardSim, ncol=length(data$JaccardSim)) + theme_classic() + xlab("MCL Inflation") + ylab("") + ggtitle("Avg. Gene Count per BGC in GCF") + scale_fill_brewer(palette="Set2") #+ scale_y_log10()
boxpl4 <- ggplot(data, aes(x=as.factor(Inflation), y=StdDevGeneCount, fill=as.factor(Inflation))) + geom_boxplot(show.legend=F) + facet_wrap(~JaccardSim, ncol=length(data$JaccardSim)) + theme_classic() + xlab("MCL Inflation") + ylab("") + ggtitle("Std. Deviation of Gene Count for BGCs in GCF") + scale_fill_brewer(palette="Set2") #+ scale_y_log10()

top_row <- plot_grid(boxpl1, boxpl2)
bot_row <- plot_grid(boxpl3, boxpl4)

print(plot_grid(boxpl1, boxpl2, boxpl3, boxpl4, rel_heights=c(1,1,1,1), ncol=1))

# Alluvian/Sankey diagrams
colors <- c("#000000", "#FF0000")
names(colors) <- c("cluster", "singleton")

alluvian_plot <- ggplot(alluvian_data, aes(y = PassageWeight, axis1 = reorder(Parameter_Source, -Source_Size), axis2 = reorder(Parameter_Dest, -Dest_Size))) +
  geom_flow(fill='darkgrey', aes(color = Parameter_Source), show.legend=F) +
  geom_stratum(aes(fill=Source_Type), color='white', show.legend=F)  + facet_grid(JaccardSim ~ Source_Inflation) +
  xlab("MCL Inflation") + ylab("Jaccard Similarity Cutoff") + scale_fill_manual(values=colors) +
  theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
print(alluvian_plot)

# determine number of parameter combinations tested and axes boundaries for scatter plots
params <- unique(data$Parameters)
npages <- length(params)
xmin <- min(data[!is.na(data$Samples_Offset),]$Samples_Offset)
ymin <- min(data[!is.na(data$MultiBGCSamples),]$MultiBGCSamples)
xmax <- max(data[!is.na(data$Samples_Offset),]$Samples_Offset)
ymax <- max(data[!is.na(data$MultiBGCSamples),]$MultiBGCSamples)

for(i in 1:npages) {
  curr_param = params[i]

  data_filt <- data[data$Parameters == curr_param,]
  data_filt_scatter <- data_filt[!is.na(data_filt$Samples_Offset),]
  largest_gcfs_data_filt <- largest_gcfs_data[largest_gcfs_data$Parameters == curr_param,]

  paramname <- ggplot()+theme_void() + ggtitle(curr_param) + theme(plot.title = element_text(hjust = 0.5, size=32, face='bold'))

  gscat <- ggplot()  +
    scale_fill_manual(values=annot_colors) + theme_classic() +
    xlab("Number of Samples with GCF") + ylab("Number of Samples with\n>= 2 BGCs for GCF") +
    geom_point(aes(x=Samples_Offset, y=MultiBGCSamples_Offset, shape=as.factor(CoreExists)), data=data_filt_scatter, size=1.3, show.legend=F, alpha=0.5) +
    geom_scatterpie(aes(x=Samples_Offset, y=MultiBGCSamples_Offset), data=data_filt_scatter, alpha=0.4, cols=annot_classes, show.legend=F) +
    coord_fixed(xlim=c(xmin, xmax), ylim=c(ymin, ymax), expand=T)


  ghist_1 <- ggplot(data_filt, aes(x=Samples, fill=CoreExists)) + geom_histogram(color="black", size=0.3) + theme_classic() +
    theme(legend.position = c(0.8, 0.8)) + scale_fill_manual(values=core_colors) + xlab("Sample Count") + guides(fill=guide_legend(title=""))

  ghist_2 <- ggplot(data_filt, aes(x=StdDevGeneCount)) + xlab("Std. Deviation in Gene Count") + geom_histogram(color="black", fill="black") + theme_classic() + scale_y_log10()

  gbar <- ggplot(largest_gcfs_data_filt, aes(x=reorder(GCF, Total.BGCs), y=Count, fill=Annotation)) +
    geom_bar(stat='identity', color='black', show.legend=F) + theme_classic() + scale_fill_manual(values=annot_colors) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("Top 50 GCFs with Most BGCs")

  mid_row <- plot_grid(ghist_1, ghist_2)
  print(plot_grid(paramname, mid_row, gscat, gbar, rel_heights=c(1,2,4,2), ncol=1))
}

dev.off()
