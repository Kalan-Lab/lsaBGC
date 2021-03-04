library(ape)
library(pegas)
library(adegenet)

args = commandArgs(trailingOnly=TRUE)

# read in codon alignment and store as a DNAbin object
codon.msa <- fasta2DNAbin(args[1], quiet=T, snpOnly=F)

# Calculate Tajima's D
tajima.results <- tajima.test(codon.msa)

# Write results to file
write.table(tajima.results, file=args[2], sep='\n')