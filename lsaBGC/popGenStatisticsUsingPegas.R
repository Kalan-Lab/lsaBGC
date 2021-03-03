library(ape)
library(pegas)
library(adegenet)

args = commandArgs(trailingOnly=TRUE)

# read in codon alignment and store as a DNAbin object
codon.msa <- fasta2DNAbin(args[1], quiet=T, snpOnly=F)

# read in loci file data
loci.dat <-

# domain  domain_index    min_pos max_pos
positions.dat <- read.table(args[2], header=T, sep='\t')
# pos     num_seqs        num_alleles     num_gaps        maj_allele_freq
popgen.dat <- read.table(args[3], header=T, sep='\t')
# pos     type    effective
