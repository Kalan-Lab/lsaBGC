library(pegas)

args = commandArgs(trailingOnly=TRUE)

# read in loci file data
pos_loci = as.numeric(args[2])
loci.dat <-  read.loci(args[1], loci.sep = "\t", col.pop=(pos_loci+1), col.loci = 2:pos_loci, row.names = 1)

# Calculate Fst
fst.results <- Fst(loci.dat)

# Write results to file
write.table(fst.results, file=args[3], sep='\t')