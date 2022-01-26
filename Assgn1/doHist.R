#!/bin/env Rscript

data=readLines('geneLens.txt')
data=as.numeric(data)
pdf("Egrandis_gene_lengths_GT.pdf") 
hist(data,main="Gene Sequence Lengths",xlab = "Length(bp)",ylab="Number of Genes",col = "yellow",border 
= "blue")

invisible(dev.off())

