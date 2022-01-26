#!/bin/bash

data <- readLines("Egrandis_297_v2.0.geneonlyGT.gff3")

new_data <- sub("(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)", 
"\\9\t\\1\t\\7\t\\4\t\\5", data, perl=TRUE)

new_data2 <- sub("ID=\\S*;Name=","", new_data, perl=TRUE )
write(new_data2, file = "Egrandis_gene_locations_GT.txt")
