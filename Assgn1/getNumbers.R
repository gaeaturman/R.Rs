#!/bin/env Rscript

dataVector=readLines("Egrandis_genelocationsGT.txt")

onPosStrand=grep("\\+", dataVector)
noMatches=length(onPosStrand)
print(paste(noMatches))


justStartLoc=sub("(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)", "\\4", dataVector, perl=TRUE)
justEndLoc=sub("(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)", "\\5", dataVector, perl=TRUE)

justLengths=as.numeric(justEndLoc)-as.numeric(justStartLoc)
#putting data into lengths file to make hist later
numCol=1
numRow= length(justLengths)/numCol

matrixLens=matrix(as.numeric(justLengths), numRow, numCol)
write.table(matrixLens, "geneLens.txt", row.names = FALSE, col.names=FALSE)

meanLen=mean(justLengths)
print(paste(meanLen))

shortestLen=min(justLengths)
print(paste(shortestLen))

greatestLen=max(justLengths)
print(paste(greatestLen))


