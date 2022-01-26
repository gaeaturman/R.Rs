#!/bin/env Rscript

handleFile = function(filepath) {
  curLine = 0
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    curLine = curLine + 1
    if ( length(line) == 0 ) {
    	break
    }
    Ggrep=grep("[AT]GATA[GA]", line)
    lenGrep=length(Ggrep)

    Ggrepx=gregexpr("[AT]GATA[GA]", line)
    lenGrepx=length(Ggrepx)

    #checks for if no grep/grepgexpr matches, if so pos == N/A	

    print(paste("--------------------------------------"))
    print(paste("For Gene ", curLine, ", result of length(grep) is ", lenGrep))
    print(paste("The output of grep() is ", Ggrep))    
    print(paste("---------------------------------------------------------------------"))
    print(paste("For Gene ", curLine, ", result of length(gregexpr) is", lenGrepx))
    print(paste("The output of gregexpr is", Ggrepx))
    print(paste("--------------------------------------"))
  }
  close(con)
}

aFile='grape_promoters.txt'
handleFile(aFile)
