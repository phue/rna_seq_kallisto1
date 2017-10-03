#!/usr/bin/env Rscript

library(rtracklayer)
library(AnnotationDbi)

args <- commandArgs(TRUE)

db=args[2]
print(db)
mydb=loadDb(db)

file=args[1]
print(file)
wigToBigWig(x=file,seqinfo=seqinfo(mydb) )


