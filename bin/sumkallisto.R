#!/usr/bin/env Rscript

##############################################
library(readr)
library(tximport)


args <- commandArgs(TRUE)

sample_id <-  dir(args[1])
sample_dir <- args[1]
txtype <- args[2]
design <- args[3]

if (txtype == "tair10"){
	library(TxDb.Athaliana.BioMart.plantsmart28)
	txdb = TxDb.Athaliana.BioMart.plantsmart28
} else{
	library(AnnotationDbi)
	txdb = loadDb(txtype)
}

kal_dirs <- sapply(sample_id, function(id) file.path(sample_dir, id))
print(kal_dirs)
s2c <- read.table(design, header = TRUE, stringsAsFactors=FALSE)
print(s2c)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
print(s2c)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)
s2c <- dplyr::mutate(s2c, file = paste0(kal_dirs,"/abundance.tsv"))
print(s2c)
s2c <- s2c[order(s2c$condition), ]

keys <- keys(txdb)
print(head(keys))

df = select(txdb, keys=keys,columns=c("TXCHROM", "TXSTART", "TXEND","TXNAME","TXSTRAND"), keytype="GENEID")
tx2gene=df[,2:1]

txi <- tximport(s2c$file, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
counts=txi$counts
colnames(counts)=s2c$sample
countsToUse = round(counts)
colD=data.frame(group=s2c$condition)

save(countsToUse,colD,file="kallistoData.Rdata")
write.table(countsToUse, file = "kallisto_counts.tab")


