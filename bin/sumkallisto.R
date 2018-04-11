#!/usr/bin/env Rscript

##############################################
library(readr)
library(tximport)


args <- commandArgs(TRUE) 
sample_id <-  dir(args[1]) ## all kallisto output folders  that exist in collected channel 
sample_dir <- args[1] ## folder name of collected channel
txtype <- args[2] ## which tx db should be used
design <- args[3] ## file with conditions, samples etc

## load appropriate package and txdb 
if (txtype == "tair10"){
	library(TxDb.Athaliana.BioMart.plantsmart28)
	txdb = TxDb.Athaliana.BioMart.plantsmart28
} else{
	library(AnnotationDbi)
	txdb = loadDb(txtype)
}

## whole path to kallisto output folders 
kal_dirs <- sapply(sample_id, function(id) file.path(sample_dir, id))
print(kal_dirs)

## generate sample 2 condition table for design 
s2c <- read.table(design, header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
## path generated from sample in table
s2c <- dplyr::mutate(s2c, path = paste(sample_dir,sample,sep="/kallisto_"))
## check that the files in the path column exist
if (length(setdiff(s2c$path,kal_dirs))>0) stop('samples and folders do not match')
## files to read 
s2c <- dplyr::mutate(s2c, file = paste0(path,"/abundance.h5"))
s2c <- s2c[order(s2c$condition), ]

### do the tximport thing 
keys <- keys(txdb)
print(head(keys))
df = select(txdb, keys=keys,columns=c("TXCHROM", "TXSTART", "TXEND","TXNAME","TXSTRAND"), keytype="GENEID")
tx2gene=df[,2:1]
txi <- tximport(s2c$file, type = "kallisto", tx2gene = tx2gene )
counts=txi$counts
colnames(counts)=s2c$sample
countsToUse = round(counts)
colD=data.frame(group=s2c$condition)

## write and save
save(countsToUse,colD,file="kallistoData.Rdata")
write.table(countsToUse, file = "kallisto_counts.tab")


