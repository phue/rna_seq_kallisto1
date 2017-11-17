#!/usr/bin/env Rscript

##############################################

args <- commandArgs(TRUE)
sample_id <-  dir(args[1]) ## all "XReadsPerGene.out.tab" files that exist in collected channel 
sample_dir <- args[1] ## folder name of collected channel
strand <- args[2] ## strand specificity in samples
design <- args[3] ## file with conditions, samples etc

## set which column to read in XReadsPerGene.out.tab (depends on strand)

if (strand =="fr-stranded"){
	st=3-1 # first column id
} else if (strand == "rf-stranded"){
	st=4-1 # first column id
} else {
	st=2-1 # first column id
}

## whole path to XReadsPerGene.out.tab files
star_dirs <- sapply(sample_id, function(id) file.path(sample_dir, id))

## generate sample 2 condition table for design 
s2c <- read.table(design, header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = paste0(paste(sample_dir,sample,sep="/"),"ReadsPerGene.out.tab"))
## check that the files in the path column exist
if (length(setdiff(s2c$path,star_dirs))>0) stop('samples and folders do not match')
s2c <- s2c[order(s2c$condition), ]

## read in and combine into one matrix and colData
counts=list()
for (i in s2c$path){
	tmp = read.table(i,row.names=1)
	counts[[i]]=tmp[,st,drop=F]
}
mat = do.call("cbind",counts)
mat = as.matrix(mat)
colnames(mat)=s2c$sample
countsToUse = mat  
colD = data.frame("group"=s2c$condition)

## write and save objects
save(countsToUse,colD,file="starData.Rdata")
write.table(countsToUse, file = "star_counts.tab") 

