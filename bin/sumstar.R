#!/usr/bin/env Rscript

##############################################

args <- commandArgs(TRUE)


sample_id <-  dir(args[1])
sample_dir <- args[1]
strand <- args[2]
design <- args[3]


if (strand =="fr-stranded"){
	st=3-1 # first column id
} else if (strand == "rf-stranded"){
	st=4-1
} else {
	st=2-1
}
print(strand)
print(st)
star_dirs <- sapply(sample_id, function(id) file.path(sample_dir, id))
print(star_dirs)
s2c <- read.table(design, header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = star_dirs)
s2c <- s2c[order(s2c$condition), ]
print(s2c)

#sample_table = read.table(design)
counts=list()
for (i in s2c$path){
	tmp = read.table(i,row.names=1)
	print(dim(tmp))
	counts[[i]]=tmp[,st,drop=F]
}
mat = do.call("cbind",counts)
colnames(mat)=s2c$sample
torm=grep("N_",rownames(mat))
countsToUse = mat  # [-torm,] keep for stats
colD = data.frame("group"=s2c$condition)

save(countsToUse,colD,file="starData.Rdata")
write.table(countsToUse, file = "star_counts.tab") 

