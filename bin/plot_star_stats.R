#!/usr/bin/env Rscript


args <- commandArgs(TRUE)
star_stat=read.table("star_stats.tab",sep="\t",header=T,stringsAsFactors = F)
star_stat = dplyr::mutate(star_stat,ID=sapply(strsplit(ID,"Log"),"[[",1) )



s2c <- read.table(args[1], header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, name = paste(condition,sample,sep="_"))
s2c <- dplyr::select(s2c, sample = run_accession, condition, name)

#star_stat dplyr::mutate(star_stat,name = sc2$name[s2c$sample==])



#star_mat = star_stat[,-1]
#rownames(star_mat) = sapply(strsplit(star_stat[,1],"_"),"[[",1) 
#star_mat = star_mat[,!is.na(colSums(star_mat))]

pdf("star_stats.pdf")
for (i in 1:ncol(star_mat)){
  barplot(star_mat[,i],names.arg=rownames(star_mat),las=2,main=colnames(star_mat)[i])
}
dev.off()
