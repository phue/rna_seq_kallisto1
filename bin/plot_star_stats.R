#!/usr/bin/env Rscript

star_stat=read.table("star_stats.tab",sep="\t",header=T,stringsAsFactors = F)

star_mat = star_stat[,-1]
rownames(star_mat) = sapply(strsplit(star_stat[,1],"_"),"[[",1) 
star_mat = star_mat[,!is.na(colSums(star_mat))]

pdf("star_stats.pdf")
for (i in 1:ncol(star_mat)){
  barplot(star_mat[,i],names.arg=rownames(star_mat),las=2,main=colnames(star_mat)[i])
}
dev.off()
