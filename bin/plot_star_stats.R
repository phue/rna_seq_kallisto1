#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)
star_stat=read.table("star_stats.tab",sep="\t",header=T,stringsAsFactors = F)
star_stat = mutate(star_stat,sample=sapply(strsplit(ID,"Log"),"[[",1) )
star_stat = star_stat %>% select (-c(ID,X))


s2c <- read.table(args[1], header = TRUE, stringsAsFactors=FALSE)
s2c <- mutate(s2c, name = paste(condition,sample,sep="_"))
s2c <-select(s2c, sample = run_accession, condition, name)

star_stat = inner_join(star_stat,s2c)

star_stat = star_stat %>% select (-c(sample,condition))

stat_mat = melt(star_stat, id = c("name"))

pdf("star_stats.pdf")
ggplot(stat_mat,aes(x=name,y=value,fill=name))+ geom_bar(stat="identity")+facet_wrap(~variable,scales = "free")+ guides(fill=FALSE)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())
dev.off()

