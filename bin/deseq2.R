#!/usr/bin/env Rscript

##############################################

library(DESeq2)
library(readr)
library(tximport)
library(TxDb.Athaliana.BioMart.plantsmart28)

###############################################
## functions

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
run_DESeq=function(dds,dds_noBeta,contrast,cutoff=0.1)
{
    res2 = results(dds_noBeta,contrast = contrast)
    res = results(dds,contrast = contrast)
    res$padj[is.na(res$padj)]=1
    res$log2_mle = res2$log2FoldChange
    res_df = (res[order(res$padj<cutoff,abs(res$log2FoldChange),decreasing = T),])
    return(res_df)
}

#################################################

args <- commandArgs(TRUE)
sample_id <- dir(args[1])
kal_dirs <- sapply(sample_id, function(id) file.path(args[1], id))

s2c <- read.table(args[2], header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c <- dplyr::mutate(s2c, file = paste0(kal_dirs,"/abundance.tsv"))
s2c <- s2c[order(s2c$condition), ]

txdb = TxDb.Athaliana.BioMart.plantsmart28
keys <- keys(txdb)
df = select(txdb, keys=keys,columns=c("TXCHROM", "TXSTART", "TXEND","TXNAME","TXSTRAND"), keytype="GENEID")
tx2gene=df[,2:1]

txi <- tximport(s2c$file, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
counts=txi$counts
colnames(counts)=s2c$sample
countsToUse = round(counts)

colD=data.frame(group=s2c$condition)
dds = DESeqDataSetFromMatrix(countsToUse,colData=colD,design~group)

#dds = DESeq(dds)
ntd <- normTransform(dds)
pdf("pairs.pdf")
pairs(assay(ntd),pch=".",labels = colData(dds)$group,upper.panel = panel.cor)
dev.off()


### PCA plot, pairs, MA

