#!/usr/bin/env Rscript

##############################################
library(AnnotationDbi)
library(DESeq2)
library(readr)
library(tximport)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(Cairo)

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

run_DESeq=function(dds,dds_noBeta,contrast,cutoff)
{
    res2 = results(dds_noBeta,contrast = contrast)
    res = results(dds,contrast = contrast)
    res$padj[is.na(res$padj)]=1
    res$log2_mle = res2$log2FoldChange
    res = as.data.frame(res)
    res_df = (res[order(res$padj<cutoff,abs(res$log2FoldChange),decreasing = T),])
    return(res_df)
}

mybarplot = function(run, p){
    sigs = dplyr::filter(as.data.frame(run),padj<p)
    counts = dplyr::transmute(sigs, log_0 = sign(log2FoldChange) ,log_1=ifelse(abs(log2FoldChange)>1,1,0)*sign(log2FoldChange) , log_2= ifelse(abs(log2FoldChange)>2,1,0)*sign(log2FoldChange))
    stats = cbind(apply(counts,2,function(x) y = sum(x==1)), apply(counts,2,function(x) y = sum(x==0)), apply(counts,2,function(x) y = sum(x==-1)))
   colnames(stats)=c("up","none","down")
   barplot(t(stats[,-2]),col=c("yellow","blue"),names.arg=c("abs(log2)>0","abs(log2)>1","abs(log2)>2" ),ylab="number significant genes")
}


myplotMA=function(dds,contrast,p){
  res = results(dds,contrast = contrast)
  plotMA(res,main=contrast,ylim=c(-2,2),alpha=p)
}
add_norm_counts=function(dds,contrast,res){
  n_counts=counts(dds,normalized=T)
  g = grep(contrast[1],colnames(colData(dds)))
  g1 = grep(contrast[2],colData(dds)[,g])
  g2 = grep(contrast[3],colData(dds)[,g])
  mean1 = rowMeans(n_counts[,g1])
  mean2 = rowMeans(n_counts[,g2])
  numb = cbind(mean1,mean2,n_counts[,g1],n_counts[,g2])
  colnames(numb)[1:2]=paste0("mean",contrast[2:3])
  res=as.data.frame(res)
  new=cbind(res,numb[rownames(res),])
  new
}
clean_up_df=function(res){
  remove=c("baseMean","lfcSE","stat","pvalue")
  new = res[,!colnames(res)%in%remove]
  new
}

write.resfile=function(res,filename){
  write.csv(res,file=filename)
}



#################################################
### set up 

args <- commandArgs(TRUE)
sample_id <- dir(args[1])
kal_dirs <- sapply(sample_id, function(id) file.path(args[1], id))

s2c <- read.table(args[2], header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, name = paste(condition,sample,sep="_"))
s2c <- dplyr::select(s2c, sample = run_accession, condition, name)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c <- dplyr::mutate(s2c, file = paste0(kal_dirs,"/abundance.tsv"))
s2c <- s2c[order(s2c$condition), ]


txdb_choice = args[5]
if(txdb_choice == "tair10"){
	txdb = TxDb.Athaliana.BioMart.plantsmart28
} else { 
txdb = loadDb(args[5])
}
keys <- keys(txdb)
df = select(txdb, keys=keys,columns=c("TXCHROM", "TXSTART", "TXEND","TXNAME","TXSTRAND"), keytype="GENEID")
tx2gene=df[,2:1]

txi <- tximport(s2c$file, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
counts=txi$counts
colnames(counts)=s2c$sample
countsToUse = round(counts)


pval <- as.numeric(args[4])


##################################################
### deseq2 - first part 

colD=data.frame(group=s2c$condition,name=s2c$name)
dds = DESeqDataSetFromMatrix(countsToUse,colData=colD,design=~group)
dds_noBeta=DESeq(dds,betaPrior=FALSE)
dds=DESeq(dds)
save(dds,dds_noBeta,file="dds.Rdata")


### initial plots
ntd <- normTransform(dds)
CairoPNG("pairs.png")
pairs(assay(ntd),pch=".",labels = colData(dds)$group,upper.panel = panel.cor)
dev.off()

rld <- rlog(dds, blind=FALSE)
CairoPNG("pca.png")
plotPCA(rld, intgroup=c("name"))
dev.off()


### analysis + MA plots
co = read.table(args[3],header=T)
runs=list()
for ( i in 1:nrow(co)){
  cont=c("group",colnames(co)[c(which(co[i,]==1),which(co[i,]==-1))])
runs[[i]]=run_DESeq(dds,dds_noBeta,contrast=cont,cutoff=pval) 
 runs[[i]]=add_norm_counts(dds,cont,runs[[i]])
runs[[i]]=clean_up_df(runs[[i]])
  CairoPNG(paste(paste("maplot",paste(cont,collapse="_"),sep="_"),"png",sep=".")) 
  myplotMA(dds,cont,p=pval)
  dev.off()
  CairoPNG(paste(paste("barplot",paste(cont,collapse="_"),sep="_"),"png",sep="."))
  mybarplot(runs[[i]],p=pval)
  dev.off()  
  write.resfile(runs[[i]], paste(paste("contrast",paste(cont,collapse="_"),sep="_"),"csv",sep="."))
}

writeLines(capture.output(sessionInfo()), "sessionInfo_deseq2.txt")

