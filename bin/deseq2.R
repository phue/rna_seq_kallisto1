#!/usr/bin/env Rscript

##############################################
library(AnnotationDbi)
library(DESeq2)
library(readr)
library(tximport)
library(TxDb.Athaliana.BioMart.plantsmart28)

###############################################
## functions


# function that adds correlation to pairs plot. Used with or without contrast file!
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

# function that runs DE analysis for one contrast (mle and map log2, p-value, padj...) # only when contrast file is provided
run_DESeq=function(dds,contrast,cutoff)
{
    res2 = results(dds,contrast = contrast)##MLE
    res = lfcShrink(dds,contrast = contrast)##MAP
    res$padj[is.na(res$padj)]=1
    res$log2_mle = res2$log2FoldChange
    res = as.data.frame(res)
    res_df = (res[order(res$padj<cutoff,abs(res$log2FoldChange),decreasing = T),])
    return(res_df)
}

# function to plot barplot with nr significant genes. # only when contrast file is provided
mybarplot = function(run, p){
    sigs = dplyr::filter(as.data.frame(run),padj<p)
    counts = dplyr::transmute(sigs, log_0 = sign(log2FoldChange) ,log_1=ifelse(abs(log2FoldChange)>1,1,0)*sign(log2FoldChange) , log_2= ifelse(abs(log2FoldChange)>2,1,0)*sign(log2FoldChange))
    stats = cbind(apply(counts,2,function(x) y = sum(x==1)), apply(counts,2,function(x) y = sum(x==0)), apply(counts,2,function(x) y = sum(x==-1)))
   colnames(stats)=c("up","none","down")
   barplot(t(stats[,-2]),col=c("yellow","blue"),names.arg=c("abs(log2)>0","abs(log2)>1","abs(log2)>2" ),ylab="number significant genes")
}


# function to plot maplot. # only when contrast file is provided
myplotMA=function(dds,contrast,p){
  res = lfcShrink(dds,contrast = contrast)
  plotMA(res,main=contrast,ylim=c(-2,2),alpha=p)
}

# function to add normalized counts and mean to table. Use with or without contrast
add_norm_counts=function(dds,contrast,res){
    n_counts=counts(dds,normalized=T)
    g = grep(paste0("^",contrast[1],"$"),colnames(colData(dds))) # column with group
    means = list()
    counts = list()
    for(i in 2:length(contrast)){
        cond_columns = grep(paste0("^",contrast[i],"$"),colData(dds)[,g])
        counts[[i-1]] = n_counts[,cond_columns]
        means[[i-1]] = rowMeans(counts[[i-1]])
  }
  a_mean = do.call("cbind",means)
  colnames(a_mean) = paste0("mean",unique(contrast[2:length(contrast)]))
  a_counts = do.call("cbind",counts)
  tot = cbind(a_mean, a_counts)
  colnames(tot) =paste("norm.counts",  colnames(tot),sep="_")
  if(res!=FALSE){
      res=as.data.frame(res)
      tot=cbind(res,tot[rownames(res),])
  }
  tot
}

# function to add tpm and mean to table. Re-write using groups..?
add_mean_tpm=function( dds, tpm, contrast){
    g = grep(paste0("^",contrast[1],"$"),colnames(colData(dds))) # column with group
    means = list()
    tpms = list()
    for(i in 2:length(contrast)){
        cond_columns = grep(paste0("^",contrast[i],"$"),colData(dds)[,g])
        cond_names = rownames(colData(dds))[cond_columns]
        tpms[[i-1]] = tpm[,cond_names]
        means[[i-1]] = rowMeans(tpms[[i-1]])
    }
    a_mean = do.call("cbind",means)
    colnames(a_mean) = paste0("mean",unique(contrast[2:length(contrast)]))
    a_tpm = do.call("cbind",tpms)
    tot = cbind(a_mean, a_tpm)
    colnames(tot) =paste("tpm",  colnames(tot),sep="_")
    tot
}

clean_up_df=function(res){
  remove=c("baseMean","lfcSE","stat")
  new = res[,!colnames(res)%in%remove]
  new
}




#################################################
### set up 

args <- commandArgs(TRUE)
## explicitly setting the arguments to make it easier to read/edit/debug
dirname <- args[1]
design_file <- args[2]
contrast_file <- args[3]
pval <- as.numeric(args[4])
txdb_choice <- args[5]
n_filter <-as.numeric(args[6])
sessID <-args[7]

sample_id <- dir(dirname)
kal_dirs <- sapply(sample_id, function(id) file.path(dirname, id))

s2c <- read.table(design_file, header = TRUE, stringsAsFactors=FALSE,sep=",",colClasses=c("character","character","character"))
mn_condition = make.names(s2c$condition)
if (length(unique(mn_condition))!=length(unique(s2c$condition)))
	stop("Condition names contain characters that can not be resolved,please use only numbers, letteers, dot and underscore")
s2c$condition = mn_condition
s2c <- dplyr::mutate(s2c, name = paste(condition,sample,sep="_"))
s2c <- dplyr::select(s2c, sample = run_accession, condition, name)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c <- dplyr::mutate(s2c, file = paste0(kal_dirs,"/abundance.h5"))
s2c <- s2c[order(s2c$condition), ]



if(txdb_choice == "tair10"){
	txdb = TxDb.Athaliana.BioMart.plantsmart28
} else { 
txdb = loadDb(txdb_choice)
}
keys <- keys(txdb)
df = select(txdb, keys=keys,columns=c("TXCHROM", "TXSTART", "TXEND","TXNAME","TXSTRAND"), keytype="GENEID")
tx2gene=df[,2:1]

txi <- tximport(s2c$file, type = "kallisto", tx2gene = tx2gene)

counts=txi$counts
tpm = txi$abundance
colnames(counts)=s2c$sample
colnames(tpm) = s2c$sample
countsToUse = round(counts)





## tpms 


##################################################
### deseq2 - first part, run always

colD=data.frame(group=s2c$condition,name=s2c$name)
if (length(unique(colD$group))==1){
        dds = DESeqDataSetFromMatrix(countsToUse,colData=colD,design=~1)
} else {
dds = DESeqDataSetFromMatrix(countsToUse,colData=colD,design=~group)
}

## in case want to filter out genes with few reads
if(n_filter>0){ # does not make sense to do if no contrast! (hence design=~group)
    n_counts = counts(dds)
    n_counts = subset(n_counts, apply(n_counts,1,max)>=n_filter)
    dds = DESeqDataSetFromMatrix(n_counts,colData=colD,design=~group)
}

dds=DESeq(dds)
save(dds,file="dds.Rdata")


### initial plots
ntd <- normTransform(dds)
png("pairs.png")
pairs(assay(ntd),pch=".",labels = colData(dds)$group,upper.panel = panel.cor)
dev.off()

rld <- rlog(dds, blind=FALSE)
png("pca.png")
plotPCA(rld, intgroup=c("name"))
dev.off()


### second part, only if contrast file
### analysis + MA plots
if(contrast_file!='NULL'){ # NULL means no contrasts should be analyized
    contrast_to_run = read.table(contrast_file,header=T,sep=",")
    runs=list()
    for ( i in 1:nrow(contrast_to_run)){
        my_contrast=c("group",colnames(contrast_to_run)[c(which(contrast_to_run[i,]==1),which(contrast_to_run[i,]==-1))])
        runs[[i]]=run_DESeq(dds,contrast=my_contrast,cutoff=pval)
        tpm_df = add_mean_tpm(dds, tpm, my_contrast)
        runs[[i]]=add_norm_counts(dds,my_contrast,runs[[i]])
        runs[[i]] =cbind(runs[[i]],tpm_df[rownames(runs[[i]]),])
        runs[[i]]=clean_up_df(runs[[i]])
        png(paste(paste("maplot",paste(my_contrast,collapse="__"),sep=""),"png",sep="."))
        myplotMA(dds,my_contrast,p=pval)
        dev.off()
        png(paste(paste("barplot",paste(my_contrast,collapse="__"),sep=""),"png",sep="."))
        mybarplot(runs[[i]],p=pval)
        dev.off()
        com = paste("#",sessID)
        tab=runs[[i]]
        tab=cbind("GeneId"=rownames(tab),tab)
        file=paste(paste("contrast",paste(my_contrast,collapse="__"),sep="_"),"csv",sep=".")
        write.table(com, file = file,sep=",",quote = FALSE,row.names=FALSE,col.names=FALSE)
        write.table(tab, file = file, append=T, sep="," , quote=FALSE,row.names=FALSE)
    }
} else {
    my_contrast = c("group",as.character(unique(colD$group)))
    com = paste("#",sessID)
    tpm_df = add_mean_tpm(dds, tpm, my_contrast)
    runs=add_norm_counts(dds,my_contrast,FALSE)
    runs = cbind(runs,tpm_df[rownames(runs),])
    runs=cbind("GeneId"=rownames(runs),runs)
    write.table(com, file = "table.csv",sep=",",quote = FALSE,row.names=FALSE,col.names=FALSE)
    write.table(runs,file="table.csv",append=T,sep=",",quote = FALSE,row.names=FALSE)
}


writeLines(capture.output(sessionInfo()), "sessionInfo_deseq2.txt")
writeLines(print(args),"Rarguments.txt")

