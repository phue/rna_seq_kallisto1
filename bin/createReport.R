#!/usr/bin/env Rscript


library("knitr")
library("rmarkdown")
args = commandArgs(TRUE)

kall_folder = args[1]
design_file = args[2]
pvalue = args[3]
stats = args[4]
report = args[5]
param = args[6]

file.copy(report,"report.Rmd")
rmarkdown::render("report.Rmd",params=list(kal_folder=kall_folder,design = design_file, p =pvalue,stat=stats, para = param))


