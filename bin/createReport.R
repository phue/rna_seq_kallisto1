#!/usr/bin/env Rscript


library("knitr")
library("rmarkdown")
args = commandArgs(TRUE)
print(args)

kall_folder = args[1]
design_file = args[2]
pvalue = args[3]
stats = args[4]
report = args[5]
param = args[6]

file.copy(args[9],"Child.Rmd")
rmarkdown::render(report,params=list(kal_folder=kall_folder,design = design_file, p =pvalue,stat=stats, para = param,contrast_file = args[8]),knit_root_dir = args[7], output_dir = args[7], intermediates_dir = args[7])


