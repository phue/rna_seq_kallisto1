#!/usr/bin/env Rscript


library("knitr")
library("rmarkdown")
args = commandArgs(TRUE)
print(args)

kall_folder = args[1]
design_file = args[2]
pvalue = args[3]
stats = args[4]
contrast = args[5]

rmarkdown::render("report.Rmd",params=list(kal_folder=kall_folder,design = design_file, p =pvalue,stat=stats,contrast_file = contrast, sessionId=args[6]))
#,knit_root_dir = args[7], output_dir = args[7], intermediates_dir = args[7])


