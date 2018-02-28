#!/usr/bin/env Rscript


library("knitr")
library("rmarkdown")
args = commandArgs(TRUE)

design_file = args[1]
pvalue = args[2]
contrast = args[3]
session = args[4]

rmarkdown::render("report.Rmd",params=list(design = design_file, p =pvalue,contrast_file = contrast, sessionId=session))


