#!/usr/bin/env Rscript


library("knitr")
library("rmarkdown")
args = commandArgs(TRUE)

design_file = args[1]
pvalue = args[2]
session = args[3]

rmarkdown::render("report.Rmd",params=list(design = design_file, p =pvalue, sessionId=session))


