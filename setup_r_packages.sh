module load R/3.3.1-foss-2016b

R -e 'source("https://bioconductor.org/biocLite.R");if(!require("AnnotationDbi")) biocLite("AnnotationDbi");if(!require("DESeq2")) biocLite("DESeq2");if(!require("readr")) biocLite("readr");if(!require("TxDb.Athaliana.BioMart.plantsmart28")) biocLite("TxDb.Athaliana.BioMart.plantsmart28");if(!require("tximport")) biocLite("tximport");library("devtools");if(compareVersion(installed.packages()["rmarkdown","Version"],"1.8.10")==-1) install_github("rstudio/rmarkdown")'

module unload R/3.3.1-foss-2016b

