library(data.table)
library(tidyverse)
library(zip)
NAFLD <- fread("/home/couchr02/Mendel_Commun/Christian/GWAS/40_diseases/meta_clean_NAFLD.txt")
# Remove all variants on chromosome 6 in the region 26MB to 34MB (the MHC region). chr6:28,477,797-33,448,354
NAFLD <- NAFLD[!(CHR == 6 & as.numeric(POS) > 28000000 & as.numeric(POS) < 34000000), ]
colnoms <- c("snpid",	"A1", "A2", 	"Zscore", "N", "P-value")

NAFLD <- NAFLD[, .(MarkerName, Allele1, Allele2, Effect, StdErr, Pvalue)]
NAFLD[,Zscore := Effect/StdErr]
NAFLD[,N := 8434 + 190120]
NAFLD <- NAFLD[,.(MarkerName, Allele1, Allele2, Zscore, N, Pvalue)]
colnames(NAFLD)<-colnoms


namefile <- "/home/gagelo01/workspace/Projects/Nooshin_NAFLD/Data/Modified/NAFLD_sumstat_LDHub"
fwrite(NAFLD, file = namefile,  col.names = TRUE, sep = "\t")

my_wd<-getwd() # save your current working directory path
dest_path<-"/home/gagelo01/workspace/Projects/Nooshin_NAFLD/Data/Modified/"
setwd(dest_path)
files<-"NAFLD_sumstat_LDHub"
named<-paste0(files,".zip")
mapply(zip,zipfile=named,files=files)
setwd(my_wd) # reset working directory path


zip(zipfile=paste0(namefile,".zip"), files=namefile)

  #Then go to LDHUB and transfer the file