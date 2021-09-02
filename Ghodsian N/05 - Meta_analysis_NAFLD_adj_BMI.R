#!/usr/bin/env Rscript

library(data.table)

wd_data = "/home/couchr02/workspace/meta"

setwd(wd_data)

trait = "NAFLD"

gwas_ukb = "/home/couchr02/Mendel_Commun/Christian/GWAS/UKB/GWAS_UKB3_NAFLD_adj_BMI_SAIGE.txt.gz" 

gwas_emerge = "/home/couchr02/Mendel_Commun/Christian/GWAS/40_diseases/NamjouB_31311600_NAFLD_mod2.txt.gz"

gwas_EB = "/home/couchr02/Mendel_Commun/Christian/GWAS/Estonian_Biobank/NAFLD_EstBB_2817Cases_GWAS.txt.gz"

mytemp = "/home/couchr02/workspace/mytemp"

start_time = Sys.time()

### Pour GWAS de UKB

ukb = as.data.frame(fread(gwas_ukb,tmpdir=mytemp))
ukb_clean = ukb[grep("rs",ukb$rsid),]
ukb_clean = ukb_clean[which(ukb_clean$AF_Allele2>=0.001 & ukb_clean$AF_Allele2<=0.999),]

for (i in 1:22) {

  ukb_chr = ukb_clean[which(ukb_clean$CHR==i),]
  fwrite(as.data.frame(ukb_chr),paste0("UKB/ukb_",trait,"_",i,".txt"),sep="\t",row.names=F)

}

rm(ukb)
rm(ukb_clean)


### Pour GWAS de eMerge

emerge = as.data.frame(fread(gwas_emerge,tmpdir=mytemp))

for (i in 1:22) {

  emerge_chr = emerge[which(emerge$CHR==i),]
  fwrite(as.data.frame(emerge_chr),paste0("eMerge/emerge_",trait,"_",i,".txt"),sep="\t",row.names=F)

}

rm(emerge)


### Pour GWAS de Estonian Biobank

EB = as.data.frame(fread(gwas_EB,tmpdir=mytemp))
EB_clean = EB[grep("rs",EB$SNPID),]
EB_clean = EB_clean[which(EB_clean$AF_Allele2>=0.001 & EB_clean$AF_Allele2<=0.999),]

for (i in 1:22) {

  EB_chr = EB_clean[which(EB_clean$CHR==i),]
  fwrite(as.data.frame(EB_chr),paste0("EB/EB_",trait,"_",i,".txt"),sep="\t",row.names=F)

}

rm(EB)
rm(EB_clean)

### Création des scripts METAL et méta-analyse

for (i in 1:22) {

script = paste0("scripts/METAL_",trait,"_",i,".txt")
logfile = paste0("logfiles/METAL_",trait,"_",i,".log")

sink(script)

cat("SCHEME STDERR\n\n")

cat("MARKER   rsid\n")
cat("WEIGHT   donotusecolumn\n")
cat("DEFAULTWEIGHT 397799\n")
cat("ALLELE   Allele2 Allele1\n")
cat("EFFECT   BETA\n")
cat("STDERR   SE\n")
cat("PVAL     p.value\n")
cat("SEPARATOR TAB\n\n")

cat(paste0("PROCESS UKB/ukb_",trait,"_",i,".txt\n\n"))

cat("MARKER   SNP\n")
cat("ALLELE   A1 A2\n")
cat("WEIGHT   N\n")
cat("EFFECT   Beta\n")
cat("STDERR   SE\n")
cat("PVAL     P\n")
cat("SEPARATOR TAB\n\n")

cat(paste0("PROCESS eMerge/emerge_",trait,"_",i,".txt\n\n"))

cat("MARKER   SNPID\n")
cat("WEIGHT   N\n")
cat("ALLELE   Allele2 Allele1\n")
cat("EFFECT   BETA\n")
cat("STDERR   SE\n")
cat("PVAL     p.value\n")
cat("SEPARATOR TAB\n\n")

cat(paste0("PROCESS EB/EB_",trait,"_",i,".txt\n\n"))

cat(paste0("OUTFILE New_Results/",trait,"_",i,"_stderr_het_ _adj_BMI.tbl\n\n"))

cat("ANALYZE RANDOM\n")

cat("QUIT\n")

sink()

word1 = "./metal"
args1 = c(script,">",logfile)
system2(word1,args1)

### Ajout du chromosome et de la position

ukb = as.data.frame(fread(paste0("UKB/ukb_",trait,"_",i,".txt")))
snps = ukb[,c("rsid","CHR","POS")]
snps = snps[!duplicated(snps$rsid),]

meta = as.data.frame(fread(paste0("New_Results/",trait,"_",i,"_stderr_het_1_adj_BMI.tbl")))
meta = meta[grep("\\?",meta$Direction,invert=T),]
meta_new = merge(meta,snps,by.x="MarkerName",by.y="rsid")
fwrite(as.data.frame(meta_new),paste0("New_Results/",trait,"_",i,"_adj_BMI.txt"),row.names=F)

} # fin de la boucle i

# Concaténation des chromosomes

word2 = "head"
args2 = c("-1",paste0("New_Results/",trait,"_1_adj_BMI.txt > New_Results/meta_clean_",trait,"_adj_BMI.txt"))
word3 = "tail"
args3 = c("-n +2 -q",paste0("New_Results/",trait,"*_adj_BMI.txt >> New_Results/meta_clean_",trait,"_adj_BMI.txt"))
system2(word2,args2)
system2(word3,args3)

# Effacement des fichiers temporaires

word5 = "rm"
args5 = c(paste0("New_Results/",trait,"*"))
args6 = c("scripts/*")
args7 = c("UKB/*")
args9 = c("eMerge/*")
args10 = c("EB/*")
args11 = c("logfiles/*")
system2(word5,args5)
system2(word5,args6)
system2(word5,args7)
system2(word5,args8)
system2(word5,args9)
system2(word5,args10)
system2(word5,args11)

stop_time = Sys.time()

duree = stop_time - start_time
print(duree)

