#!/usr/bin/env Rscript

library(data.table)

wd_data = "/home/couchr02/workspace/meta"

setwd(wd_data)

trait = "NAFLD"

gwas_ukb = "/home/couchr02/Mendel_Commun/Christian/GWAS/UKB/GWAS_UKB3_NAFLD_SAIGE.txt.gz" 

gwas_finn = "/home/couchr02/Mendel_Commun/FinnGen_r4/GWAS/finngen_R4_NAFLD.gz"

gwas_emerge = "/home/couchr02/Mendel_Commun/Christian/GWAS/40_diseases/NamjouB_31311600_NAFLD_mod2.txt.gz"

gwas_EB = "/home/couchr02/Mendel_Commun/Christian/GWAS/Estonian_Biobank/NAFLD_EstBB_4119Cases_GWAS.txt.gz"

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


### Pour GWAS de FinnGen

finn = as.data.frame(fread(gwas_finn,tmpdir=mytemp))
finn_clean = finn[grep("rs",finn$rsids),]
finn_clean = finn_clean[which(finn_clean$maf>=0.001 & finn_clean$maf<=0.999),]

for (i in 1:22) {

  finn_chr = finn_clean[which(finn_clean$"#chrom"==i),]
  fwrite(as.data.frame(finn_chr),paste0("FinnGen/finn_",trait,"_",i,".txt"),sep="\t",row.names=F)

}

rm(finn)
rm(finn_clean)

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

#cat("GENOMICCONTROL ON\n")

cat("MARKER   rsid\n")
cat("WEIGHT   donotusecolumn\n")
cat("DEFAULTWEIGHT 397799\n")
cat("ALLELE   Allele2 Allele1\n")
cat("EFFECT   BETA\n")
cat("STDERR   SE\n")
cat("PVAL     p.value\n")
cat("SEPARATOR TAB\n\n")

cat(paste0("PROCESS UKB/ukb_",trait,"_",i,".txt\n\n"))

cat("MARKER   rsids\n")
cat("ALLELE   alt ref\n")
cat("WEIGHT   donotusecolumn\n")
cat("DEFAULTWEIGHT 176900\n")
cat("EFFECT   beta\n")
cat("STDERR   sebeta\n")
cat("PVAL     pval\n")
cat("SEPARATOR TAB\n\n")

cat(paste0("PROCESS FinnGen/finn_",trait,"_",i,".txt\n\n"))

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

cat(paste0("OUTFILE New_Results/",trait,"_",i,"_stderr_het_ _new.tbl\n\n"))

cat("ANALYZE RANDOM\n")

cat("QUIT\n")

sink()

word1 = "./metal"
args1 = c(script,">",logfile)
system2(word1,args1)

### Ajout du chromosome et de la position et nettoyage du fichier de résultats (nom "rs" seulement)

ukb = as.data.frame(fread(paste0("UKB/ukb_",trait,"_",i,".txt")))
ukb = ukb[,c(3,1,2)]
colnames(ukb) = c("SNP","CHR","POS")
 
finn=as.data.frame(fread(paste0("FinnGen/finn_",trait,"_",i,".txt")))
finn=finn[,c(5,1,2)]
colnames(finn) = c("SNP","CHR","POS") # attention, la position est sur le build38 ici!!!

emerge=as.data.frame(fread(paste0("eMerge/emerge_",trait,"_",i,".txt")))
emerge=emerge[,c(1,4,5)]
colnames(emerge) = c("SNP","CHR","POS")

EB=as.data.frame(fread(paste0("EB/EB_",trait,"_",i,".txt")))
EB=EB[,c(3,1,2)]
colnames(EB) = c("SNP","CHR","POS")

snps = rbind(ukb,finn,emerge,EB)
#snps = rbind(ukb,finn)
snps = snps[!duplicated(snps$SNP),]

meta=as.data.frame(fread(paste0("New_Results/",trait,"_",i,"_stderr_het_1_new.tbl")))
meta_new=merge(meta,snps,by.x=1,by.y=1)
fwrite(as.data.frame(meta_new),paste0("New_Results/",trait,"_",i,"_new.txt"),row.names=F)

} # fin de la boucle i

# Concaténation des chromosomes
# La méta-analyse "brut" risque d'avoir des mélanges de positions entre le build37 et le build38, Ne pas l'utiliser (à moins
# de corriger les positions).

word2 = "head"
args2 = c("-1",paste0("New_Results/",trait,"_1_new.txt > New_Results/meta_",trait,".txt"))
word3 = "tail"
args3 = c("-n +2 -q",paste0("New_Results/",trait,"*_new.txt >> New_Results/meta_",trait,".txt"))
system2(word2,args2)
system2(word3,args3)

# Création du fichier contenant seulement les SNPs communs à toutes les GWAS...
# Ici, on est correct. Le "duplicated" de la ligne 173 va garder la première instance des SNPs communs aux 4 GWAS. Donc, la
# position gardée sera celle de UKB, qui est sur le build37. On n'aura pas de mélanges de positions entre le build37 et le build38.

meta_comm = as.data.frame(fread(paste0("New_Results/meta_",trait,".txt")))
meta_comm = meta_comm[grep("\\?",meta_comm$Direction,invert=T),]
fwrite(as.data.frame(meta_comm),paste0("New_Results/meta_clean_",trait,".txt"),row.names=F)

# Compression du fichier brut

word4 = "gzip"
args4 = c(paste0("New_Results/meta_",trait,".txt"))
system2(word4,args4)

# Effacement des fichiers temporaires

word5 = "rm"
args5 = c(paste0("New_Results/",trait,"*"))
args6 = c("scripts/*")
args7 = c("UKB/*")
args8 = c("FinnGen/*")
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

