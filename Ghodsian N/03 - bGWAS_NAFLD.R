#!/usr/bin/env Rscript

### Afin de ne pas avoir à spécifier l'emplacement de ZMatrices à chaque fois, j'ai créé un lien symbolique dans mon répertoire home
### pointant vers le répertoire en question (~/workspace/bGWAS/ZMatrices)

library(bGWAS)
library(data.table)
library(ggplot2)  # afin de sauvegarder les graphiques avec ggsave

options(tibble.width = Inf)  # ne tronque pas les valeurs et affiche toutes les colonnes

wd_data = "/home/couchr02/workspace/bGWAS/"

setwd(wd_data)

nafld = as.data.frame(fread("/mnt/sde/couchr02/meta/New_Results/meta_clean_NAFLD.txt"))
nafld = nafld [,c("MarkerName","Allele1","Allele2","Effect","StdErr")]
colnames(nafld) = c("snp","a1","a2","beta","se")
nafld$a1 = toupper(nafld$a1)  # important d'avoir des majuscules sinon les SNPs ne seront pas reconnus
nafld$a2 = toupper(nafld$a2)  # important d'avoir des majuscules sinon les SNPs ne seront pas reconnus

# AllStudies = list_priorGWASs()  # pas nécessaire si on sait déjà quelles "prior" GWAS utiliser

### On veut ces "prior" GWAS pour ces analyses:

### BMI (GIANT) => ID=1
### Triglycerides (GLGC) => ID=16

MyStudies = c(1,16)

name_bGWAS = "NAFLD"

A = bGWAS(name = name_bGWAS,
          GWAS = nafld,
		  prior_studies = MyStudies,
		  stepwise_threshold = 0.05,
		  res_pruning_LD = 0.1,  # important, sinon le pruning se fait par défaut dans une fenêtre de 1000Kb (±500Kb)
         )

fwrite(as.data.frame(A$all_BFs),paste0(wd_data,name_bGWAS,"/NAFLD_all_BFs.txt"),row.names=F)		 
		 
BF_hits = extract_results_bGWAS(A, SNPs = "significant")
write.csv2(BF_hits,paste0(wd_data,name_bGWAS,"/NAFLD_BF_hits.csv"),row.names=F)

post_hits = extract_results_bGWAS(A, SNPs = "significant", results = "posterior")
write.csv2(post_hits,paste0(wd_data,name_bGWAS,"/NAFLD_post_hits.csv"),row.names=F) 

dir_hits = extract_results_bGWAS(A, SNPs = "significant", results = "direct")
write.csv2(dir_hits,paste0(wd_data,name_bGWAS,"/NAFLD_dir_hits.csv"),row.names=F)

ggsave(coefficients_plot_bGWAS(A),file="NAFLD/NAFLD_coeff_plot.png",type="cairo")
ggsave(manhattan_plot_bGWAS(A),file="NAFLD/NAFLD_manh_BF_plot.png",width=25,height=15,type="cairo")
ggsave(manhattan_plot_bGWAS(A, results="posterior"),file="NAFLD/NAFLD_manh_post_plot.png",width=25,height=15,type="cairo")
ggsave(manhattan_plot_bGWAS(A, results="direct"),file="NAFLD/NAFLD_manh_direct_plot.png",width=25,height=15,type="cairo")
