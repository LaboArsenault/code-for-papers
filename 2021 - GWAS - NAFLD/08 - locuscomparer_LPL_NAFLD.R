#!/usr/bin/env Rscript

wd_data = "/home/couchr02/workspace/coloc/"

setwd(wd_data)

library(data.table)
library(locuscomparer)
library(ggplot2)

gene = "LPL"
tissue = "Adipose_Subcutaneous"

nafld = as.data.frame(fread("/home/couchr02/Mendel_Commun/Christian/GWAS/40_diseases/meta_clean_NAFLD.txt"))
snps = as.data.frame(fread("/home/couchr02/workspace/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt"))
gencode = as.data.frame(fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt"))

gencode_small = gencode[which(gencode$gene_name==gene),]
gene_id = unique(gencode_small$gene_id)
chr = unique(gencode_small$chr)
mystart = min(gencode_small$start)
myend = max(gencode_small$end)

QTL = as.data.frame(fread(paste0("/home/couchr02/workspace/GTEx_v8/QTLtools/All_tissues_",gene,"_QTLtools_nominal_1_window_1000kb.txt")))

QTL = QTL[which(QTL$tissue==tissue),]

merged = merge(QTL,snps,by.x=c("chr_top_variant","start_top_variant"),by.y=c("chr","pos_b38"),sort=F)
merged = merged[which(merged$rsid %in% nafld$MarkerName),]
merged = merged[,c("rsid","p_val")]
colnames(merged) = c("rsid","pval")

nafld = nafld[which(nafld$MarkerName %in% merged$rsid),]
nafld = nafld[,c("MarkerName","Pvalue")]
colnames(nafld) = c("rsid","pval")

locus_zoom <- locuscompare(
  in_fn1 = merged,
  in_fn2 = nafld,
  snp = "rs13702"
)

ggsave(locus_zoom,file=paste0("Results/",gene,"_",tissue,"_NAFLD.png"),width=20, height=15,units="cm",dpi=200,type="cairo")
