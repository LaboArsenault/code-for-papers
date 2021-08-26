#!/bin/bash

# Step 1

Rscript step1_fitNULLGLMM.R \
	--plinkFile=./input_E03/ukb_cal_allchr_pruned_93ksnps \
	--phenoFile=./input_NAFLD/pheno_NAFLD_ukb3.txt \
	--phenoCol=NAFLD \
	--covarColList=Sex,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
	--sampleIDColinphenoFile=ID \
	--traitType=binary \
	--outputPrefix=./output_NAFLD/UKB_93ksnps_408kSuj_binary_NAFLD_step1 \
	--nThreads=5 \
	--LOCO=FALSE

# Step 2

for i in {1..22}
do
time Rscript step2_SPAtests.R \
       --bgenFile=/home/couchr02/Mendel_UKB/Source/Imputed_data_V3/BGEN/ukb_imp_chr${i}_v3.bgen \
	   --bgenFileIndex=/home/couchr02/Mendel_UKB/Source/Imputed_data_V3/BGEN/ukb_imp_chr${i}_v3.bgen.bgi \
	   --chrom=${i} \
	   --minMAF=0.001 \
	   --minMAC=3 \
	   --sampleFile=./input_E03/ukb_v3_samples.txt \
	   --GMMATmodelFile=./output_NAFLD/UKB_93ksnps_408kSuj_binary_NAFLD_step1.rda \
	   --varianceRatioFile=./output_NAFLD/UKB_93ksnps_408kSuj_binary_NAFLD_step1.varianceRatio.txt \
	   --SAIGEOutputFile=./output_NAFLD/UKB_V3_400kSuj_NAFLD.chr${i}.SAIGE.dropmissing.txt \
	   --numLinesOutput=2 \
	   --IsOutputAFinCaseCtrl=TRUE \
	   --IsDropMissingDosages=TRUE \
	   > ./output_NAFLD/UKB_V3_NAFLD_chr${i}.log &
done
