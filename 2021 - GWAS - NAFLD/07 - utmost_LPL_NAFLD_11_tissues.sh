#!/bin/bash

TISSUE_GTEx=(Adipose_Subcutaneous Artery_Aorta Brain_Cerebellum Heart_Left_Ventricle Lung Nerve_Tibial \
	Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Spleen Thyroid Whole_Blood)
for TISSUE in ${TISSUE_GTEx[@]}
do
python ./single_tissue_association_test.py \
  --model_db_path /home/couchr02/Mendel_Commun/Christian/GTEx8_EUR_models/GTEx8_eur_dbs/GTEx8_${TISSUE}_HapMap_alpha0.5_window1e6_filtered.db \
  --covariance /home/couchr02/Mendel_Commun/Christian/GTEx8_EUR_models/GTEx8_eur_cov/GTEx8_${TISSUE}_HapMap_alpha0.5_window1e6_cut.txt.gz \
  --gwas_folder . \
  --gwas_file_pattern meta_NAFLD_LPL_1000KB.txt \
  --snp_column MarkerName \
  --effect_allele_column Allele1 \
  --non_effect_allele_column Allele2 \
  --beta_column Effect \
  --pvalue_column Pvalue \
  --output_file ./results/LPL_${TISSUE}_NAFLD_1000KB.csv
done

