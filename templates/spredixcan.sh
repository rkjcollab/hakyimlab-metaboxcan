
python ${params.MetaXcan}/software/SPrediXcan.py \
  --gwas_file ${gwas_file} \
  --snp_column variant_id \
  --effect_allele_column effect_allele \
  --non_effect_allele_column non_effect_allele \
  --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${covs} \
  --keep_non_rsid \
  --additional_output \
  --model_db_snp_key rsid \
  --throw \
  --output_file ${outname}
