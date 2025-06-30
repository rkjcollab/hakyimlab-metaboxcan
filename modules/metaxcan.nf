#! /usr/bin/env nextflow

process impute_gene {
    tag "PrediXcan"
    label "process_medium"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/predicted_gene" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'move'

    input:
    path(model_db)
    path(geno)
    val(prefix)

    output:
    path(predicted)
    path(summary)

    script:
    predicted = prefix + '_predictions.txt'
    summary = prefix + '_predictions_summary'
    // use template to change parameters easily
    """
    python ${params.MetaXcan}/software/Predict.py --help > log.txt
      --model_db_path ${model_db} \
      --model_db_snp_key rsid \
      --vcf_genotypes ${geno} \
      --vcf_mode genotyped \
      --prediction_output ${predicted} \
      --prediction_summary_output ${summary} \
      --skip_palindromic \
      --verbosity 9 \
      --throw
    """
}


process spredixcan_raw {
  tag "individual tissue association raw"
  label "process_medium"
  publishDir path: { params.keepIntermediate ? "${params.outdir}/gene-trait" : params.outdir },
             saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

  input:
  tuple val(tissue), path(db_cov)
  tuple val(gwas_name), path(gwas_file)

  output:
  path(outname), emit: outname
  path(db), emit: db

  script:
  outname = gwas_name + "_GT_" + tissue + ".csv"
  db = db_cov[0]
  covs = db_cov[1]
  
  """
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
  """
}

process spredixcan_adj {
  tag "individual tissue association adj"
  label "process_medium"
  publishDir path: { params.keepIntermediate ? "${params.outdir}/gene-trait" : params.outdir },
             saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

  input:
  tuple val(tissue), path(db_cov)
  tuple val(gwas_name), path(gwas_file)
  val gwas_N
  val gwas_h2

  output:
  path(outname), emit: outname
  path(db), emit: db

  script:
  outname = gwas_name + "_GT_" + tissue + ".csv"
  db = db_cov[0]
  covs = db_cov[1]
  
  """
  python ${params.MetaXcan}/software/SPrediXcan.py \
  --gwas_file ${gwas_file} \
  --gwas_N ${gwas_N} \
  --gwas_h2 ${gwas_h2} \
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
  """
}

process smultixcan {
  tag "meta analyzed associations"
  label "process_medium"
  publishDir path: { params.keepIntermediate ? "${params.outdir}/mgene-trait" : params.outdir },
             saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

  input:
  path(gene_models)
  path(gene_traits)
  tuple val(gwas_name), path(gwas_file)
  output:
  path(mgene_trait), emit: multixcan
  script:
  mgene_trait = "multixcan_" + gwas_name + "_results.txt"
  metax_filt = gwas_name + "_(.*).csv"
  metax_parse = "(.*)_GT_(.*).csv"

  template 'smultixcan.sh'
}

process smetaboxcan_raw {
  tag "individual metabolite association raw"
  label "process_medium"
  publishDir path: { params.keepIntermediate ? "${params.outdir}/metabo-trait" : params.outdir },
             saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

  input:
  tuple val(tissue), path(db_cov)
  tuple val(gwas_name), path(gwas_file)

  output:
  path(outname), emit: outname

  script:
  outname = gwas_name + "_" + tissue + ".csv"
  db = db_cov[0]
  covs = db_cov[1]
  
  template 'smetaboxcan.sh'
}

process smetaboxcan_adj {
  tag "individual metabolite association adj"
  label "process_medium"
  publishDir path: { params.keepIntermediate ? "${params.outdir}/metabo-trait" : params.outdir },
             saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

  input:
  tuple val(tissue), path(db_cov)
  tuple val(gwas_name), path(gwas_file)
  val gwas_N
  val gwas_h2

  output:
  path(outname), emit: outname

  script:
  outname = gwas_name + "_" + tissue + ".csv"
  db = db_cov[0]
  covs = db_cov[1]
  
  template 'smetaboxcan_adj.sh'
}

process gwas_database {
  tag "covert gwas sumstat into a database"
  label "process_medium"
  input:
  tuple val(gwas_name), path(gwas_file)

  output:
  tuple val(gwas_name), path(out_db)

  script:
  out_db = gwas_name + ".db"

  """
  trait2db.R --input_file ${gwas_file} --output_db ${out_db}
  """

}


process html_report {
  tag "generate report"
  label "process_medium"
  publishDir path: "${params.outdir}/report",
             saveAs: it, mode: 'copy'

  input:
  path spred_results
  path multixcan
  path smetaboxcan
  tuple val(tissue), path(db_cov)
  tuple val(gwas_name), path(gwas_db)
  path gene2metabo
  path gene_annot
  path mqtl_db
  path metabo_map
  path ld_blocks
  path omics_network
  path gene_info

  output:
  path(report), emit: viz_report
  path "nets/*", emit: full_network
  path "summaries/*", emit: table_summary

  script:
  report = gwas_name + "_report.html"
  metabo_db = db_cov[0]
  
  """
  report_compilation_db.R \
    --metaxcan_dir ./ \
    --smultixcan_file ${multixcan} \
    --smetaboxcan_file ${smetaboxcan} \
    --predict_db ${metabo_db} \
    --phenotype_name ${gwas_name} \
    --gwas_file ${gwas_db} \
    --gene_metabo_assoc_file ${gene2metabo} \
    --gene_annotation ${gene_annot} \
    --mqtl_db ${mqtl_db} \
    --metabo_pathways ${metabo_map} \
    --ld_blocks ${ld_blocks} \
    --omics_network ${omics_network} \
    --gene_info ${gene_info} \
    --ntop 20 \
    --rmd_template ${baseDir}/bin/report_template_dbs.Rmd \
    --output_html ${report}
  """
}

