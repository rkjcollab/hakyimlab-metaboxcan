#! /usr/bin/env nextflow

/*
==================================================
                    MetaboXcan
==================================================
 MetaboXcan Analysis Pipeline.
 #### Homepage / Documentation
 https:://github.com/hakyimlab/metaboxcan
---------------------------------------------------
*/

// Enable DSL 2 syntax
nextflow.enable.dsl = 2


def helpMessage(){
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run main.nf

    Mandatory arguments:
      --gene_models_folder [path]
      --metabolite_models_folder [path]
      --gwas_file [path]
      --gwas_name [char]
      --covariance [file]

    Options:
      --gwas_N [num]
      --gwas_h2 [num]
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// catching both -name and --name if specified by user
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

//Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']              = custom_runName ?: workflow.runName
summary['GWAS name']             = params.gwas_name
summary['GWAS file']             = params.gwas_file
if (params.gwas_N) {
    summary['GWAS sample size']  = params.gwas_N
    summary['GWAS heritability'] = params.gwas_h2
}
summary['Gene models']           = params.gene_models_folder
summary['Metabolite models']     = params.metabolite_models_folder
summary['Max Resources']         = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output dir']            = params.outdir
summary['Launch dir']            = workflow.launchDir
summary['Working dir']           = workflow.workDir
summary['Script dir']            = workflow.projectDir
summary['User']                  = workflow.userName

//log it out
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

//set channels
gene_models = Channel.fromFilePairs(params.gene_models_folder)
metabolite_models = Channel.fromFilePairs(params.metabolite_models_folder)
gwas = Channel
            .fromPath(params.gwas_file)
            .map {gwas_file -> tuple(params.gwas_name,gwas_file)}
gene_annot = Channel.fromPath(params.gene_annot)
gene2metabo = Channel.fromPath(params.gene2metabo)
mqtl_db = Channel.fromPath(params.mqtl)
ld_blocks = Channel.fromPath(params.ld_blocks)
metabolite_map = Channel.fromPath(params.metabolite_metadata)
omics_network = Channel.fromPath(params.omics_network)
gene_info = Channel.fromPath(params.gene_info)


// Import modules
include {smultixcan;gwas_database;html_report } from './modules/metaxcan.nf'

// Import workflows
include { spredixcan; smetaboxcan } from './workflow/main_workflow.nf' addParams(outdir: "${params.outdir}")

// Create workflows
workflow OMICS_PIPELINE {
   take:
      gene_models
      metabolite_models
      gwas
      gene2metabo
      gene_annot
      mqtl_db
      metabolite_map
      ld_blocks
      omics_network
      gene_info

   main:
      // Run Summary prediXcan
      spredixcan(gene_models,gwas)

      // Run smultixcan
      smultixcan(spredixcan.out.db.collect(),spredixcan.out.outname.collect(),gwas)

      // Run Summary metaboXcan
      smetaboxcan(metabolite_models,gwas)

      // Convert sumstat to db
      gwas_database(gwas)

      // Generate a report with a locus zoom plot
      html_report(spredixcan.out.outname.collect(), smultixcan.out.multixcan, smetaboxcan.out.outname,
                  metabolite_models,gwas_database.out,gene2metabo,gene_annot,mqtl_db,metabolite_map,
                  ld_blocks,omics_network,gene_info)
}

workflow {
    main:
      //gwas.view()
      OMICS_PIPELINE(gene_models,metabolite_models,gwas,gene2metabo,gene_annot,mqtl_db,metabolite_map,
                    ld_blocks,omics_network,gene_info)
}
