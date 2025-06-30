#! /usr/bin/env nextflow

include { spredixcan_raw;spredixcan_adj;smetaboxcan_raw; smetaboxcan_adj} from '../modules/metaxcan.nf'

workflow spredixcan {
    take:
        gene_models
        gwas
    main:
    if (params.gwas_N && params.gwas_h2) {
        // Adjust for inflation
        spred = spredixcan_adj(gene_models,gwas.first(),params.gwas_N,params.gwas_h2)
    } else {
        // Run the old version
        spred = spredixcan_raw(gene_models,gwas.first())
    }

    emit:
        outname = spred.outname
        db = spred.db
}

workflow smetaboxcan {
    take:
        metabolite_models
        gwas
    main:
    if (params.gwas_N && params.gwas_h2) {
        // Adjust for inflation
        spred = smetaboxcan_adj(metabolite_models,gwas,params.gwas_N,params.gwas_h2)
    } else {
        // Run the old version
        spred = smetaboxcan_raw(metabolite_models,gwas)
    }

    emit:
        outname = spred.outname
}