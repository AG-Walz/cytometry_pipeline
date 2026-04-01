#!/usr/bin/env nextflow

/*
 * =============================================================================
 * Cytometry Analysis Pipeline — Nextflow DSL2
 * =============================================================================
 * Modular flow cytometry analysis framework supporting raw FCS files,
 * preprocessed FCS, CSV expression matrices, FlowJo workspace output,
 * and CSV gating data.
 *
 * Usage:
 *   nextflow run nextflow/main.nf --pipeline A --input_dir /data/fcs ...
 *   nextflow run nextflow/main.nf -params-file nextflow/conf/pipeline_A_params.yaml
 *
 * Pipelines:
 *   A — Raw FCS -> Full Analysis (load, QC, transform, cluster, analyze, visualize)
 *   B — CSV Expression -> Statistics (marker stats, heatmap)
 *   C — CSV Gating -> Cluster Analysis (enrichment or abundance)
 *   D — FlowJo Comparison (gating comparison, Fleiss' kappa)
 *   E — Preprocessed FCS -> UMAP/Clustering (skip QC, cluster, visualize)
 * =============================================================================
 */

nextflow.enable.dsl = 2

include { PIPELINE_A } from './subworkflows/pipeline_A'
include { PIPELINE_B } from './subworkflows/pipeline_B'
include { PIPELINE_C } from './subworkflows/pipeline_C'
include { PIPELINE_D } from './subworkflows/pipeline_D'
include { PIPELINE_E } from './subworkflows/pipeline_E'

// Validate required params
def validateParams() {
    if (!params.pipeline) {
        error "ERROR: --pipeline is required. Choose A, B, C, D, or E."
    }
    if (!params.input_dir && params.pipeline != 'D') {
        error "ERROR: --input_dir is required for pipeline ${params.pipeline}."
    }
}

workflow {
    validateParams()

    log.info """
    =============================================================
    Cytometry Pipeline ${params.pipeline}
    =============================================================
    input_dir          : ${params.input_dir ?: 'N/A'}
    output_dir         : ${params.output_dir}
    markers            : ${params.markers_of_interest ?: 'auto'}
    clustering_method  : ${params.clustering_method}
    batch_method       : ${params.batch_method}
    dimred_method      : ${params.dimred_method}
    =============================================================
    """.stripIndent()

    switch (params.pipeline) {
        case 'A':
            PIPELINE_A()
            break
        case 'B':
            PIPELINE_B()
            break
        case 'C':
            PIPELINE_C()
            break
        case 'D':
            PIPELINE_D()
            break
        case 'E':
            PIPELINE_E()
            break
        default:
            error "Unknown pipeline '${params.pipeline}'. Choose A, B, C, D, or E."
    }
}
