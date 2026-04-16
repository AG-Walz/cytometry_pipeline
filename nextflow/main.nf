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
 *   nextflow run nextflow/main.nf -params-file config/pipeline_A_params.yaml
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
    // Normalize case (accept both 'a' and 'A')
    params.pipeline = params.pipeline.toUpperCase()

    if (!params.input_dir && params.pipeline != 'D') {
        error "ERROR: --input_dir is required for pipeline ${params.pipeline}."
    }
    if (params.pipeline in ['A', 'B', 'E'] && !params.markers_of_interest) {
        error "ERROR: --markers_of_interest is required for pipeline ${params.pipeline}. Provide comma-separated marker names, e.g. --markers_of_interest 'CD3,CD4,CD8'"
    }
    if (params.pipeline == 'D') {
        if (!params.wsp_paths) {
            error "ERROR: --wsp_paths is required for pipeline D. Provide comma-separated paths to .wsp files."
        }
        if (!params.wsp_names) {
            error "ERROR: --wsp_names is required for pipeline D. Provide comma-separated workspace names."
        }
    }
    // Validate method choices
    def valid_clustering = ['leiden', 'flowsom', 'phenograph', 'dbscan']
    if (!(params.clustering_method in valid_clustering)) {
        error "ERROR: --clustering_method '${params.clustering_method}' is not valid. Choose from: ${valid_clustering.join(', ')}"
    }
    def valid_batch = ['none', 'harmony', 'cytonorm', 'combat']
    if (!(params.batch_method in valid_batch)) {
        error "ERROR: --batch_method '${params.batch_method}' is not valid. Choose from: ${valid_batch.join(', ')}"
    }
    def valid_dimred = ['umap', 'tsne', 'diffusion_map']
    if (!(params.dimred_method in valid_dimred)) {
        error "ERROR: --dimred_method '${params.dimred_method}' is not valid. Choose from: ${valid_dimred.join(', ')}"
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
