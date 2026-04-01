/*
 * Pipeline E: Preprocessed FCS -> UMAP/Clustering
 * Flow: Load FCS (no QC) -> Homogenize -> Downsample -> Cluster -> UMAP
 */

include { LOAD_FCS_RAW }       from '../modules/loading/load_fcs_raw'
include { HOMOGENIZE_PANELS }  from '../modules/harmonization/homogenize'
include { DOWNSAMPLE }         from '../modules/preprocessing/downsample'
include { BUILD_EXPR_MATRIX }  from '../modules/clustering/build_expr_matrix'
include { CLUSTER }            from '../modules/clustering/cluster'
include { DIMRED }             from '../modules/visualization/dimred'

workflow PIPELINE_E {

    // Step 1: Scatter — load each FCS file independently
    Channel
        .fromPath("${params.input_dir}/*.fcs")
        .map { fcs -> tuple(fcs.baseName.replaceAll(/\.fcs$/, ''), fcs) }
        .set { ch_fcs }

    ch_loaded = LOAD_FCS_RAW(ch_fcs)

    // Step 2: Gather — homogenize panels (needs all samples to find common markers)
    ch_loaded.flowframe
        .map { name, rds -> rds }
        .collect()
        .set { ch_all_rds }

    ch_homogenized = HOMOGENIZE_PANELS(ch_all_rds)

    // Step 3: Re-scatter — downsample each sample independently
    ch_homogenized.flowframes
        .flatten()
        .map { rds -> tuple(rds.baseName.replaceAll(/\.rds$/, ''), rds) }
        .set { ch_per_sample }

    ch_downsampled = DOWNSAMPLE(ch_per_sample)

    // Step 4: Gather — build combined expression matrix
    ch_downsampled.flowframe
        .map { name, rds -> rds }
        .collect()
        .set { ch_all_ds }

    ch_expr = BUILD_EXPR_MATRIX(ch_all_ds)

    // Step 5: Cluster (whole dataset)
    ch_clusters = CLUSTER(ch_expr.expr_data)

    // Step 6: Dimensionality reduction visualization
    DIMRED(ch_expr.expr_data, ch_clusters.cluster_ids)
}
