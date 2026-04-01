/*
 * Pipeline A: Raw FCS -> Full Analysis
 * Flow: Load -> Preprocess -> Homogenize -> [Batch Correct] -> Downsample
 *       -> Cluster -> Analyze -> Visualize
 */

include { LOAD_FCS_RAW }       from '../modules/loading/load_fcs_raw'
include { PREPROCESS_SAMPLE }  from '../modules/preprocessing/preprocess_sample'
include { HOMOGENIZE_PANELS }  from '../modules/harmonization/homogenize'
include { BATCH_CORRECT }      from '../modules/harmonization/batch_correction'
include { DOWNSAMPLE }         from '../modules/preprocessing/downsample'
include { BUILD_EXPR_MATRIX }  from '../modules/clustering/build_expr_matrix'
include { CLUSTER }            from '../modules/clustering/cluster'
include { CLUSTERS_TO_CSV }    from '../modules/clustering/clusters_to_csv'
include { SPLIT_CLUSTERS }     from '../modules/analysis/cluster_splitting'
include { DIFF_ABUNDANCE }     from '../modules/analysis/differential_abundance'
include { DIMRED }             from '../modules/visualization/dimred'

workflow PIPELINE_A {

    // Step 1: Scatter — load each FCS file
    Channel
        .fromPath("${params.input_dir}/*.fcs")
        .map { fcs -> tuple(fcs.baseName.replaceAll(/\.fcs$/, ''), fcs) }
        .set { ch_fcs }

    ch_loaded = LOAD_FCS_RAW(ch_fcs)

    // Step 2: Scatter — preprocess each sample independently
    ch_preprocessed = PREPROCESS_SAMPLE(ch_loaded.flowframe)

    // Step 3: Gather — homogenize panels
    ch_preprocessed.flowframe
        .map { name, rds -> rds }
        .collect()
        .set { ch_all_pp }

    ch_homogenized = HOMOGENIZE_PANELS(ch_all_pp)

    // Step 3b: Optional batch correction (conditional on params.batch_method)
    if (params.batch_method != 'none') {
        ch_after_batch = BATCH_CORRECT(ch_homogenized.flowframes.collect())
        ch_for_ds = ch_after_batch.flowframes
            .flatten()
            .map { rds -> tuple(rds.baseName.replaceAll(/\.rds$/, ''), rds) }
    } else {
        ch_for_ds = ch_homogenized.flowframes
            .flatten()
            .map { rds -> tuple(rds.baseName.replaceAll(/\.rds$/, ''), rds) }
    }

    // Step 4: Re-scatter — downsample each sample
    ch_downsampled = DOWNSAMPLE(ch_for_ds)

    // Step 5: Gather — build expression matrix
    ch_downsampled.flowframe
        .map { name, rds -> rds }
        .collect()
        .set { ch_all_ds }

    ch_expr = BUILD_EXPR_MATRIX(ch_all_ds)

    // Step 6: Cluster
    ch_clusters = CLUSTER(ch_expr.expr_data)

    // Step 7: Generate gating CSVs and split
    ch_csvs = CLUSTERS_TO_CSV(ch_expr.expr_data, ch_clusters.cluster_ids)
    SPLIT_CLUSTERS(ch_csvs.csv_dir)

    // Step 8: Differential abundance (if conditions file provided)
    if (params.conditions_csv) {
        ch_conditions = Channel.fromPath(params.conditions_csv)
        DIFF_ABUNDANCE(ch_csvs.csv_dir, ch_conditions)
    }

    // Step 9: Visualization
    DIMRED(ch_expr.expr_data, ch_clusters.cluster_ids)
}
