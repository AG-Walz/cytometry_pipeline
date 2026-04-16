/*
 * Pipeline C: CSV Gating -> Cluster Analysis Pipeline
 * Flow: Load gating CSV -> Split clusters -> Enrichment/Abundance -> Barplots
 */

include { LOAD_CSV_GATING }     from '../modules/loading/load_csv_gating'
include { SPLIT_CLUSTERS }      from '../modules/analysis/cluster_splitting'
include { CLUSTER_ENRICHMENT }  from '../modules/analysis/cluster_enrichment'
include { DIFF_ABUNDANCE }      from '../modules/analysis/differential_abundance'
include { BARPLOT }             from '../modules/visualization/barplot'
include { ABUNDANCE_PLOT }      from '../modules/analysis/abundance_plotting'

workflow PIPELINE_C {

    ch_input_dir = Channel.fromPath(params.input_dir, type: 'dir')

    // Step 1: Load CSV gating data
    ch_gating = LOAD_CSV_GATING(ch_input_dir)

    // Step 2: Split by cluster (use original input CSVs)
    ch_split = SPLIT_CLUSTERS(ch_input_dir)

    // Step 3: Conditional analysis
    if (params.conditions_csv) {
        ch_conditions = Channel.fromPath(params.conditions_csv)

        if (params.analysis_method == 'enrichment') {
            CLUSTER_ENRICHMENT(ch_input_dir, ch_conditions)
        } else {
            DIFF_ABUNDANCE(ch_input_dir, ch_conditions)
        }
    }

    // Step 4: Composition barplot (runs when split produces output)
    BARPLOT(ch_split.split_result)

    // Step 5: Abundance analysis
    if (params.conditions_csv) {
        ABUNDANCE_PLOT(ch_input_dir)
    }
}
