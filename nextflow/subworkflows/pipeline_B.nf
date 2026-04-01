/*
 * Pipeline B: CSV Expression -> Statistics Pipeline
 * Flow: Load CSV expression -> Marker statistics -> Heatmap
 */

include { LOAD_CSV_EXPRESSION } from '../modules/loading/load_csv_expression'
include { MARKER_STATISTICS }   from '../modules/analysis/marker_statistics'
include { HEATMAP }             from '../modules/visualization/heatmap'

workflow PIPELINE_B {

    // Step 1: Load CSV expression data
    ch_csv_dir = Channel.fromPath(params.input_dir, type: 'dir')

    ch_loaded = LOAD_CSV_EXPRESSION(ch_csv_dir)

    // Step 2: Marker statistics
    MARKER_STATISTICS(ch_csv_dir)

    // Step 3: Heatmap visualization
    HEATMAP(ch_loaded.expr_list)
}
