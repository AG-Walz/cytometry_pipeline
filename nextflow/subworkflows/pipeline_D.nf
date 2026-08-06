/*
 * Pipeline D: FlowJo Workspace -> Gating Comparison
 * Flow: Load WSP -> Compare gatings -> Overlap statistics
 */

include { LOAD_FLOWJO_WSP }  from '../modules/loading/load_flowjo_wsp'
include { COMPARE_GATINGS }  from '../modules/metadata_qc/compare_gatings'
include { OVERLAP_STATS }    from '../modules/metadata_qc/overlap_statistics'

workflow PIPELINE_D {

    // Step 1: Load FlowJo workspaces — stage WSP files into work dir
    ch_wsp_files = Channel
        .of(params.wsp_paths.split(','))
        .map { it.trim() }
        .map { file(it) }
        .collect()

    ch_pop_strings = LOAD_FLOWJO_WSP(ch_wsp_files, params.wsp_names)

    // Step 2: Compare gating strategies
    ch_comparison = COMPARE_GATINGS(ch_pop_strings.pop_strings)

    // Step 3: Overlap statistics (Fleiss' kappa)
    OVERLAP_STATS(ch_comparison.label_lists)
}
