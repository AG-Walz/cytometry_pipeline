/*
 * Pipeline F: FlowJo workspaces -> one shared UMAP of a gated population
 *
 * Flow: Extract FlowJo population -> Dimred -> Facet plots
 *
 * The same shape as pipeline E, entered one step later. E starts from raw FCS
 * and has to homogenize, downsample and build a matrix before it can embed
 * anything; here the workspaces already say which cells are wanted, so
 * EXTRACT_FLOWJO_POPULATION does the gating and hands on the SAME
 * list(combined_matrix, sample_labels) that BUILD_EXPR_MATRIX produces. Every
 * step after that is a step the pipeline already had.
 *
 * Panel homogenization and downsampling happen inside the extraction rather
 * than as separate processes, because both have to be decided per SAMPLE while
 * the GatingSets are open — the marker names are recovered from the channel
 * metadata of the very files being gated, and downsampling has to happen after
 * gating to be meaningful. Doing them as later processes would mean writing
 * every cell of every sample to disk first.
 *
 * The FlowJo subgate of each cell is passed to DIMRED as its cluster labels.
 * A gate IS a cluster as far as the plotting is concerned, so the standard
 * by-cluster figure comes out correct with no extra code, and UMAP_FACETS adds
 * the colourings that only a gating tree can produce.
 */

include { EXTRACT_FLOWJO_POPULATION } from '../modules/loading/extract_flowjo_population'
include { DIMRED }                    from '../modules/visualization/dimred'
include { UMAP_FACETS }               from '../modules/visualization/umap_facets'

workflow PIPELINE_F {

    // Step 1: workspaces + FCS -> gated expression matrix, subgate per cell,
    // and the load report the GUI reads.
    ch_pop = EXTRACT_FLOWJO_POPULATION(
        Channel.value(params.wsp),
        Channel.value(params.fcs_dir)
    )

    // Step 2: embed once. --output_rds is what lets the facet step re-colour
    // this embedding instead of computing a second, different one.
    ch_dimred = DIMRED(ch_pop.expr_data, ch_pop.subpop_ids)

    // Step 3: the FlowJo-specific figures.
    UMAP_FACETS(
        ch_pop.expr_data,
        ch_pop.labels,
        ch_pop.subpop_ids,
        ch_dimred.coords
    )
}
