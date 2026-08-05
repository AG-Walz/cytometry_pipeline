process EXTRACT_FLOWJO_POPULATION {
    tag "flowjo_${params.population}"
    label 'process_high'

    publishDir "${params.output_dir}/flowjo_population", mode: 'copy'

    input:
    val wsp
    val fcs_dir

    output:
    path "flowjo_population/expr_data.rds",     emit: expr_data
    path "flowjo_population/subpop_ids.rds",    emit: subpop_ids
    path "flowjo_population/flowjo_labels.rds", emit: labels
    path "flowjo_population/load_report.json",  emit: report

    script:
    // Optional arguments are built up rather than always passed: an empty
    // --marker_aliases or --population_map is not the same as an absent one on
    // the R side, and passing "null" would be read as a value.
    def alias_arg   = params.marker_aliases   ? "--marker_aliases '${params.marker_aliases}'"   : ''
    def popmap_arg  = params.population_map   ? "--population_map '${params.population_map}'"   : ''
    def markers_arg = params.markers          ? "--markers '${params.markers}'"                 : ''
    def subsel_arg  = params.subpop_select    ? "--subpop_select '${params.subpop_select}'"     : ''
    def meta_arg    = params.sample_metadata  ? "--sample_metadata '${params.sample_metadata}'" : ''
    def pattern_arg = params.sample_pattern   ? "--sample_pattern '${params.sample_pattern}'"   : ''
    def annot_arg   = params.annotate         ? "--annotate '${params.annotate}'"               : ''
    def allow_arg   = params.allow_missing_samples ? '--allow_missing_samples' : ''
    def norecover   = params.recover_markers ? '' : '--no_recover_markers'
    """
    extract_flowjo_population_cli.R \
        --wsp '${wsp}' \
        --fcs_dir '${fcs_dir}' \
        --population '${params.population}' \
        --group '${params.group}' \
        --match_mode ${params.match_mode} \
        --subpop_mode ${params.subpop_mode} \
        --normalize ${params.normalize} \
        --downsample ${params.downsample} \
        --seed ${params.seed} \
        --output_dir flowjo_population \
        ${allow_arg} ${norecover} ${alias_arg} ${popmap_arg} ${markers_arg} \
        ${subsel_arg} ${meta_arg} ${pattern_arg} ${annot_arg}
    """
}
