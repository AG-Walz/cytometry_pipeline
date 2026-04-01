process ABUNDANCE_PLOT {
    tag "abundance_plot"
    label 'process_low'

    publishDir "${params.output_dir}/analysis/abundance", mode: 'copy'

    input:
    path csv_dir

    output:
    path "abundance/", emit: plots

    script:
    def metadata_arg = params.cytolution_metadata_file ? "--metadata_file ${params.cytolution_metadata_file}" : ""
    def exclude_arg = params.exclude_patterns ? "--exclude_patterns '${params.exclude_patterns}'" : ""
    """
    abundance_plot_cli.R \
        --input_dir ${csv_dir} \
        --output_dir abundance \
        --cluster_regex '${params.cluster_regex}' \
        ${exclude_arg} \
        ${metadata_arg}
    """
}
