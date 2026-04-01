process MARKER_STATISTICS {
    tag "marker_stats"
    label 'process_medium'

    publishDir "${params.output_dir}/statistics", mode: 'copy'

    input:
    path csv_dir

    output:
    path "marker_stats/", emit: results

    script:
    def markers_arg = params.markers_of_interest ? "--markers ${params.markers_of_interest}" : ""
    """
    marker_stats_cli.R \
        --input_folder ${csv_dir} \
        --output_folder marker_stats \
        ${markers_arg} \
        --test_type ${params.test_type}
    """
}
