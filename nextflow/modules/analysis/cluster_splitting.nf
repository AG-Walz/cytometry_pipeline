process SPLIT_CLUSTERS {
    tag "split_clusters"
    label 'process_medium'

    publishDir "${params.output_dir}/analysis/split_clusters", mode: 'copy'

    input:
    path csv_dir

    output:
    path "split_clusters/",          emit: split_dir
    path "split_clusters/split_result.rds", emit: split_result, optional: true

    script:
    def markers_arg = params.markers_of_interest ? "--markers_of_interest ${params.markers_of_interest}" : ""
    """
    split_clusters_cli.R \
        --input_dir ${csv_dir} \
        --output_dir split_clusters \
        --cluster_regex '${params.cluster_regex}' \
        ${markers_arg}
    """
}
