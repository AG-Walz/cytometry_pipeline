process DIFF_ABUNDANCE {
    tag "diff_abundance"
    label 'process_medium'

    publishDir "${params.output_dir}/analysis/differential_abundance", mode: 'copy'

    input:
    path csv_dir
    path conditions_csv

    output:
    path "diff_abundance/", emit: results

    script:
    """
    diff_abundance_cli.R \
        --data_dir ${csv_dir} \
        --conditions_csv ${conditions_csv} \
        --output_dir diff_abundance \
        --cluster_regex '${params.cluster_regex}' \
        --condition_levels '${params.condition_levels}' \
        --min_cluster_size ${params.min_cluster_size} \
        --p_adjust_method ${params.p_adjust_method} \
        --seed ${params.seed}
    """
}
