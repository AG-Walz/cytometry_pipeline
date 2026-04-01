process CLUSTER_ENRICHMENT {
    tag "enrichment"
    label 'process_medium'

    publishDir "${params.output_dir}/analysis/enrichment", mode: 'copy'

    input:
    path csv_dir
    path conditions_csv

    output:
    path "enrichment/", emit: results

    script:
    """
    cluster_enrichment_cli.R \
        --data_dir ${csv_dir} \
        --conditions_csv ${conditions_csv} \
        --output_dir enrichment \
        --cluster_regex '${params.cluster_regex}' \
        --condition_levels '${params.condition_levels}' \
        --stat_method ${params.stat_method} \
        --p_adjust_method ${params.p_adjust_method} \
        --seed ${params.seed}
    """
}
