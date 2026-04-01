process CLUSTERS_TO_CSV {
    tag "clusters_to_csv"
    label 'process_medium'

    publishDir "${params.output_dir}/analysis/cluster_data", mode: 'copy'

    input:
    path expr_data_rds
    path cluster_ids_rds

    output:
    path "gating_csvs/", emit: csv_dir

    script:
    """
    clusters_to_csv_cli.R \
        --expr_rds ${expr_data_rds} \
        --cluster_rds ${cluster_ids_rds} \
        --output_dir gating_csvs
    """
}
