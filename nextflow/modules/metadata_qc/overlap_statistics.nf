process OVERLAP_STATS {
    tag "overlap_stats"
    label 'process_low'

    publishDir "${params.output_dir}/overlap_stats", mode: 'copy'

    input:
    path label_lists_rds

    output:
    path "overlap_stats/", emit: results

    script:
    """
    overlap_stats_cli.R \
        --label_lists_rds ${label_lists_rds} \
        --output_dir overlap_stats
    """
}
