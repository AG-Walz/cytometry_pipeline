process BARPLOT {
    tag "barplot"
    label 'process_low'

    publishDir "${params.output_dir}/visualization", mode: 'copy'

    input:
    path split_result_rds

    output:
    path "barplots/", emit: plots

    script:
    """
    barplot_cli.R \
        --split_result_rds ${split_result_rds} \
        --output_dir barplots
    """
}
