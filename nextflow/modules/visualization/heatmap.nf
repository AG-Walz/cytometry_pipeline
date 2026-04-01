process HEATMAP {
    tag "heatmap"
    label 'process_low'

    publishDir "${params.output_dir}/visualization", mode: 'copy'

    input:
    path expr_rds

    output:
    path "heatmap.png", emit: plot

    script:
    """
    heatmap_cli.R \
        --expr_rds ${expr_rds} \
        --output_path heatmap.png \
        --max_cells ${params.max_cells_per_file ?: 500}
    """
}
