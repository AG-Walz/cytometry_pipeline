process UMAP_FACETS {
    tag "facets_${params.population}"
    label 'process_medium'

    publishDir "${params.output_dir}/visualization", mode: 'copy'

    input:
    path expr_data_rds
    path labels_rds
    path subpop_rds
    path coords_rds

    output:
    path "umap_facets/", emit: plots

    script:
    """
    umap_facets_cli.R \
        --expr_rds ${expr_data_rds} \
        --labels_rds ${labels_rds} \
        --subpop_rds ${subpop_rds} \
        --coords_rds ${coords_rds} \
        --output_dir umap_facets \
        --normalize ${params.normalize} \
        --point_size ${params.point_size} \
        --point_alpha ${params.point_alpha} \
        --plot_width ${params.plot_width} \
        --plot_height ${params.plot_height} \
        --output_formats ${params.output_formats} \
        --n_cores ${task.cpus}
    """
}
