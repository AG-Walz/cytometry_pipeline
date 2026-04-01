process DIMRED {
    tag "dimred_${params.dimred_method}"
    label 'process_high'

    publishDir "${params.output_dir}/visualization", mode: 'copy'

    input:
    path expr_data_rds
    path cluster_ids_rds

    output:
    path "dimred/", emit: plots

    script:
    """
    dimred_cli.R \
        --expr_rds ${expr_data_rds} \
        --cluster_rds ${cluster_ids_rds} \
        --output_dir dimred \
        --method ${params.dimred_method} \
        --n_neighbors ${params.n_neighbors} \
        --umap_min_dist ${params.umap_min_dist} \
        --tsne_perplexity ${params.tsne_perplexity} \
        --tsne_max_iter ${params.tsne_max_iter} \
        --output_formats ${params.output_formats} \
        --plot_width ${params.plot_width} \
        --plot_height ${params.plot_height} \
        --seed ${params.seed}
    """
}
