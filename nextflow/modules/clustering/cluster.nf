process CLUSTER {
    tag "cluster_${params.clustering_method}"
    label 'process_high'

    publishDir "${params.output_dir}/clustering", mode: 'copy'

    input:
    path expr_data_rds

    output:
    path "cluster_ids.rds", emit: cluster_ids

    script:
    def eps_arg = params.dbscan_eps ? "--dbscan_eps ${params.dbscan_eps}" : ""
    """
    cluster_cli.R \
        --input_rds ${expr_data_rds} \
        --output_rds cluster_ids.rds \
        --method ${params.clustering_method} \
        --n_neighbors ${params.n_neighbors} \
        --resolution ${params.resolution} \
        --n_clusters ${params.n_clusters} \
        --grid_xdim ${params.grid_xdim} \
        --grid_ydim ${params.grid_ydim} \
        --dbscan_min_pts ${params.dbscan_min_pts} \
        ${eps_arg} \
        --seed ${params.seed}
    """
}
