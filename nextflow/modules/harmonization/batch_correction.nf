process BATCH_CORRECT {
    tag "batch_${params.batch_method}"
    label 'process_high'

    when:
    params.batch_method != 'none'

    input:
    path rds_files    // collected per-sample RDS files

    output:
    path "corrected/*.rds", emit: flowframes

    script:
    """
    mkdir -p input_rds
    cp *.rds input_rds/ 2>/dev/null || true

    batch_correct_cli.R \
        --input_rds_dir input_rds \
        --method ${params.batch_method} \
        --batch_names '${params.batch_names}' \
        --markers_of_interest ${params.markers_of_interest} \
        --batch_do_pca ${params.batch_do_pca} \
        --combat_method ${params.combat_method} \
        --cytonorm_n_clusters ${params.cytonorm_n_clusters} \
        --cytonorm_goal ${params.cytonorm_goal ?: 'mean'} \
        --seed ${params.seed} \
        --output_dir corrected
    """
}
