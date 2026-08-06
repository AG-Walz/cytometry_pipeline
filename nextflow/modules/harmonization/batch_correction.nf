process BATCH_CORRECT {
    tag "batch_${params.batch_method}"
    label 'process_high'

    publishDir "${params.output_dir}/harmonization/batch_corrected", mode: 'copy',
               enabled: params.save_intermediates ?: false

    when:
    params.batch_method != 'none'

    input:
    path rds_files    // collected per-sample RDS files

    output:
    path "corrected/*.rds", emit: flowframes

    script:
    def batch_names_arg = params.batch_names ? "--batch_names '${params.batch_names}'" : ""
    """
    mkdir -p input_rds
    cp *.rds input_rds/
    n_files=\$(ls input_rds/*.rds 2>/dev/null | wc -l)
    if [ "\$n_files" -eq 0 ]; then
        echo "ERROR: No RDS files found for batch correction" >&2
        exit 1
    fi

    batch_correct_cli.R \
        --input_rds_dir input_rds \
        --method ${params.batch_method} \
        ${batch_names_arg} \
        --markers_of_interest '${params.markers_of_interest}' \
        --batch_do_pca ${params.batch_do_pca} \
        --combat_method ${params.combat_method} \
        --cytonorm_n_clusters ${params.cytonorm_n_clusters} \
        --cytonorm_goal ${params.cytonorm_goal ?: 'mean'} \
        --seed ${params.seed} \
        --output_dir corrected
    """
}
