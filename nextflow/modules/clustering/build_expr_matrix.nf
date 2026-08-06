process BUILD_EXPR_MATRIX {
    tag "build_matrix"
    label 'process_high'

    publishDir "${params.output_dir}/clustering", mode: 'copy',
               enabled: params.save_intermediates ?: false

    input:
    path rds_files    // collected per-sample RDS files

    output:
    path "expr_data.rds", emit: expr_data

    script:
    """
    mkdir -p input_rds
    cp *.rds input_rds/
    n_files=\$(ls input_rds/*.rds 2>/dev/null | wc -l)
    if [ "\$n_files" -eq 0 ]; then
        echo "ERROR: No RDS files found for expression matrix construction" >&2
        exit 1
    fi

    build_expr_matrix_cli.R \
        --input_rds_dir input_rds \
        --markers_of_interest '${params.markers_of_interest}' \
        --output_rds expr_data.rds
    """
}
