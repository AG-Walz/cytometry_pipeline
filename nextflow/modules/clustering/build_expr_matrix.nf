process BUILD_EXPR_MATRIX {
    tag "build_matrix"
    label 'process_high'

    input:
    path rds_files    // collected per-sample RDS files

    output:
    path "expr_data.rds", emit: expr_data

    script:
    """
    mkdir -p input_rds
    cp *.rds input_rds/ 2>/dev/null || true

    build_expr_matrix_cli.R \
        --input_rds_dir input_rds \
        --markers_of_interest ${params.markers_of_interest} \
        --output_rds expr_data.rds
    """
}
