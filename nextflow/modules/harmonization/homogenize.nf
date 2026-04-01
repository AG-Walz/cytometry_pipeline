process HOMOGENIZE_PANELS {
    tag "homogenize"
    label 'process_medium'

    input:
    path rds_files    // collected per-sample RDS files

    output:
    path "homogenized/*.rds", emit: flowframes

    script:
    def replacing = params.markername_replacing ?: '{}'
    """
    mkdir -p input_rds
    cp *.rds input_rds/ 2>/dev/null || true

    homogenize_cli.R \
        --input_rds_dir input_rds \
        --markers_of_interest ${params.markers_of_interest} \
        --markername_replacing '${replacing}' \
        --output_dir homogenized
    """
}
