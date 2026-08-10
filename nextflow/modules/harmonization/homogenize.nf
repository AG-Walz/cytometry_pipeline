process HOMOGENIZE_PANELS {
    tag "homogenize"
    label 'process_medium'

    publishDir "${params.output_dir}/harmonization", mode: 'copy',
               enabled: params.save_intermediates ?: false

    input:
    path rds_files    // collected per-sample RDS files

    output:
    path "homogenized/*.rds", emit: flowframes

    script:
    def replacing = params.markername_replacing ?: '{}'
    """
    mkdir -p input_rds
    cp *.rds input_rds/
    n_files=\$(ls input_rds/*.rds 2>/dev/null | wc -l)
    if [ "\$n_files" -eq 0 ]; then
        echo "ERROR: No RDS files found for panel homogenization" >&2
        exit 1
    fi

    homogenize_cli.R \
        --input_rds_dir input_rds \
        --markers_of_interest '${params.markers_of_interest}' \
        --markername_replacing '${replacing}' \
        --output_dir homogenized
    """
}
