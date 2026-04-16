process LOAD_CSV_GATING {
    tag "load_csv_gating"
    label 'process_low'

    publishDir "${params.output_dir}/loading", mode: 'copy'

    input:
    path input_dir

    output:
    path "loaded_gating/", emit: gating_dir

    script:
    """
    load_csv_gating_cli.R \
        --input_fcs_path ${input_dir} \
        --input_csv_path ${input_dir} \
        --output_dir loaded_gating \
        --gating_pattern '${params.gating_pattern}' \
        --processing_mode ${params.processing_mode}
    """
}
