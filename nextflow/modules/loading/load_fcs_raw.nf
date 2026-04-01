process LOAD_FCS_RAW {
    tag "${sample_name}"
    label 'process_low'

    input:
    tuple val(sample_name), path(fcs_file)

    output:
    tuple val(sample_name), path("${sample_name}.rds"), emit: flowframe

    script:
    """
    load_fcs_raw_cli.R \
        --fcs_file ${fcs_file} \
        --sample_name ${sample_name} \
        --output_rds ${sample_name}.rds \
        --truncation ${params.truncation}
    """
}
