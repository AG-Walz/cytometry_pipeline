process PREPROCESS_SAMPLE {
    tag "${sample_name}"
    label 'process_medium'

    publishDir "${params.output_dir}/preprocessing", mode: 'copy',
               enabled: params.save_intermediates ?: false

    input:
    tuple val(sample_name), path(input_rds)

    output:
    tuple val(sample_name), path("${sample_name}_pp.rds"), emit: flowframe

    script:
    def markers_arg = params.markers_of_interest ? "--markers_of_interest '${params.markers_of_interest}'" : ""
    """
    preprocess_sample_cli.R \
        --input_rds ${input_rds} \
        --sample_name ${sample_name} \
        --output_rds ${sample_name}_pp.rds \
        --transform_type ${params.transform_type} \
        --qc_type ${params.qc_type} \
        --pp_order ${params.pp_order} \
        ${markers_arg} \
        --quant_norm ${params.quant_norm} \
        --scaling ${params.scaling}
    """
}
