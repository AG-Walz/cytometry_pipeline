process DOWNSAMPLE {
    tag "${sample_name}"
    label 'process_low'

    input:
    tuple val(sample_name), path(input_rds)

    output:
    tuple val(sample_name), path("${sample_name}_ds.rds"), emit: flowframe

    script:
    """
    downsample_cli.R \
        --input_rds ${input_rds} \
        --sample_name ${sample_name} \
        --output_rds ${sample_name}_ds.rds \
        --target_cells ${params.downsample_n} \
        --method ${params.downsample_method} \
        --seed ${params.seed}
    """
}
