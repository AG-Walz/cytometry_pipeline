process LOAD_CSV_EXPRESSION {
    tag "load_csv_expression"
    label 'process_low'

    publishDir "${params.output_dir}/loading", mode: 'copy'

    input:
    path csv_dir

    output:
    path "expr_data_list.rds", emit: expr_list

    script:
    """
    load_csv_expression_cli.R \
        --input_dir ${csv_dir} \
        --output_rds expr_data_list.rds
    """
}
