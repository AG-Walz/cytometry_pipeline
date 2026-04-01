process COMPARE_GATINGS {
    tag "compare_gatings"
    label 'process_medium'

    publishDir "${params.output_dir}/gating_comparison", mode: 'copy'

    input:
    path pop_strings_rds

    output:
    path "gating_comparison/",                    emit: plots
    path "gating_comparison/label_lists.rds",     emit: label_lists

    script:
    def label_dict_str = params.label_dict ?: '{}'
    def names_arg = params.wsp_names ? "--custom_name '${params.wsp_names}'" : ""
    def target_arg = params.target_fcs ? "--target_fcs ${params.target_fcs}" : ""
    """
    compare_gatings_cli.R \
        --pop_strings_rds ${pop_strings_rds} \
        ${target_arg} \
        --label_dict '${label_dict_str}' \
        ${names_arg} \
        --output_dir gating_comparison
    """
}
