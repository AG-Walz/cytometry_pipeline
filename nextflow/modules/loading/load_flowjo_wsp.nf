process LOAD_FLOWJO_WSP {
    tag "load_wsp"
    label 'process_low'

    input:
    val wsp_paths
    val wsp_names

    output:
    path "pop_strings.rds", emit: pop_strings

    script:
    """
    load_flowjo_wsp_cli.R \
        --wsp_paths '${wsp_paths}' \
        --wsp_names '${wsp_names}'
    """
}
