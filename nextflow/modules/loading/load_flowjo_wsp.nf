process LOAD_FLOWJO_WSP {
    tag "load_wsp"
    label 'process_low'

    input:
    path wsp_files    // staged WSP files
    val wsp_names

    output:
    path "pop_strings.rds", emit: pop_strings

    script:
    // Build comma-separated paths from staged files
    def staged_paths = wsp_files.collect { it.name }.join(',')
    """
    load_flowjo_wsp_cli.R \
        --wsp_paths '${staged_paths}' \
        --wsp_names '${wsp_names}'
    """
}
