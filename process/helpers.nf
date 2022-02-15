process publishResults {
    input:
        path f

    output:
        path f

    cpus 1
    publishDir path: params.outputFolder, mode: 'copy', overwrite: 'true'

    when:
        params.publishResults

    script:
    """
        echo '${f}'
    """
}