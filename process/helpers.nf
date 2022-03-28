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

process dot2svg {

    input:
        val f
        val temp

    output:
        path 'dag.svg'

    cpus 1
    container params.dagContainer

    when:
        params.dag

    script:
    """
        dot ${f} -Tsvg -o dag.svg
    """


}