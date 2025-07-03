process concat_linking_files {

    publishDir "${params.outdir}/linking_file", mode: 'copy'

    input:
    path linking_files

    output:
    path("linking_file.tsv")

    script:
    """
    cat ${linking_files.join(' ')} > linking_file.tsv
    """
}