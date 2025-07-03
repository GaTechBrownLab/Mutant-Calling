process concat_and_reformat {
    input:
        path fasta_files

    output:
        path "all_host_genomes.fna", emit: all_hosts

    script:
    """
    cat ${fasta_files.join(' ')} > all_host_genomes.fna
    sed -i '/>/ s/ .*\$//' all_host_genomes.fna
"""
}