process CONCAT_AND_REFORMAT {

    tag "CONCAT_AND_REFORMAT"

    input:
        path fasta_files

    output:
        path "all_host_genomes.fna", emit: all_hosts

    script:
    """
    # Concatenate fasta files
    cat ${fasta_files.join(' ')} > all_host_genomes.fna

    # Reformat fasta headers
    sed -i '/>/ s/ .*\$//' all_host_genomes.fna
"""
}