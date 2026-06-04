process BLAST_DB_P {
    
    tag "BLAST_DB_P: ${ref}"

    input:
        tuple val(ID), path(ref)

    output:
        tuple val(ID), path("${ref}_blast_db*"), emit: blast_db_path_prot

    script:
    """
    # Make blast database of each protein
    makeblastdb -in $ref -parse_seqids -out ${ref}_blast_db/${ref}_blast_db -dbtype prot

    """
}

