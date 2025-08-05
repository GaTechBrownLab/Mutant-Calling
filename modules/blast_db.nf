process blast_db {
    
    input:
        tuple val(ID), path(ref)

    output:
        tuple val(ID), path("${ref}_blast_db*"), emit: blast_db_path

    script:
    """
    makeblastdb -in $ref -parse_seqids -out ${ref}_blast_db/${ref}_blast_db -dbtype nucl

    """
}

