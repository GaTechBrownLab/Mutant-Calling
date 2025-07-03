process blast_db_p {

    conda './envs/mutant_calling_genomics.yaml'

    input:
        tuple val(ID), path(ref)

    output:
        tuple val(ID), path("${ref}_blast_db*"), emit: blast_db_path_prot

    script:
    """
    makeblastdb -in $ref -parse_seqids -out ${ref}_blast_db/${ref}_blast_db -dbtype prot

    """
}

