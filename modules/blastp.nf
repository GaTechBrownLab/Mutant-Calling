process blastp {
    label 'blast'
    tag "${gene_ID}"

    memory '20 GB'

    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}", mode: 'copy'

    input:
        tuple val(gene_ID), path(query), path(blast_db_path)

    output:
        tuple val(gene_ID), path("${gene_ID}_prot_blast.txt"), emit: blast_output

    script:
    """  
    blastp \
        -query $query \
        -task blastp \
        -db  "${blast_db_path}/${blast_db_path.name}" \
        -outfmt "6 qseqid sseqid sacc pident nident qlen length evalue slen qstart qend sstart send" \
        -evalue 0.01 \
        -num_threads ${task.cpus} \
        > "${gene_ID}_prot_blast.txt"
    
    """
}

