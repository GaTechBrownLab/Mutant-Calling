process BLASTP {
    
    tag "BLASTP: ${query}-${gene_ID}"
    
    cpus "${params.blast_threads}"

    memory '20 GB'

    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}", mode: 'copy'

    input:
        tuple val(gene_ID), path(query), path(blast_db_path)

    output:
        tuple val(gene_ID), path("${gene_ID}_prot_blast.txt"), emit: blast_output

    script:
    """  
    # Run blastp on all protein sequences
    blastp \
        -query $query \
        -task blastp \
        -db  "${blast_db_path}/${blast_db_path}" \
        -outfmt "6 qseqid sseqid sacc pident nident qlen length evalue slen qstart qend sstart send" \
        -evalue 0.01 \
        -num_threads $params.blast_threads \
        >> "${gene_ID}_prot_blast.txt"
    
    """
}

