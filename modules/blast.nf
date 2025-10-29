process blast {
    
    cpus "${params.blast_threads}"

    memory '20 GB'

    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}/${gene_ID}_gene_blast_results", mode: 'copy'

    input:
        tuple val(gene_ID), path(query), val(genome_ID), path(blast_db_path)

    output:
        tuple val(gene_ID), path("${genome_ID}_${gene_ID}.txt"), emit: blast_output

    script:
    """  
    blastn \
        -query $query \
        -task blastn \
        -db "${blast_db_path}/${blast_db_path}" \
        -outfmt "6 qseqid sseqid sacc pident nident qlen length evalue slen qstart qend sstart send sstrand" \
        -evalue 0.01 \
        -num_threads $params.blast_threads \
        >> "${genome_ID}_${gene_ID}.txt"

    """
}

