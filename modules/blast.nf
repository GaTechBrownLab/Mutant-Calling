process blast {

    conda './envs/mutant_calling_genomics.yaml'

    cpus "${params.blast_threads}"

    memory '20 GB'

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(query), val(genome_ID), path(blast_db_path)

    output:
        tuple val(gene_ID), path("${gene_ID}/${gene_ID}_gene_blast_results/${genome_ID}_${gene_ID}.txt"), emit: blast_output

    script:
    """  
    mkdir -p "${gene_ID}"
    mkdir -p "${gene_ID}/${gene_ID}_gene_blast_results"
    blastn -query $query -task blastn -db "${blast_db_path}/${blast_db_path}" -outfmt "6 qseqid sseqid sacc pident nident qlen length evalue slen qstart qend sstart send sstrand" -evalue 0.01 -num_threads $params.blast_threads >> "${gene_ID}/${gene_ID}_gene_blast_results/${genome_ID}_${gene_ID}.txt"
    """
}

