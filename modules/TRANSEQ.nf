process TRANSEQ {
    tag "${gene_ID}"
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}", mode: 'copy'

    input:
        tuple val(gene_ID), path(fasta)

    output:
        tuple val(gene_ID), path("${gene_ID}_protein_seq.faa"), emit: translated_prot

    script:
    """  
    # Translate ORF into protein
    transeq $fasta ${gene_ID}_protein_seq.faa -frame=1 -trim=TRUE
    
    """
}

