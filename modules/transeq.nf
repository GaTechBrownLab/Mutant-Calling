process transeq {

    conda './envs/mutant_calling_genomics.yaml'

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(fasta)

    output:
        tuple val(gene_ID), path("${gene_ID}/${gene_ID}_protein_seq.faa"), emit: translated_prot

    script:
    """  
    mkdir -p ${gene_ID}

    transeq $fasta ${gene_ID}/${gene_ID}_protein_seq.faa -frame=1 -trim=TRUE
    """
}

