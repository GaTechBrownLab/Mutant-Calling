process combine_aln {

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(table)

    output:
        tuple val(gene_ID), path("${gene_ID}/alignments/${gene_ID}_combined_aln.csv"), emit: complete_blast

        
    script:
    """
    mkdir -p ${gene_ID}/alignments
    cat ${table.join(' ')} > "${gene_ID}/alignments/${gene_ID}_combined_aln.csv"

    """
}