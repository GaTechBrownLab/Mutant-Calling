process combine_blast {

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(table)

    output:
        tuple val(gene_ID), path("${gene_ID}/${gene_ID}_complete_combined_blast.txt"), optional: true, emit: complete_blast

        
    script:
    """
    mkdir -p ${gene_ID}
    cat ${table.join(' ')} > "${gene_ID}/${gene_ID}_complete_combined_blast.txt"

    """
}