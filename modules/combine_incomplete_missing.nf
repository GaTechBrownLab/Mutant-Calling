process combine_incomplete_missing {

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(missing), path(incomplete)

    output:
        tuple val(gene_ID), path("${gene_ID}/${gene_ID}_incomplete.txt"), emit: incomplete_final_list
        tuple val(gene_ID), path("${gene_ID}/alignments/${gene_ID}_combined_aln.csv"), emit: complete_blast

        
    script:
    """
    mkdir -p ${gene_ID}

    touch ${gene_ID}/${gene_ID}_incomplete.txt
    touch "${gene_ID}/${gene_ID}_missing.txt"
    
    cat ${incomplete.join(' ')} > "${gene_ID}/${gene_ID}_incomplete.txt"
    cat ${missing.join(' ')} > "${gene_ID}/${gene_ID}_missing.txt"
    """
}