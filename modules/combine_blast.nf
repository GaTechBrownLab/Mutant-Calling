process combine_blast {

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(table)

    output:
        tuple val(gene_ID), path("${gene_ID}/${gene_ID}_complete_combined_blast.txt"), optional: true, emit: complete_blast

    script:
    """
    mkdir -p ${gene_ID}

    for f in ${table.join(' ')}; do
        awk -v fname=\$(basename \$f) -F '\\t' 'BEGIN {OFS="\\t"} {print fname, \$0}' \$f
    done > "${gene_ID}/${gene_ID}_complete_combined_blast.txt"

    #cat ${table.join(' ')} > "${gene_ID}/${gene_ID}_complete_combined_blast.txt"

    """
}