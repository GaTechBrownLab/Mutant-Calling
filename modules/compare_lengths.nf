process compare_lengths {

    publishDir "${params.outdir}/mutant_calling_output", pattern: "**/Blast_gene_status/**",mode: 'copy'

    input:
        tuple val(gene_ID), path(table), path(gene_file)

    output:
        path "**"
        tuple val(gene_ID), path("${gene_ID}/Blast_gene_status/Complete/**"), optional: true, emit: Complete
        tuple val(gene_ID), path("${gene_ID}/Blast_gene_status/Incomplete/**"), optional: true, emit: Incomplete
        tuple val(gene_ID), path("${gene_ID}/Blast_gene_status/Missing/**"), optional: true, emit: Missing
        tuple val(gene_ID), path("${gene_ID}/Complete_output.txt"), optional: true, emit: Complete_list
        tuple val(gene_ID), path("${gene_ID}/Incomplete_output.txt"), optional: true, emit: Incomplete_list
        tuple val(gene_ID), path("${gene_ID}/Missing_output.txt"), optional: true, emit: Missing_list

    script:
    """
    mkdir -p ${gene_ID}/Blast_gene_status/Complete ${gene_ID}/Blast_gene_status/Incomplete ${gene_ID}/Blast_gene_status/Missing

    length=\$(sed -n '2p' ${gene_file} | tr -d '\n' | wc -c)

    compare_lengths.py "\$length" "$table" "${gene_ID}" "$params.gene_proportion" "$params.gene_difference" "$params.tolerance"

"""
}