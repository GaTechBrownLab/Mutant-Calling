process final_table {
    
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}/final_mut_tables", mode: 'copy'

    input:
    tuple val(gene_ID), path(csv), path(cd_hit_table)

    output:
    tuple val(gene_ID), path("${gene_ID}_final_mut_table.csv"), emit: final_mut_table
    tuple val(gene_ID), path("${gene_ID}_final_mut_table_refs.csv"), emit: final_mut_table_refs

    script:
    """
    final_table.py "$gene_ID" "$csv" "$cd_hit_table"
    
    """
}