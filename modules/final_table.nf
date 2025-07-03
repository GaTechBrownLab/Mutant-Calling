process final_table {

    conda './envs/python_scripts.yaml'

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
    tuple val(gene_ID), path(csv), path(cd_hit_table), path(scripts)

    output:
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_final_mut_table.csv"), emit: final_mut_table
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_final_mut_table_refs.csv"), emit: final_mut_table_refs

    script:
    """
    mkdir -p ${gene_ID}
    python "${scripts}/final_table.py" "$gene_ID" "$csv" "$cd_hit_table"
    """
}