process graphing_r {
    label 'r_graphing'
    publishDir "${params.outdir}/graphing_results/${gene_ID}", mode: 'copy'

    input:
        tuple val(gene_ID), path(functions_incomplete), path(functions_pres_abs_incomplete)

    output:
        path "*"

    script:
    """

    Cheat_code_pruned_lineage_loop_no_BiSSE.R \
        $params.tree \
        $params.tree_root \
        $gene_ID \
        $functions_incomplete \
        $functions_pres_abs_incomplete \
        $params.fastani \
        $params.env_data \
        $params.lineage_probability
    
    overall_graphs_clean.R $gene_ID "${gene_ID}_env.csv" "${gene_ID}_lin_env.csv"

    """
}

