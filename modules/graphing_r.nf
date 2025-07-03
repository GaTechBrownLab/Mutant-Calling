process graphing_r {

    conda './envs/graphing_r.yaml'

    publishDir 'results/graphing_results', mode: 'copy'

    input:
        tuple val(gene_ID), path(functions_incomplete), path(functions_pres_abs_incomplete), path(tree), val(tree_root), path(fastani), path(env_data), val(lineage_probability), path(scripts)

    output:
        path "*"

    script:
    """
    mkdir -p $gene_ID

    Rscript ${scripts}/Cheat_code_pruned_lineage_loop_no_BiSSE.R $tree $tree_root $gene_ID $functions_incomplete $functions_pres_abs_incomplete $fastani $env_data $lineage_probability

    Rscript ${scripts}/overall_graphs_clean.R $gene_ID "${gene_ID}/${gene_ID}_env.csv" "${gene_ID}/${gene_ID}_lin_env.csv"

    """
}

