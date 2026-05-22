process ANI_CLUSTER_FUNCTION {
    
    tag "ANI_CLUSTER_FUNCTION: ${gene_ID}"
    
    publishDir "${params.outdir}/ani_cluster_function/${gene_ID}", mode: 'copy'

    input:
    tuple val(gene_ID), path(functions), path(ani_clusters)

    output:
    tuple val(gene_ID), path("${gene_ID}_cluster_status_test.csv"), emit: cluster_status_test
    tuple val(gene_ID), path("${gene_ID}_ANI_cluster_env_counts_status.csv"), emit: final_mut_table
    tuple val(gene_ID), path("${gene_ID}_ANI_distributions.svg"), emit: ani_distributions
    tuple val(gene_ID), path("${gene_ID}_ANI_distributions_prop.svg"), emit: ani_distributions_prop
    tuple val(gene_ID), path("${gene_ID}_ANI_clusters_chi.csv"), emit: ani_chi
    tuple val(gene_ID), path("${gene_ID}_environment_heatmap.svg"), emit: env_heatmap
    tuple val(gene_ID), path("${gene_ID}_environment_heatmap_funct.svg"), emit: env_heatmap_funct
    tuple val(gene_ID), path("${gene_ID}_no_funct_clusters_env.csv"), emit: no_funct_clusters_env
    
    script:
    """
    # Graph ANI clusters with associated gene status (function versus non-functional)
    ani_cluster_function.py "$gene_ID" "$functions" "$ani_clusters" "$params.env_data"

    ani_cluster_graphs.R "$gene_ID" "${gene_ID}_ANI_cluster_env_counts_status.csv" "$ani_clusters" "$params.env_data" "$functions"

    """
}