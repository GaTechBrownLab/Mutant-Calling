process split_files {

    conda './envs/python_scripts.yaml'

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(cd_hit_clust), val(scripts)

    output:
        tuple val(gene_ID), path("${gene_ID}/${gene_ID}_cd_hit_table.csv"), emit: cd_hit_table
        tuple val(gene_ID), path("${gene_ID}/${gene_ID}_muts.txt"), emit: muts_list
        
    script:
    """
    mkdir -p ${gene_ID}

    csplit -z ${cd_hit_clust} /Cluster/ '{*}'

    for infile in xx*;
    do
        sed -i '1d' \${infile}
    done

    python "${scripts}/reformat_cd_hit.py" "$gene_ID"

    # Reformat CD-Hit table
    sed -i 's/\\.\\.\\.//g' "${gene_ID}/${gene_ID}_cd_hit_table.csv"

    # Reformat CD-Hit references file
    sed -i 's/\\.\\.\\.\$//' "${gene_ID}/${gene_ID}_muts.txt"

    """
}