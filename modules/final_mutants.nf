process final_mutants {
    
    publishDir "${params.outdir}/ouputs_for_graphing", mode: 'copy'

    input:
    tuple val(gene_ID), path(incomplete_files), path(missing_files), path(final_mut_table), path(linking), path(muts_graphing), path(input_muts)

    output:
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_incomplete.txt"), emit: incomplete
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_missing.txt"), emit: missing
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_functions_pres_abs_incomplete.csv"), emit: functions_pres_abs_incomplete
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_all_functions_incomplete.csv"), emit: all_functions_incomplete
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_functions_incomplete.csv"), emit: functions_incomplete

    script:
    """
    mkdir -p ${gene_ID}

    touch ${gene_ID}/${gene_ID}_incomplete.txt
    touch "${gene_ID}/${gene_ID}_missing.txt"
    
    cat ${incomplete_files.join(' ')} > "${gene_ID}/${gene_ID}_incomplete.txt"
    cat ${missing_files.join(' ')} > "${gene_ID}/${gene_ID}_missing.txt"

    reformat_muts_new_2.py "$final_mut_table" "$gene_ID" ""${gene_ID}/${gene_ID}_incomplete.txt"" "${gene_ID}/${gene_ID}_missing.txt" "$linking" "$muts_graphing" "$input_muts"

    echo "1"
    """
}