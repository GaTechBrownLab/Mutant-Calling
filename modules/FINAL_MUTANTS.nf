process FINAL_MUTANTS {

    tag "FINAL_MUTANTS: ${gene_ID}"
    
    publishDir "${params.outdir}/ouputs_for_graphing/${gene_ID}", mode: 'copy'

    input:
    tuple val(gene_ID), path(incomplete_files), path(missing_files), path(split_files), path(final_mut_table), path(linking), path(muts_graphing), path(input_muts)

    output:
    tuple val(gene_ID), path("${gene_ID}_functions_pres_abs_incomplete.csv"), emit: functions_pres_abs_incomplete
    tuple val(gene_ID), path("${gene_ID}_all_functions_incomplete.csv"), emit: all_functions_incomplete
    tuple val(gene_ID), path("${gene_ID}_functions_incomplete.csv"), emit: functions_incomplete

    script:
    """
    # For identification of mutants using specific target amino acid subsitutions
    reformat_muts_new_2.py \
        "$final_mut_table" \
        "$gene_ID" \
        "$incomplete_files" \
        "$missing_files" \
        "$split_files" \
        "$linking" \
        "$muts_graphing" \
        "$input_muts"

    """
}