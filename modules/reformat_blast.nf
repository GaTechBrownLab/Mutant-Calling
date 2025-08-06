process reformat_blast {
    
    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
    tuple val(gene_ID), path(blast_table)

    output:
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_muts_and_sig_diff.csv"), emit: muts_and_sig_diff
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_WT.csv"), emit: WT
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_muts_accessions.txt"), emit: muts_accessions
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_potential_deletions_accessions.txt"), emit: potential_deletions_accessions
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_muts_graphing.csv"), emit: muts_graphing
    

    script:
    """
    mkdir -p ${gene_ID}
    
    reformat_blast_new.py "$blast_table" "$gene_ID" "$params.WT_cutoff" "$params.mut_cutoff"

    """
}