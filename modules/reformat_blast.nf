process reformat_blast {

    conda './envs/python_scripts.yaml'

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
    tuple val(gene_ID), path(blast_table), path(scripts), val(WT_cutoff), val(mut_cutoff)

    output:
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_muts_and_sig_diff.csv"), emit: muts_and_sig_diff
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_WT.csv"), emit: WT
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_muts_accessions.txt"), emit: muts_accessions
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_potential_deletions_accessions.txt"), emit: potential_deletions_accessions
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_muts_graphing.csv"), emit: muts_graphing
    

    script:
    """
    mkdir -p ${gene_ID}
    python "${scripts}/reformat_blast_new.py" "$blast_table" "$gene_ID" "$WT_cutoff" "$mut_cutoff"
    """
}