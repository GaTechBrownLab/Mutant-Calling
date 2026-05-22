process REFORMAT_BED {

    tag "REFORMAT_BED: ${gene_ID}"
    
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}", mode: 'copy'

    input:
    tuple val(gene_ID), path(blast_table)

    output:
    tuple val(gene_ID), path("${gene_ID}_bed.bed"), emit: bed_file

    script:
    """
    # Reformat blast output for bedtools with reformat_bed.py
    reformat_bed.py "$blast_table" "$gene_ID"
    
    """
}