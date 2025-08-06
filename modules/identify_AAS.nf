process identify_AAS {
    
    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
    tuple val(gene_ID), path(file), val(genome)

    output:
    tuple val(gene_ID), path("${gene_ID}/alignments/csv/${genome}.csv")

    script:
    """
    mkdir -p ${gene_ID}/alignments/csv

    identify_AAS_deletions.py "$gene_ID" "$genome" "$file"
    
    """
}