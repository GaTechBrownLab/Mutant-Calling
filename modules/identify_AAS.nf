process identify_AAS {
    
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}/alignments/csv", mode: 'copy'

    input:
    tuple val(gene_ID), path(file), val(genome)

    output:
    tuple val(gene_ID), path("${genome}.csv")

    script:
    """
    identify_AAS_deletions.py "$gene_ID" "$genome" "$file"
    
    """
}