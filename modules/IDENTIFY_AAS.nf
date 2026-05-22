process IDENTIFY_AAS {

    tag "IDENTIFY_AAS: ${gene_ID}-${genome}"
    
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}/alignments/csv", mode: 'copy'

    input:
    tuple val(gene_ID), path(file), val(genome)

    output:
    tuple val(gene_ID), path("${genome}.csv")

    script:
    """
    # Identify amino acid subsitutions, insertions, and deletions with identify_AAS_deletions.py
    identify_AAS_deletions.py "$gene_ID" "$genome" "$file"
    
    """
}