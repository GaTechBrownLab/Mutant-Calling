process MAKE_LINKING_FILE {

    tag "MAKE_LINKING_FILE: ${gene_ID}"
    
    input:
    tuple val(gene_ID), path(fasta_file)

    output:
    path("${gene_ID}_linking.tsv")

    script:
    """
    # Make file linking protein accessions to filenames
    grep '^>' ${fasta_file} | sed 's/^>//; s/ .*//' | awk -v id=${gene_ID} '{print id "\\t" \$0}' > ${gene_ID}_linking.tsv
    """
}